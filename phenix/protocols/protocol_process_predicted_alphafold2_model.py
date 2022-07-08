# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

from pyworkflow import Config
from pyworkflow import utils as pwutils
from pwem.protocols import EMProtocol
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, BooleanParam, EnumParam,
                                        StringParam, FloatParam, IntParam, PathParam)
from phenix.constants import PROCESS, PHENIX_HOME
from pwem.convert.atom_struct import fromCIFToPDB, fromPDBToCIF, \
    fromCIFTommCIF, AtomicStructHandler, retry

try:
    from pwem.objects import AtomStruct, Sequence
except:
    from pwem.objects import PdbFile as AtomStruct
from phenix import Plugin


class PhenixProtProcessPredictedAlphaFold2Model(EMProtocol):
    """Process Predicted Model
     Replace values in b-factor field with estimated B values.
     Optionally remove low-confidence residues and split into domains."""
    _label = 'process predicted model'
    _program = ""
    # _version = VERSION_1_2
    PROCESSPREDICTEDFILE = '_processed.pdb'
    SEQREMAINDER = '_remainder.seq'
    ContentOfBvalueField = ['LDDT (AlphaFold2)', 'RMSD', 'B-value']

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPredictedModel', PointerParam,
                      pointerClass="AtomStruct",
                      label='Predicted AlphaFold2 model', important=True,
                      help="Atom structure model (PDBx/mmCIF) retrieved from AlphaFold2.")
        form.addParam('contentBvalueField', EnumParam,
                      choices=self.ContentOfBvalueField, important=True,
                      label="Contents of B-value field:", default=0,
                      help="The B-value field in most predicted models represents \n"
                      "one of three possible values:\nA confidence (LDDT) on a scale "
                      "of either 0 to 1 or 0 to 100.\nAn estimate of error in A (rmsd)\n"
                      "An actual B-value (atomic displacement parameter)\n"
                      "In process_predicted_model, confidence values or error estimates "
                      "in A or are first converted in new pseudo B-values.")
        form.addParam('minLDDT', IntParam, default=70,
                      label='Minimun LDDT value',
                      condition=('contentBvalueField==%d ' % 0),
                      help="""Cutoff value to remove low-confidence residues. 
                       Values of LDDT range between 0 and 100. A minimum LDDT 
                       of 70 corresponds to a maximum RMSD of 1.5.\nModel 
                       Confidence:\nVery high (pLDDT > 90)\n
                       Confident (90 > pLDDT > 70)\nLow (70 > pLDDT > 50)\n
                       Very low (pLDDT < 50): Probably unstructured in isolation""")
        form.addParam('maxRMSD', FloatParam, default=1.5,
                      label='Maximum RMSD value',
                      condition=('contentBvalueField==%d ' % 1),
                      help="""Cutoff value to remove low-confidence residues.\n 
                      A maximum RMSD of 1.5 A corresponds to a minimum LDDT 
                      of 70.\nModel Confidence:\nVery high (pLDDT > 90)\n
                      Confident (90 > pLDDT > 70)\nLow (70 > pLDDT > 50)\n
                      Very low (pLDDT < 50): Probably unstructured in isolation""")
        form.addParam('paeFile', PointerParam,
                      pointerClass="PAE", allowsNull=True,
                      label='PAE file', condition=('contentBvalueField!=%d ' % 2),
                      help="Optional input .json file with matrix of inter-residue"
                           " estimated errors.")
        form.addParam('removeLowConfidenceResidues', BooleanParam, default=True,
                      label='Remove low-confidence residues',
                      help="""For AlphaFold2 models, low-confidence corresponds 
                      approximately to an LDDT value of about 0.7 (on a scale of 
                      0 to 1, or 70 on a scale of 0 to 100), or to an RMSD value 
                      of about 1.5, or to a B-value of about 60.""")
        form.addParam('splitModel', BooleanParam, default=True,
                      label='Split model into compact regions',
                      help="""group the pieces from your trimmed model into compact 
                      domains, or even to split some pieces into compact domains""")
        form.addParam('maximumDomains', IntParam, default=3,
                      label='Processing option: Maximum domains',
                      help="""Maximum domains to obtain. You can use this to merge
                      the closest domains at the end of splitting the model. Make it
                      bigger to get more domains.""")
        form.addParam('minimumDomainLength', IntParam, default=10,
                      label='Processing option: Minimum domain length (residues)',
                      help="""Minimum length of a domain to keep (reject at end if
                      smaller).""")
        form.addParam('extraParams', StringParam,
                      label="Extra Params ",
                      default="",
                      expertLevel=LEVEL_ADVANCED,
                      help="This string will be added to the phenix command.\n"
                           "Syntax: paramName1=value1 paramName2=value2 ")

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):

        self._insertFunctionStep('runProcessPredictedModel')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def runProcessPredictedModel(self):
        atomStruct = os.path.abspath(
            self.inputPredictedModel.get().getFileName())

        atomStruct_localPath = os.path.abspath(
            self._getExtraPath(atomStruct.split('/')[-1]))

        if str(atomStruct) != str(atomStruct_localPath):
            pwutils.path.createLink(atomStruct, atomStruct_localPath)
            atomStruct = atomStruct_localPath
        args = self._writeArgsProcessAlphaFold(atomStruct)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(PROCESS),
              args, cwd=cwd,
              listAtomStruct=[atomStruct], log=self._log)
    def createOutputStep(self):
        pdb = AtomStruct()
        for fileName in os.listdir(self._getExtraPath()):
            if fileName.endswith(self.PROCESSPREDICTEDFILE):
                pdb.setFileName(self._getExtraPath(fileName))
            else:
                print("fileName: ", fileName)
        self.output = pdb.getFileName().split('/')[-1]
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputPredictedModel.get(), pdb)

        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram(PROCESS)
        if not os.path.exists(program):
            errors.append("Cannot find " + program)

            # If there is any error at this point it is related to config variables
            errors.append("Check configuration file: " +
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set PHENIX_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % Plugin.getVar(PHENIX_HOME))
                errors.append("PROCESS = %s" % PROCESS)

        # Check if the B-factor column contains common LDDT, RMSD or B-factor
        # values and avoid running the protocol if the average of those values
        # don't follow the expected value
        atomStruct = os.path.abspath(
            self.inputPredictedModel.get().getFileName())
        self.structureHandler = AtomicStructHandler()
        self.structureHandler.read(atomStruct)
        bFactorValues, listOfResiduesBFactors = \
            self.structureHandler. getStructureBFactorValues()
        if self._average(bFactorValues) > 20.0 and self.contentBvalueField == 1:
            errors.append("WARNING!!!: average B-factor column > 20.0\n"
                          "Review your input prediction file\n"
                          "(check the values of the B-factor column)")

        # Check  if there are at least 5 sequential residues satisfying the threshold
        if self.contentBvalueField == 0:
            threshold = self.minLDDT
        elif self.contentBvalueField == 1:
            threshold = self.maxRMSD
        result = \
            self._minSequentialResidues(
                listOfResiduesBFactors, self.contentBvalueField, threshold)
        if not result:
            errors.append("WARNING!!!:\n"
                          "Less than five sequential residues matching params")
        return errors

    def _summary(self):
        summary = []
        # try:
        #     summary.append("protocol finished with results")
        # except:
        #     summary.append("processed predicted model not yet docked")
        summary.append(
            "https://phenix-online.org/version_docs/dev-4380/reference/process_predicted_model.html")
        # summary.append("Tom Terwilliger, Claudia Millan Nebot, Tristan Croll")

        return summary

    def _citations(self):
        return ['Terwilliger_2022']

    # --------------------------- UTILS functions --------------------------

    def _writeArgsProcessAlphaFold(self, atomStruct):
        args = " "
        args += "%s " % atomStruct
        args += " "
        args += "b_value_field_is="
        if self.contentBvalueField == 0:
            args += "lddt"
            if self.minLDDT != 70:
                args += " minimum_lddt=%d maximum_rmsd=None  " % self.minLDDT
            if self.paeFile.get():
                args += " pae_file=%s  " % self.paeFile.get()
        elif self.contentBvalueField == 1:
            args += "rmsd"
            if self.maxRMSD != 1.5:
                args += " minimum_lddt=None maximum_rmsd=%.2f  " % self.maxRMSD
            if self.paeFile.get():
                args += " pae_file=%s  " % self.paeFile.get()
        else:
            args += "b_value"
        if self.removeLowConfidenceResidues != True:
            args += " remove_low_confidence_residues=False"
        if self.splitModel != True:
            args += " split_model_by_compact_regions=False"
        if self.maximumDomains != 3:
            args += " maximum_domains=" + str(self.maximumDomains)
        if self.minimumDomainLength != 10:
            args += " minimum_domain_length=" + str(self.minimumDomainLength)
        if len(str(self.extraParams)) > 0:
            args += " %s " % self.extraParams.get()
        return args

    def _average(self, list):
        sum_list = 0.00
        for item in list:
            sum_list += item
        avg = sum_list / len(list)
        return avg

    def _minSequentialResidues(self, listOfResiduesBFactors, contentBvalueField, threshold):
        """This method retuns True when there exist 5 consecutive residues matching the
        threshold condition (example: LDDT > 70.00 or RMSD < 1.50). It returns False if there
         aren't 5 consecutive residues matching the threshold condition"""
        list_res = []
        for item in listOfResiduesBFactors:
            # item example LDDT (Bfactor colum, residue average(atoms)) == 70.0:
            # [0, 'A', 103, 'PHE', 96.01]
            # adding 1 (TRUE) or 0 (FALSE) to each item list if it satisfies
            # the threshold condition
            # item example with LDDT (Bfactor colum, residue average(atoms))==
            # 70.0 [0, 'A', 103, 'PHE', 96.01, 1]
            if contentBvalueField == 0:
                if item[4] > threshold:
                    item.append(1)
                else:
                    item.append(0)
            elif contentBvalueField == 1:
                if item[4] < threshold:
                    item.append(1)
                else:
                    item.append(0)
            list_res.append(item)
        size = len(list_res)
        # minimum_sequential_residues = 5
        for i in range(size - 4):
            # same model, same chain and same threshold condition (TRUE)
            if list_res[i][0] == list_res[i + 1][0] == list_res[i + 2][0] \
                    == list_res[i + 3][0] == list_res[i + 4][0] and \
                    list_res[i][1] == list_res[i + 1][1] == list_res[i + 2][1] \
                    == list_res[i + 3][1] == list_res[i + 4][1] and \
                    list_res[i][1] == list_res[i + 1][1] == list_res[i + 2][1] \
                    == list_res[i + 3][1] == list_res[i + 4][1] and \
                    list_res[i][5] == list_res[i + 1][5] == list_res[i + 2][5] \
                    == list_res[i + 3][5] == list_res[i + 4][5] == 1:
                return True
        return False
