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
                                        StringParam, FloatParam, IntParam)
from phenix.constants import PROCESS, PHENIX_HOME
from pwem.convert.atom_struct import fromCIFToPDB, fromPDBToCIF, fromCIFTommCIF, AtomicStructHandler, retry

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
                      label='Minimun LDDT value (max. RMSD)',
                      help="""Cutoff value to remove low-confidence residues. 
                       Values of LDDT range between 0 and 100. A minimum LDDT 
                       of 70 corresponds to a maximum RMSD of 1.5.""")
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
        form.addParam('maximumDomains', FloatParam, default=3.0,
                      label='Processing option: Maximum domains',
                      help="""Maximum domains to obtain. You can use this to merge
                      the closest domains at the end of splitting the model. Make it
                      bigger to get more domains.""")
        form.addParam('minimumDomainLength', FloatParam, default=10.0,
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
        elif self.contentBvalueField == 1:
            args += "rmsd"
        else:
            args += "b_value"
        if self.minLDDT != 70:
            args += " minimum_lddt=%d maximum_rmsd=None  " % self.minLDDT
        if self.removeLowConfidenceResidues != True:
            args += " remove_low_confidence_residues=False"
        if self.splitModel != True:
            args += " split_model_by_compact_regions=False"
        if self.maximumDomains != 3.0:
            args += " maximum_domains=" + str(self.maximumDomains)
        if self.minimumDomainLength != 10.0:
            args += " minimum_domain_length=" + str(self.minimumDomainLength)
        if len(str(self.extraParams)) > 0:
            args += " %s " % self.extraParams.get()
        return args