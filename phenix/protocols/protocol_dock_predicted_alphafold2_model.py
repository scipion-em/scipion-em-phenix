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
from pwem.convert import Ccp4Header
from pwem.protocols import EMProtocol
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, BooleanParam, EnumParam,
                                        StringParam, FloatParam, IntParam)
from phenix.constants import DOCKPREDICTEDMODEL, PHENIX_HOME
from pwem.convert.atom_struct import fromCIFToPDB, fromPDBToCIF, \
    fromCIFTommCIF, AtomicStructHandler, retry

try:
    from pwem.objects import AtomStruct, Sequence
except:
    from pwem.objects import PdbFile as AtomStruct
from phenix import Plugin


class PhenixProtDockPredictedAlphaFold2Model(EMProtocol):
    """Dock Predicted Model.
     Dock predicted model docks the domains from a model produced by AlphaFold,
     RoseTTAFold and other prediction software into a cryo EM map. It uses the
     connectivity of the model as a restraint in the docking process so that
     the docked domains normally are in a reasonable arrangement. It can take
     map symmetry into account."""
    _label = 'dock predicted model'
    _program = ""
    # _version = VERSION_1_2
    DOCKINMAPFILE = 'docking_map.mrc'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPredictedModel', PointerParam,
                      pointerClass="AtomStruct",
                      label='Predicted AlphaFold2 model', important=True,
                      help="Atom structure model (PDBx/mmCIF) retrieved from AlphaFold2.")
        form.addParam('inputProcessedPredictedModel', PointerParam,
                      pointerClass="AtomStruct",
                      label='Processed AlphaFold2 model', important=True,
                      help="Atom structure model (PDBx/mmCIF) retrieved from AlphaFold2 "
                           "and processed by Phenix process_predicted_model protocol.")
        # form.addParam('modelCopies', IntParam,
        #               default=1,
        #               label='Number of copies',
        #               help="Write here the number of symmetric copies of your atom "
        #                    "structure. ")
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input map', important=True,
                      help="Set the starting density map.")
        # form.addParam('asymmetricMap', BooleanParam, default=True,
        #               label='Assymetric map:',
        #               help="If your map has symmetry be sure to set this param No."
        #                    " Otherwise symmetry will be automatically determined.")
        form.addParam('resolution', FloatParam, default=3.0,
                      label='High-resolution limit (A):',
                      help="Map resolution (Angstroms).")
        form.addParam('numberOfThreads', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1,
                      label='Number of threads',
                      help="Write here the number of threads to run the protocol. ")
        form.addParam('extraParams', StringParam,
                      label="Extra Params ",
                      default="",
                      expertLevel=LEVEL_ADVANCED,
                      help="This string will be added to the phenix command.\n"
                           "Syntax: paramName1=value1 paramName2=value2 ")

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runDockPredictedModel')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def convertInputStep(self):
        """ convert 3D maps to MRC '.mrc' format
        """
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        newFn = self._getExtraPath(self.DOCKINMAPFILE)
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        Ccp4Header.fixFile(inVolName, newFn, origin, sampling, Ccp4Header.START)  # ORIGIN

    def runDockPredictedModel(self):
        predictedAtomStruct = os.path.abspath(
            self.inputPredictedModel.get().getFileName())
        processedAtomStruct = os.path.abspath(
            self.inputProcessedPredictedModel.get().getFileName())
        mapFile = self.DOCKINMAPFILE
        vol = os.getcwd() + "/" + self._getExtraPath(mapFile)

        predictedAtomStruct_localPath = os.path.abspath(
            self._getExtraPath(predictedAtomStruct.split('/')[-1]))
        if str(predictedAtomStruct) != str(predictedAtomStruct_localPath):
            pwutils.path.createLink(predictedAtomStruct, predictedAtomStruct_localPath)
            predictedAtomStruct = predictedAtomStruct_localPath
        processedAtomStruct_localPath = os.path.abspath(
            self._getExtraPath(processedAtomStruct.split('/')[-1]))
        if str(processedAtomStruct) != str(processedAtomStruct_localPath):
            pwutils.path.createLink(processedAtomStruct, processedAtomStruct_localPath)
            processedAtomStruct = processedAtomStruct_localPath
        prefix = os.path.abspath(self._getExtraPath(processedAtomStruct))
        args = self._writeArgsDockAlphaFold(
            predictedAtomStruct, processedAtomStruct, vol, prefix)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(DOCKPREDICTEDMODEL),
              args, cwd=cwd,
              listAtomStruct=[predictedAtomStruct, processedAtomStruct],
              log=self._log)
    def createOutputStep(self):
        pdb = AtomStruct()
        for fileName in os.listdir(self._getExtraPath()):
            if (fileName.endswith(".cif.pdb") or fileName.endswith(".pdb.pdb")):
                pdb.setFileName(self._getExtraPath(fileName))
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputPredictedModel.get(), pdb)
        self._defineSourceRelation(self.inputProcessedPredictedModel.get(), pdb)

        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram(DOCKPREDICTEDMODEL)
        if not os.path.exists(program):
            errors.append("Cannot find " + program)

            # If there is any error at this point it is related to config variables
            errors.append("Check configuration file: " +
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set PHENIX_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % Plugin.getVar(PHENIX_HOME))
                errors.append("PROCESS = %s" % DOCKPREDICTEDMODEL)

        return errors

    def _summary(self):
        summary = []
        # try:
        #     summary.append("protocol finished with results")
        # except:
        #     summary.append("processed predicted model not yet docked")
        summary.append(
            "https://phenix-online.org/version_docs/dev-4380/reference/dock_predicted_model.html")

        return summary

    def _citations(self):
        return ['Terwilliger_2022']

    # --------------------------- UTILS functions --------------------------

    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputPredictedModel.get().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol

    def _writeArgsDockAlphaFold(
            self, predictedAtomStruct, processedAtomStruct, vol, prefix):
        args = " "
        args += "model=%s " % predictedAtomStruct
        args += "processed_model_file=%s " % processedAtomStruct
        # if self.modelCopies > 1:
        #     args += " search_model_copies=%d" % self.modelCopies
        #     args += " use_symmetry=True "
        args += "full_map=%s " % vol
        # if not self.asymmetricMap:
        #     args += " asymmetric_map=False "
        args += "resolution=%f" % self.resolution
        args += " "
        args += "output_model_prefix=%s" % prefix
        args += " "
        if self.numberOfThreads > 1:
            print("self.numberOfThreads: ", self.numberOfThreads)
            args += "nproc=%d " % self.numberOfThreads
        if len(str(self.extraParams)) > 0:
            args += " %s " % self.extraParams.get()
        return args

    def _runChangingCifFormatSuperpose(self, list_args):
        cwd = os.getcwd() + "/" + self._getExtraPath()
        try:
            if list_args[0].endswith(".cif") and list_args[1].endswith(".cif"):
                try:
                    # upgrade cifs
                    list_args1 = []
                    for i in range(0, 2):
                        list_args1.append(fromCIFTommCIF(list_args[i], list_args[i]))
                    args1 = list_args1[0] + " " + list_args1[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args1,
                                            extraEnvDict=None, cwd=cwd)
                except:
                    # convert cifs to pdbs
                    list_args2 = []
                    for i in range(0, 2):
                        list_args2.append(fromCIFToPDB(
                            list_args[i], list_args[i].replace('.cif', '.pdb')))
                    args2 = list_args2[0] + " " + list_args2[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args2,
                                            extraEnvDict=None, cwd=cwd)
            elif list_args[0].endswith(".cif") and list_args[1].endswith(".pdb"):
                try:
                    # pdbs: convert cif to pdb
                    list_args1 = []
                    list_args1.append(fromCIFToPDB(
                        list_args[0], list_args[0].replace('.cif', '.pdb')))
                    args1 = list_args1[0] + " " + list_args[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args1,
                                            extraEnvDict=None, cwd=cwd)
                except:
                    try:
                        # cifs: convert pdb to cif
                        list_args2 = []
                        list_args2.append(fromPDBToCIF(
                            list_args[1], list_args[1].replace('.pdb', '.cif')))
                        args2 = list_args[0] + " " + list_args2[0]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args2,
                                                extraEnvDict=None, cwd=cwd)
                    except:
                        # upgrade cif
                        list_args3 = []
                        list_args0 = args2.split()
                        for i in range(0, 2):
                            list_args3[i].append(fromCIFTommCIF(
                                list_args0[i], list_args0[i]))
                        args3 = list_args3[0] + " " + list_args3[1]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE),
                                                args3, extraEnvDict=None, cwd=cwd)
            elif list_args[0].endswith(".pdb") and list_args[1].endswith(".cif"):
                try:
                    # pdbs: convert cif to pdb
                    list_args1 = []
                    list_args1.append(fromCIFToPDB(
                        list_args[1], list_args[1].replace('.cif', '.pdb')))
                    args1 = list_args[0] + " " + list_args1[0]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args1,
                                            extraEnvDict=None, cwd=cwd)
                except:
                    try:
                        # cifs: convert pdb to cif
                        list_args2 = []
                        list_args2.append(fromPDBToCIF(
                            list_args[0], list_args[0].replace('.pdb', '.cif')))
                        args2 = list_args2[0] + " " + list_args[1]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args2,
                                                extraEnvDict=None, cwd=cwd)
                    except:
                        # upgrade cifs
                        list_args3 = []
                        list_args0 = args2.split()
                        for i in range(0, 2):
                            list_args3.append(fromCIFTommCIF(
                                list_args0[i], list_args0[i]))
                        args3 = list_args3[0] + " " + list_args3[1]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE),
                                                args3, extraEnvDict=None, cwd=cwd)
        except:
            # biopython conversion
            aSH = AtomicStructHandler()
            try:
                for i in range(0, 2):
                    aSH.read(list_args[i])
                    aSH.write(list_args[i])
                    args = list_args[0] + " " + list_args[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE),
                                            args, extraEnvDict=None, cwd=cwd)
            except:
                print("CIF file standarization failed.")