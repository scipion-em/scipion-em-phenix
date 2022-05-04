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
from phenix.constants import DOCKANDREBUILD, PHENIX_HOME
from pwem.convert.atom_struct import fromCIFToPDB, fromPDBToCIF, \
    fromCIFTommCIF, AtomicStructHandler, retry

try:
    from pwem.objects import AtomStruct, Sequence
except:
    from pwem.objects import PdbFile as AtomStruct
from phenix import Plugin


class PhenixProtDockAndRebuildAlphaFold2Model(EMProtocol):
    """Rebuild Docked Predicted Model.
     Rebuild predicted model morphs and rebuilds a model produced by AlphaFold,
     RoseTTAFold and other prediction software into a cryo EM map, using a set
     of docked domains from the predicted model as a template.
    """
    _label = 'dock and rebuild predicted model'
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
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input map', important=True,
                      help="Set the starting density map.")
        # form.addParam('asymmetricMap', BooleanParam, default=True,
        #               label='Asymmetric map:',
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
        self._insertFunctionStep('runDockAndBuildModel')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def runDockAndBuildModel(self):
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        mapFile = self.DOCKINMAPFILE
        localVolName = os.path.abspath(self._getExtraPath(mapFile))
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        Ccp4Header.fixFile(inVolName, localVolName, origin, sampling, Ccp4Header.START)  # ORIGIN
        predictedAtomStruct = os.path.abspath(
            self.inputPredictedModel.get().getFileName())

        predictedAtomStruct_localPath = os.path.abspath(
            self._getExtraPath(predictedAtomStruct.split('/')[-1]))
        if str(predictedAtomStruct) != str(predictedAtomStruct_localPath):
            pwutils.path.createLink(predictedAtomStruct, predictedAtomStruct_localPath)
            predictedAtomStruct = predictedAtomStruct_localPath

        self.prefix = os.path.abspath(self._getExtraPath(predictedAtomStruct))
        args = self._writeArgsDockAlphaFold(
            predictedAtomStruct, localVolName, self.prefix)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(DOCKANDREBUILD),
              args, cwd=cwd,
              listAtomStruct=[predictedAtomStruct],
              log=self._log)

    def createOutputStep(self):
        pdb = AtomStruct()
        nameProcessed = ''
        nameProcessed = self.prefix.split("/")[-1]
        for fileName in os.listdir(self._getExtraPath()):
            if (fileName.endswith("pdb") and
                    len(fileName.split(".")) > len(nameProcessed.split("."))):
                print("fileName: ", fileName)
                pdb.setFileName(self._getExtraPath(fileName))
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputPredictedModel.get(), pdb)

        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram(DOCKANDREBUILD)
        print("program: ", program)
        if not os.path.exists(program):
            errors.append("Cannot find " + program)

            # If there is any error at this point it is related to config variables
            errors.append("Check configuration file: " +
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set PHENIX_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % Plugin.getVar(PHENIX_HOME))
                errors.append("PROCESS = %s" % DOCKANDREBUILD)

        return errors

    def _summary(self):
        summary = []
        # try:
        #     summary.append("protocol finished with results")
        # except:
        #     summary.append("processed predicted model not yet docked")
        summary.append(
            "https://phenix-online.org/version_docs/dev-4380/reference/dock_and_rebuild.html")

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
            self, predictedAtomStruct, vol, prefix):
        args = " "
        args += "model=%s " % predictedAtomStruct
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
