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

from pwem.convert import Ccp4Header
from pwem.convert.atom_struct import retry, fromPDBToCIF, fromCIFTommCIF

from pyworkflow import Config
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, FloatParam, IntParam, LEVEL_ADVANCED
from phenix.constants import  PHENIX_HOME, DOCKINMAP, DISPLAY
from os.path import relpath

try:
    from pwem.objects import AtomStruct
except:
    from pwem.objects import PdbFile as AtomStruct
from phenix import Plugin


class PhenixProtRunDockInMap(EMProtocol):
    """Docking of a PDB (one or several copies) into a map """
    _label = 'dock in map'
    _program = ""
    # _version = VERSION_1_2
    DOCKINMAPFILE = 'dock_in_map.mrc'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume1', PointerParam, pointerClass="Volume",
                      label='Input map', important=True,
                      help="Set the starting density map.")
        form.addParam('resolution', FloatParam, default=3.0,
                      label='Resolution (A):',
                      help="Map resolution (Angstroms)."
                           "Use at least the double of the sampling rate ("
                           "Angstroms/pixel)")
        form.addParam('inputStructure', PointerParam,
                      pointerClass="AtomStruct", important=True,
                      label='Input atom structure',
                      help="PDBx/mmCIF to be fitted against the volume. ")
        form.addParam('modelCopies', IntParam,
                      default=1,
                      label='Atom structure number of copies',
                      help="Write here the number of copies of your atom structure. ")
        form.addParam('numberOfThreads', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1,
                      label='Number of threads',
                      help="Write here the number of threads to run the protocol. ")

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runDockInMapStep')
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

    def runDockInMapStep(self):
        # starting structure
        atomStructFileName = self.inputStructure.get().getFileName()
        # atomStruct = os.getcwd() + "/" + atomStructFileName
        atomStruct = atomStructFileName
        # starting map (.mrc)
        mapFile = self.DOCKINMAPFILE
        vol = os.getcwd() + "/" + self._getExtraPath(mapFile)
        args = self._writeArgsDocKInMap(vol, atomStruct)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(DOCKINMAP),
            # args, cwd=os.path.abspath(self._getExtraPath()),
            args, cwd=cwd,
            listAtomStruct=[atomStruct], log=self._log)

    def createOutputStep(self):
        self._getDockInMapOutput()
        pdb = AtomStruct()
        pdb.setFileName(relpath(self.outAtomStructName))

        if self.inputVolume1.get() is not None:
            pdb.setVolume(self.inputVolume1.get())
        else:
            pdb.setVolume(self.inputStructure.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), pdb)
        if self.inputVolume1.get() is not None:
            self._defineSourceRelation(self.inputVolume1.get(), pdb)

        # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram(DOCKINMAP)
        if program is None:
            errors.append("Missing variables DOCKINMAP and/or PHENIX_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)
        elif not self.is_tool(DISPLAY):
            errors.append("display program missing.\n "
                          "Please install imagemagick package")

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: " +
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set DOCKINMAP and PHENIX_HOME variables "
                          "properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % Plugin.getVar(PHENIX_HOME))
                errors.append("DOCKINMAP = %s" % DOCKINMAP)

        # Check that the input volume exist
        if self._getInputVolume() is None:
            errors.append("Error: You should provide a map.\n")
        if self.inputStructure is None:
            errors.append("Error: You should provide an atomic structure to fit.\n")

        return errors

    def _summary(self):
        #  Think on how to update this summary with created PDB
        summary = []
        # try:
        #     dataDict = json.loads(str(self.stringDataDict))
        #     summary.append("Optimal Threshold: %0.2f   Rotamer-Ratio: %0.2f"
        #                    % (dataDict['Optimal Threshold'],
        #                       dataDict['Rotamer-Ratio']))
        #     summary.append("Max Zscore:         %0.2f  Model Length: %d"
        #                    % (dataDict['Max Zscore'],
        #                       dataDict['Model Length']))
        #     summary.append("EMRinger Score:  %0.2f"
        #                    % dataDict['EMRinger Score'])
        # except:
        #     summary = ["EMRinger Score not yet computed"]

        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")

        return methodsMsgs

    def _citations(self):
        return ['Barad_2015']

        # --------------------------- UTILS functions --------------------------

    def _getInputVolume(self):
        if self.inputVolume1.get() is None:
            fnVol = self.inputStructure.get().getVolume()
        else:
            fnVol = self.inputVolume1.get()
        return fnVol

    def is_tool(self, name):
        """Check whether `name` is on PATH."""
        from shutil import which
        return which(name) is not None

    def _getDockInMapOutput(self):
        outAtomStructName = os.getcwd() + "/" +\
                            self._getExtraPath("placed_model.pdb")
        # convert cif to mmcif by using maxit program
        # to get the right number and name of chains
        log = self._log
        self.outAtomStructName = outAtomStructName.replace("pdb", "cif")
        fromPDBToCIF(outAtomStructName, self.outAtomStructName, log)
        fromCIFTommCIF(self.outAtomStructName, self.outAtomStructName, log)

    def _writeArgsDocKInMap(self, vol, atomStruct):
        args = ""
        args += " map_file=%s" % vol
        args += " resolution=%f" % self.resolution
        args += " search_model=%s" % atomStruct
        if self.modelCopies > 1:
            args += " search_model_copies=%d" % self.modelCopies
            args += " use_symmetry=False"
        if self.numberOfThreads > 1:
            args += " nproc=%d" % self.numberOfThreads
        return args

