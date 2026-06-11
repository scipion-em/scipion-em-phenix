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
from pwem.convert import Ccp4Header, headers
from pwem.protocols import EMProtocol
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, BooleanParam, EnumParam,
                                        StringParam, FloatParam, IntParam,
                                        MultiPointerParam)
from phenix.constants import PREDICTANDBUILD, PHENIX_HOME, VERSION
from pwem.convert.atom_struct import fromCIFToPDB, fromPDBToCIF, \
    fromCIFTommCIF, AtomicStructHandler, retry
from pwem.objects import Volume, Sequence, SetOfSequences

try:
    from pwem.objects import AtomStruct
except:
    from pwem.objects import PdbFile as AtomStruct
from phenix import Plugin


_version = VERSION
OUTPUT_BEST_MAP = 'predict_and_build_best_map.mrc'
OUTPUT_BEST_PDB = 'predict_and_build_best_pdb.pdb'
OUTPUT_BEST_SUPERPOSED_PREDICTED_MODELS = 'predict_and_build_best_superposed_predicted_models.pdb'
OUTPUT_SUPERPOSED_UNTRIMMED_CYCLE_1 = 'predict_and_build_superposed_predicted_untrimmed_cycle_1.pdb'
OUTPUT_SUPERPOSED_UNTRIMMED_CYCLE_2 = 'predict_and_build_superposed_predicted_untrimmed_cycle_2.pdb'
OUTPUT_SEQUENCE_FILE = 'seq_file.fasta'

class PhenixPredictAndBuildCryoEM(EMProtocol):
    """Predict, dock and rebuild AlphaFold models into a cryo-EM map. 
    This method works with a sequence file and 2
    half maps or one full map.
    """
    _label = 'predict and build cyoEM'
    _program = ""
    _possibleOutputs = {
        OUTPUT_BEST_MAP: Volume,
        OUTPUT_BEST_PDB: AtomStruct,
        OUTPUT_BEST_SUPERPOSED_PREDICTED_MODELS: AtomStruct,
        OUTPUT_SUPERPOSED_UNTRIMMED_CYCLE_1: AtomStruct,
        OUTPUT_SUPERPOSED_UNTRIMMED_CYCLE_2:AtomStruct,
        OUTPUT_SEQUENCE_FILE:SetOfSequences
    }


    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('useHalfMapsInsteadVol', BooleanParam, default=False,
                      label="Would you like to use half maps?",
                      help='Phenix Predict and Build cryoEM uses either half maps or the full map.\n'
                            'Please, select the type of input map(s) you will provide.\n'
                             'If you select half maps, map optimization by density modification will be carried out.\n')

        form.addParam('halfMapsAttached', BooleanParam, default=True,
                      condition='useHalfMapsInsteadVol',
                      label="Are the half maps included in the volume?",
                      help='When you import a map, you can associate half maps to it. Select *yes* if the half maps are associated'
                           'to the input volume. If half maps are not associated, select *No* and'
                           'you will be able to provide then as regular maps')
        form.addParam('inputHalf1', PointerParam, pointerClass='Volume',
                      label="Volume Half 1", important=True,
                      condition='useHalfMapsInsteadVol and not halfMapsAttached',
                      help='Select half map 1')

        form.addParam('inputHalf2', PointerParam, pointerClass='Volume',
                      label="Volume Half 2", important=True,
                      condition='useHalfMapsInsteadVol and not halfMapsAttached',
                      help='Select half map 2')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input map', important=True,
                      condition='not useHalfMapsInsteadVol or halfMapsAttached',
                      help="Set the starting full density map.")
        # form.addParam('asymmetricMap', BooleanParam, default=True,
        #               label='Asymmetric map:',
        #               help="If your map has symmetry be sure to set this param No."
        #                    " Otherwise symmetry will be automatically determined.")
        form.addParam('inputSequenceS', MultiPointerParam,
                      pointerClass="Sequence", important=True,
                      label='Protein sequences',
                      help="Include here one or more sequences to predict the AlphaFold structure\n"
                            "Sequences of all the unique chains in your molecule and how many of these are present.\n"
                            "These can be supplied as a .fasta file.")
        form.addParam('resolution', FloatParam, default=3.0,
                      label='High-resolution limit (A):',
                      help="Map resolution (Angstroms).")
        form.addParam('inputPredictedModel', PointerParam, pointerClass='AtomStruct', 
                      expertLevel=LEVEL_ADVANCED, allowsNull=True,
                      label='Predicted model. (optional)',
                      help="Set the atomic structure obtained in any way by yourself.\n"
                           "Supported formats are PDB or mmCIF; this last one"
                           " is especially useful for very large structures.")
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
        self._insertFunctionStep('runPredictAndBuildCryoEM')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def convertInputStep(self):
        """ Read the input volume."""
        self.input_half1_fn = None
        self.input_half2_fn = None
        self.input_vol_fn = None

        if self.useHalfMapsInsteadVol.get():
            if self.halfMapsAttached.get():
                vol = self.inputVolume.get()
                self.input_half1_fn, self.input_half2_fn = vol.getHalfMaps().split(',')
                self.input_half1_fn = os.path.abspath(self.input_half1_fn)
                self.input_half2_fn = os.path.abspath(self.input_half2_fn)
            else:
                half1Vol = self.inputHalf1.get()
                half2Vol = self.inputHalf2.get()
                self.input_half1_fn = os.path.abspath(self._inputVol2Mrc(half1Vol))
                self.input_half2_fn = os.path.abspath(self._inputVol2Mrc(half2Vol))

        else:
            vol = self.inputVolume.get()
            self.input_vol_fn = os.path.abspath(self._inputVol2Mrc(vol))

    def _writeArgsPredictAndBuild(self, fastafilename):
        prefix = self._getExtraPath().split("/")[-2]
        self.prefix = f'PredictAndBuid_{prefix.split("_")[0]}_rebuilt'
        print("prefix3", prefix)
        args = " "
        if self.input_half1_fn is not None:
            args += "half_map=%s " % self.input_half1_fn
            args += "half_map=%s " % self.input_half2_fn
        else:
            args += "full_map=%s " % self.input_vol_fn
        args += "seq_file=%s " % os.path.abspath(fastafilename)
        args += "crystal_info.resolution=%f " % self.resolution.get()
        # add templates
        args += "output_model_prefix=%s " % self.prefix
        if self.numberOfThreads > 1:
            print("self.numberOfThreads: ", self.numberOfThreads)
            args += "nproc=%d " % self.numberOfThreads
            if len(str(self.extraParams)) > 0:
                args += " %s " % self.extraParams.get()
        return args

    def runPredictAndBuildCryoEM(self):
        fastaFileName = self.createInputFastaFile()
        args = self._writeArgsPredictAndBuild(fastafilename=fastaFileName)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(PREDICTANDBUILD),
        args, cwd=cwd,
        listAtomStruct=[],
        log=self._log, sdterrLog = self.getLogsLastLines)

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
        program = Plugin.getProgram(PREDICTANDBUILD)
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
                    errors.append("PROCESS = %s" % PREDICTANDBUILD)

        return errors

    def _summary(self):
        summary = []
            # try:
                # summary.append("protocol finished with results")
            # except:
                # summary.append("processed predicted model not yet docked")
        summary.append(
            "https://phenix-online.org/version_docs/dev-4380/reference/dock_and_rebuild.html")

        return summary

    def _citations(self):
        return ['Terwilliger_2022']

# # --------------------------- UTILS functions --------------------------

    def _inputVol2Mrc(self, vol):
        """ convert 3D maps to MRC '.mrc' format
        """
        volFileName = vol.getFileName()
        if headers.getFileFormat(volFileName) == headers.MRC:
            mrcFileName = volFileName
        else:
            mrcFileName = os.path.basename(volFileName)
            mrcFileName = pwutils.replaceBaseExt(mrcFileName, "mrc")
            mrcFileName = self._getExtraPath(mrcFileName)
            origin = vol.getOrigin(force=True).getShifts()
            sampling = vol.getSamplingRate()
            Ccp4Header.fixFile(volFileName, mrcFileName, origin, sampling, Ccp4Header.START) # ORIGIN
        return mrcFileName

    def createInputFastaFile(self):
        """ Get sequence as string and create the corresponding fasta file. """
        inputSeqs = self.inputSequenceS
        print("inputSeqs", inputSeqs)
        fastaFileName = self._getExtraPath(OUTPUT_SEQUENCE_FILE)

        with open(fastaFileName, "w") as f:
            for seq in inputSeqs:
                s = seq.get()
                f.write(f"> {s.getId()}\n")
                f.write(f"{s.getSequence()}\n")

        return fastaFileName

# def _registerAtomStruct(self, name, path):
# if not os.path.exists(path):
# raise FileNotFoundError("Output %s not found." % path)

# output = AtomStruct(filename=path)
# self._defineOutputs(**{name: output})
# self._defineSourceRelation(self.inputVolume, output)

# seqs = self.inputSequenceS
# if seqs:
# for seq in seqs:
# self._defineSourceRelation(seq, output)
