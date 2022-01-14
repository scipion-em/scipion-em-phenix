# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
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

from pwem.objects import AtomStruct, SetOfAtomStructs
from phenix.constants import (REALSPACEREFINE,
                              MOLPROBITY2,
                              VALIDATION_CRYOEM,
                              PHENIXVERSION,
                              PHENIXVERSION19)

from pyworkflow.protocol.params import PointerParam, FloatParam, StringParam, LEVEL_ADVANCED,\
  BooleanParam,  IntParam, STEPS_PARALLEL

from pwem.convert.atom_struct import retry, fromCIFTommCIF
from .protocol_refinement_base import PhenixProtRunRefinementBase
from phenix import Plugin
import re

PDB = 0
mmCIF = 1
OUTPUT_FORMAT = ['pdb', 'mmcif']


class PhenixProtRunRSRefineSet(PhenixProtRunRefinementBase):
    """Tool for extensive real-space refinement of an atomic structure
    against the map provided. The map can be derived from X-ray or neutron
    crystallography, or cryoEM. The program obtains a model that fits the map
    as well as possible having appropriate geometry. The model should not show
    validation outliers, such as Ramachandran plot or rotamer outliers.
    """
    _label = 'Set real space refine'
    _program = ""
    # _version = VERSION_1_2
    REALSPACEFILE = 'real_space.mrc'
    if Plugin.getPhenixVersion() != PHENIXVERSION:
        VALIDATIONCRYOEMPKLFILE = 'validation_cryoem.pkl'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="")
        form.addParam('resolution', FloatParam,
                      label='Resolution (A):',
                      default=3.0,
                      help='Set the resolution of the input volume.')
        form.addParam('inputStructureSet', PointerParam,
                      pointerClass="SetOfAtomStructs", allowsNull=False,
                      label='Input atomic structures.',
                      help="Set the atomic structure to be processed.\n"
                           "Supported formats are PDB or mmCIF; this last one"
                           " is especially useful for very large structures.")
        form.addParam('extraParams', StringParam,
                      label="Extra Params ",
                      default="",
                      expertLevel=LEVEL_ADVANCED,
                      help="This string will be added to the phenix command.\n"
                           "Syntax: paramName1=value1 paramName2=value2 ")

        param = form.getParam('inputVolume')
        param.help.set("\nSet the starting volume.\nPhenix will refine the "
                       "atomic structure according to the volume density.\n"
                       "Volume and atomic structure have to be correctly fitted. "
                       "Otherwise, values of real-space correlation will indicate "
                       "not correlation at all.\n")
        form.addParam("doSecondary", BooleanParam, label="Secondary structure",
                      default=True, expertLevel=LEVEL_ADVANCED,
                      help="Set to TRUE to use secondary structure "
                           "restraints.\nOnly for PHENIX versions higher than 1.13.")
        form.addParam("macroCycles", IntParam, label="Macro cycles",
                      default=5, expertLevel=LEVEL_ADVANCED,
                      help="Number of iterations of refinement.\nAlthough 5 "
                           "macro-cycles is usually sufficient, in cases in "
                           "which model geometry or/and model-to-map fit is "
                           "poor the use of more macro-cycles could be "
                           "helpful.\n")
        group = form.addGroup('Optimization strategy options')
        group.addParam('minimizationGlobal', BooleanParam,
                       label="Global minimization: ", default=True,
                       expertLevel=LEVEL_ADVANCED,
                       help="Phenix default parameter to look for the global "
                            "minimum of the model.\nGenerally, refinement "
                            "with all defaults is "
                            "sufficient.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('rigidBody', BooleanParam,
                       label="Rigid body: ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Refinement strategy that considers groups of "
                            "atoms that move (rotate and translate) as a "
                            "single body.\n")
        group.addParam('localGridSearch', BooleanParam,
                       label="Local grid search: ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Refinement strategy that considers "
                            "local rotamer fitting.\n\n Generally, refinement "
                            "with all defaults is sufficient.\n Including "
                            "local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('morphing', BooleanParam,
                       label="Morphing ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Morphing procedure distorts a model to match an "
                            "electron density map.\n\nGenerally, refinement "
                            "with all defaults is "
                            "sufficient.\n Including local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('simulatedAnnealing', BooleanParam,
                       label="Simulated annealing ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Optimization technique known as molecular "
                            "dynamics refinement; it minimizes the energy of "
                            "the model.\n"
                            "Generally, refinement with all defaults is "
                            "sufficient.\n Including local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('adp', BooleanParam,
                       label="Atomic Displacement Parameters (ADPs) ",
                       default=True,
                       expertLevel=LEVEL_ADVANCED,
                       help="Phenix default parameter.\nGenerally, refinement "
                            "with all defaults is sufficient.\n\nADP ("
                            "B-factors) refinement against the map is "
                            "performed at the last macro-cycle only. ")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        convId = self._insertFunctionStep('convertInputStep', self.REALSPACEFILE, prerequisites=[])
        refineIds = []
        for inStructFn in self.getAtomStructFileNames():
            refineIds.append(self._insertFunctionStep('runRSrefineStep', self.REALSPACEFILE,
                                                      inStructFn, prerequisites=[convId]))
          
        molProbIds = []
        for outDir in self.getOutputDirs():
            molProbIds.append(self._insertFunctionStep('runMolprobityStep', self.REALSPACEFILE,
                                                       outDir, prerequisites=refineIds))
          
        if Plugin.getPhenixVersion() != PHENIXVERSION:
            valIds = []
            for outDir in self.getOutputDirs():
                valIds.append(self._insertFunctionStep('runValidationCryoEMStep', self.REALSPACEFILE,
                                                       outDir, prerequisites=molProbIds))
            molProbIds = valIds
        self._insertFunctionStep('createOutputStep', prerequisites=molProbIds)

    # --------------------------- STEPS functions --------------------------
    def runRSrefineStep(self, tmpMapFile, inStructFn):
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        
        outDir = self._getExtraPath(self.getAtomStructBasename(inStructFn))
        os.mkdir(outDir)

        args = self._writeArgsRSR(inStructFn, vol)

        retry(Plugin.runPhenixProgram,
              Plugin.getProgram(REALSPACEREFINE), args,
              # cwd=os.path.abspath(self._getExtraPath()),
              cwd=outDir,
              listAtomStruct=[inStructFn], log=self._log)
        self.refinedFile = False
        for item in os.listdir(outDir):
            p = re.compile('\d+')
            if p.search(item) is not None and item.endswith(".cif"):
                self.refinedFile = True
                break

        if self.refinedFile == False:
            print("WARNING!!!\nPHENIX error:\n pdb_interpretation.clash_guard" \
                  " failure: High number of nonbonded interaction distances " \
                  "< 0.5. This error has been disable by running the same " \
                  "command with the same following additional " \
                  "argument:\npdb_interpretation.clash_guard." \
                  "nonbonded_distance_threshold=None ")
            args += " pdb_interpretation.clash_guard." \
                    "nonbonded_distance_threshold=None"
            retry(Plugin.runPhenixProgram,
                  Plugin.getProgram(REALSPACEREFINE), args,
                  # cwd=os.path.abspath(self._getExtraPath()),
                  cwd=outDir,
                  listAtomStruct=[inStructFn], log=self._log)

    def runMolprobityStep(self, tmpMapFile, outDir):
        # PDBx/mmCIF
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        outRefined = self._getRSRefineOutput(outDir)
        outRefined = os.path.abspath(outRefined)

        args = self._writeArgsMolProbity(outRefined, vol)
        retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY2),
              # args, cwd=os.path.abspath(self._getExtraPath()),
              args, cwd=outDir,
              listAtomStruct=[outRefined], log=self._log)

    def runValidationCryoEMStep(self, tmpMapFile, outDir):
        # PDBx/mmCIF
        volume = os.path.abspath(self._getExtraPath(tmpMapFile))
        iniVol = self._getInputVolume()
        outRefined = self._getRSRefineOutput(outDir, mmCIF=False)
        outRefined = os.path.abspath(outRefined)

        args = self._writeArgsValCryoEM(outRefined, volume, iniVol)
        retry(Plugin.runPhenixProgram, Plugin.getProgram(VALIDATION_CRYOEM),
              # args, cwd=os.path.abspath(self._getExtraPath()),
              args, cwd=outDir,
              listAtomStruct=[outRefined], log=self._log)

    def createOutputStep(self):
        outStructs = SetOfAtomStructs.create(self._getPath())
        outVol = self._getInputVolume()

        for inStructFn in self.getRefinedAtomStructFileNames():
            outDir = os.path.dirname(inStructFn)
            newAS = AtomStruct(filename=inStructFn)
            newAS.setVolume(outVol)

            if Plugin.getPhenixVersion() == PHENIXVERSION:
                MOLPROBITYOUTFILENAME = os.path.join(outDir, self.MOLPROBITYOUTFILENAME)
                self._parseFile(MOLPROBITYOUTFILENAME)
            else:
                VALIDATIONCRYOEMPKLFILENAME = os.path.join(outDir, self.VALIDATIONCRYOEMPKLFILE)
                self._readValidationPklFile(VALIDATIONCRYOEMPKLFILENAME)
            newAS = self.addMolProbScores(newAS)
            outStructs.append(newAS.clone())

        self._defineOutputs(outputAtomStructs=outStructs)
        self._defineSourceRelation(self.inputStructureSet, outStructs)

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = self.validateBase(REALSPACEREFINE, 'REALSPACEREFINE')

        # Check that the input volume exist
        if self._getInputVolume() is None:
            errors.append("Error: You should provide a volume.\n")
        return errors

    def _citations(self):
        return ['Barad_2015']

    def _summary(self):
        summary = PhenixProtRunRefinementBase._summary(self)
        summary.append(
            "https://www.phenix-online.org/documentation/reference/"
            "real_space_refine.html")
        return summary

    # --------------------------- UTILS functions --------------------------

    def _getRSRefineOutput(self, outDir, mmCIF=True):
        list_cif = []
        for item in os.listdir(outDir):
            p = re.compile('\d+')
            if p.search(item) is not None and item.endswith(".cif"):
                list_cif.append(item)
        name = sorted(list_cif)[-1]
        outAtomStructName = os.path.join(outDir, name)
        # convert cif to mmcif by using maxit program
        # to get the right number and name of chains
        log = self._log
        outAtomStructNamemmCIF = outAtomStructName
        if mmCIF:
          fromCIFTommCIF(outAtomStructName, outAtomStructNamemmCIF, log)
          return outAtomStructNamemmCIF
        else:
          return outAtomStructName

    def _writeArgsRSR(self, atomStruct, vol):
        if Plugin.getPhenixVersion() == PHENIXVERSION19:
            args = " "
        else:
            args = " model_file="
        args += "%s " % atomStruct
        if Plugin.getPhenixVersion() == PHENIXVERSION19:
            args += " "
        else:
            args += " map_file="
        args += "%s " % vol
        args += " resolution=%f" % self.resolution
        args += " secondary_structure.enabled=%s" % self.doSecondary
        args += " run="
        if self.minimizationGlobal == True:
            args += "minimization_global+"
        if self.rigidBody == True:
            args += "rigid_body+"
        if self.localGridSearch == True:
            args += "local_grid_search+"
        if self.morphing == True:
            args += "morphing+"
        if self.simulatedAnnealing == True:
            args += "simulated_annealing+"
        if self.adp == True:
            args += "adp+"
        args = args[:-1]
        # args += " run=minimization_global+local_grid_search+morphing+simulated_annealing"
        args += " macro_cycles=%d" % self.macroCycles
        args += " model_format=pdb+mmcif"
        # args += " write_pkl_stats=True"
        args += " %s " % self.extraParams.get()
        numberOfThreads = self.numberOfThreads.get()
        if numberOfThreads > 1:
            args += " nproc=%d" % numberOfThreads
        return args

    def getAtomStructFileNames(self):
        fns = []
        for inStruct in self.inputStructureSet.get():
            fns.append(os.path.abspath(inStruct.getFileName()))
        return fns

    def getOutputDirs(self):
        outDirs = []
        for inStruct in self.getAtomStructFileNames():
            outDirs.append(self._getExtraPath(self.getAtomStructBasename(inStruct)))
        return outDirs

    def getRefinedAtomStructFileNames(self):
        fns = []
        for inStruct in self.getAtomStructFileNames():
            outDir = self._getExtraPath(self.getAtomStructBasename(inStruct))
            fns.append(self._getRSRefineOutput(outDir, mmCIF=False))
        return fns


    def getAtomStructBasename(self, atomStructFn):
        return os.path.splitext(os.path.basename(atomStructFn))[0].split('_real_space_refined')[0]

    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputStructureSet.get().getFirstItem().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol

    def addMolProbScores(self, newAS):
        '''Add the MolProbity scores to the AtomStruct'''
        newAS.molProbScore = self.overallScore
        newAS.clashScore = self.clashscore
        return newAS
