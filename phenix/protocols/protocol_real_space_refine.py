# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

from pwem.objects import AtomStruct
from pyworkflow.protocol.params import BooleanParam,  IntParam
from phenix.constants import (REALSPACEREFINE,
                              MOLPROBITY2,
                              VALIDATION_CRYOEM,
                              PHENIXVERSION,
                              PHENIXVERSION19,
                              PHENIXVERSION20)

from pyworkflow.protocol.constants import LEVEL_ADVANCED

from pwem.convert.atom_struct import retry, fromCIFTommCIF
from .protocol_refinement_base import PhenixProtRunRefinementBase
from phenix import Plugin
import re

PDB = 0
mmCIF = 1
OUTPUT_FORMAT = ['pdb', 'mmcif']


class PhenixProtRunRSRefine(PhenixProtRunRefinementBase):
    """Tool for extensive real-space refinement of an atomic structure
    against the map provided. The map can be derived from X-ray or neutron
    crystallography, or cryoEM. The program obtains a model that fits the map
    as well as possible having appropriate geometry. The model should not show
    validation outliers, such as Ramachandran plot or rotamer outliers.
    """
    _label = 'real space refine'
    _program = ""
    # _version = VERSION_1_2
    REALSPACEFILE = 'real_space.mrc'
    if Plugin.getPhenixVersion() != PHENIXVERSION:
        VALIDATIONCRYOEMPKLFILE = 'validation_cryoem.pkl'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtRunRSRefine, self)._defineParams(form)
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

        # form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', self.REALSPACEFILE)
        self._insertFunctionStep('runRSrefineStep', self.REALSPACEFILE)
        self._insertFunctionStep('runMolprobityStep', self.REALSPACEFILE)
        if Plugin.getPhenixVersion() != PHENIXVERSION:
            self._insertFunctionStep('runValidationCryoEMStep', self.REALSPACEFILE)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------
    def runRSrefineStep(self, tmpMapFile):
        atomStruct = os.path.abspath(self.inputStructure.get().getFileName())
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        args = self._writeArgsRSR(atomStruct, vol)
        cwd = os.getcwd() + "/" + self._getExtraPath()

        retry(Plugin.runPhenixProgram,
              Plugin.getProgram(REALSPACEREFINE), args,
              # cwd=os.path.abspath(self._getExtraPath()),
              cwd=cwd,
              listAtomStruct=[atomStruct], log=self._log)
        self.refinedFile = False
        for item in os.listdir(self._getExtraPath()):
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
                  cwd=cwd,
                  listAtomStruct=[atomStruct], log=self._log)

    def runMolprobityStep(self, tmpMapFile):
        # PDBx/mmCIF
        self._getRSRefineOutput()
        atomStruct = os.path.abspath(self.outAtomStructName)
        # starting volume (.mrc)
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        args = self._writeArgsMolProbity(atomStruct, vol)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY2),
              # args, cwd=os.path.abspath(self._getExtraPath()),
              args, cwd=cwd,
              listAtomStruct=[atomStruct], log=self._log)

    def runValidationCryoEMStep(self, tmpMapFile):
        # PDBx/mmCIF
        atomStruct = os.path.abspath(self.outAtomStructName)
        # starting volume (.mrc)
        volume = os.path.abspath(self._getExtraPath(tmpMapFile))
        if self.inputVolume.get() is not None:
            vol = self.inputVolume.get()
        else:
            vol = self.inputStructure.get().getVolume()

        args = self._writeArgsValCryoEM(atomStruct, volume, vol)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(VALIDATION_CRYOEM),
              # args, cwd=os.path.abspath(self._getExtraPath()),
              args, cwd=cwd,
              listAtomStruct=[atomStruct], log=self._log)

    def createOutputStep(self):
        # self._getRSRefineOutput()
        pdb = AtomStruct()
        pdb.setFileName(self.outAtomStructName)

        if self.inputVolume.get() is not None:
            pdb.setVolume(self.inputVolume.get())
        else:
            pdb.setVolume(self.inputStructure.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), pdb)
        if self.inputVolume.get() is not None:
            self._defineSourceRelation(self.inputVolume.get(), pdb)

        if Plugin.getPhenixVersion() == PHENIXVERSION:
            MOLPROBITYOUTFILENAME = self._getExtraPath(
                self.MOLPROBITYOUTFILENAME)
            self._parseFile(MOLPROBITYOUTFILENAME)
        else:
            VALIDATIONCRYOEMPKLFILENAME = self._getExtraPath(
                self.VALIDATIONCRYOEMPKLFILE)
            self._readValidationPklFile(VALIDATIONCRYOEMPKLFILENAME)
        self._store()
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

    def _getRSRefineOutput(self):
        inPdbName = os.path.basename(self.inputStructure.get().getFileName())
        list_cif = []
        for item in os.listdir(self._getExtraPath()):
            p = re.compile('\d+')
            if p.search(item) is not None and item.endswith(".cif"):
                list_cif.append(item)
        name = sorted(list_cif)[-1]
        outAtomStructName = self._getExtraPath(name)
        # convert cif to mmcif by using maxit program
        # to get the right number and name of chains
        log = self._log
        self.outAtomStructName = outAtomStructName
        fromCIFTommCIF(outAtomStructName, self.outAtomStructName, log)

    def _writeArgsRSR(self, atomStruct, vol):
        if Plugin.getPhenixVersion() == PHENIXVERSION19 or PHENIXVERSION20:
            args = " "
        else:
            args = " model_file="
        args += "%s " % atomStruct
        if Plugin.getPhenixVersion() == PHENIXVERSION19 or PHENIXVERSION20:
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
