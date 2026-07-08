# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
from pyworkflow.protocol.params import BooleanParam,  IntParam, EnumParam, StringParam, FloatParam, PointerParam
from phenix.constants import (AQUAREF, MOLPROBITY, VALIDATION_CRYOEM)
from . import PhenixProtRunRefinementBase

from pyworkflow.protocol.constants import LEVEL_ADVANCED

from pwem.convert.atom_struct import retry, fromCIFTommCIF
from phenix import Plugin
import re

PDB = 0
mmCIF = 1
OUTPUT_FORMAT = ['pdb', 'mmcif']


class PhenixProtAQuaRefRSR(PhenixProtRunRefinementBase):
    """
    Module that carries out refinement of bio-macromolecules utilizing chemical restraints QM calculations.
    """
    _label = 'AQuaRef real space refine'
    _program = ""
    ENGINES = ['aimnet2', 'mopac', 'torchani', 'xtb', 'orca']
    MODES = ['refine', 'opt', 'gtest']
    MINIMIZERS = ['lbfgsb', 'lbfgs']

    REALSPACEFILE = 'real_space.mrc'
    VALIDATIONCRYOEMPKLFILE = 'validation_cryoem.pkl'
    OUTPDB = "real_space_refined.pdb"

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtAQuaRefRSR, self)._defineParams(form)
        param = form.getParam('inputVolume')
        param.help.set("\nSet the starting volume.\nPhenix will refine the "
                       "atomic structure according to the volume density.\n"
                       "Volume and atomic structure have to be correctly fitted. "
                       "Otherwise, values of real-space correlation will indicate "
                       "not correlation at all.\n")
        form.addParam('engine',EnumParam,choices=self.ENGINES, default=0,
            label='QM engine: ',
            help='Quantum mechanical engine used during refinement.'
        )
        form.addParam('mode', EnumParam, choices=self.MODES, default=0,
            label='Mode'
        )
        form.addParam("macroCycles", IntParam, label="Macro cycles: ",
                      default=1,
                      help="Number of iterations of refinement.\nAlthough 5 "
                           "macro-cycles is usually sufficient, in cases in "
                           "which model geometry or/and model-to-map fit is "
                           "poor the use of more macro-cycles could be "
                           "helpful.\n")
        form.addParam('refineCycles',IntParam, default=5,
            label='Refinement cycles: '
        )
        form.addParam('refineSites', BooleanParam, default=True,
            label='Refine coordinates: '
        )
        form.addParam('refineADP', BooleanParam, default=False,
            label='Refine ADPs: ',
            help="Phenix default parameter.\nGenerally, refinement "
                "with all defaults is sufficient.\n\nADP ("
                "B-factors) refinement against the map is "
                "performed at the last macro-cycle only. "
        )
        form.addParam('exclude', StringParam, default='', label='Exclude selection: ',
            help='Example: element MG or name ATP. Separated by commas.'
        )

        group = form.addGroup('Optimization strategy options', expertLevel=LEVEL_ADVANCED)
        group.addParam('weightCycles', IntParam, default=10,
            label='Weight search cycles: '
        )
        group.addParam('skipWeightSearch', BooleanParam, default=False,
            label='Skip weight search: '
        )
        group.addParam('dataWeight', FloatParam, allowsNull=True,
            label='Data weight: '
        )
        group.addParam('minimizer', EnumParam, choices=self.MINIMIZERS, default=0,
            label='Minimizer'
        )

        form.addParallelSection(threads=1, mpi=0)


    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep', self.REALSPACEFILE)
        self._insertFunctionStep('runAQuaRefStep')
        self._insertFunctionStep('runMolprobityStep', self.REALSPACEFILE)
        self._insertFunctionStep('runValidationCryoEMStep', self.REALSPACEFILE)
        self._insertFunctionStep('createOutputStep')

    def runAQuaRefStep(self):
        Plugin.runQRefineProgram(
            "qr.aquaref",
            self.getExecArgs(),
            cwd=self._getExtraPath()
        )

    def runMolprobityStep(self, tmpMapFile):
        # PDBx/mmCIF
        outPdb = self._getExtraPath(self.OUTPDB)
        atomStruct = os.path.abspath(outPdb)
        # starting volume (.mrc)
        vol = os.path.abspath(self._getExtraPath(tmpMapFile))
        args = self._writeArgsMolProbity(atomStruct, vol)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY),
              args, cwd=cwd,
              listAtomStruct=[atomStruct], log=self._log,
              sdterrLog=self.getLogsLastLines)

    def runValidationCryoEMStep(self, tmpMapFile):
        # PDBx/mmCIF
        outPdb = self._getExtraPath(self.OUTPDB)
        atomStruct = os.path.abspath(outPdb)
        # starting volume (.mrc)
        volume = os.path.abspath(self._getExtraPath(tmpMapFile))
        if self.inputVolume.get() is not None:
            vol = self.inputVolume.get()
        else:
            vol = self.inputStructure.get().getVolume()

        args = self._writeArgsValCryoEM(atomStruct, volume, vol)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        retry(Plugin.runPhenixProgram, Plugin.getProgram(VALIDATION_CRYOEM),
              args, cwd=cwd,
              listAtomStruct=[atomStruct], log=self._log,
              sdterrLog=self.getLogsLastLines)

    def createOutputStep(self):
        pdb = AtomStruct()
        outPdb = self._getExtraPath(self.OUTPDB)
        pdb.setFileName(outPdb)

        if self.inputVolume.get() is not None:
            pdb.setVolume(self.inputVolume.get())
        else:
            pdb.setVolume(self.inputStructure.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructure.get(), pdb)
        if self.inputVolume.get() is not None:
            self._defineSourceRelation(self.inputVolume.get(), pdb)

        VALIDATIONCRYOEMPKLFILENAME = self._getExtraPath(
            self.VALIDATIONCRYOEMPKLFILE)
        self._readValidationPklFile(VALIDATIONCRYOEMPKLFILENAME)
        self._store()

    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = []
        return errors

    def _citations(self):
        return []

    def _summary(self):
        summary = PhenixProtRunRefinementBase._summary(self)
        summary.append(
            "https://www.phenix-online.org/documentation/reference/"
            "real_space_refine.html")
        return summary

    # --------------------------- UTILS functions --------------------------
    def getExecArgs(self):
        atomStruct = os.path.abspath(self.inputStructure.get().getFileName())
        volume = os.path.abspath(self.inputVolume.get().getFileName())

        args = f"{atomStruct} {volume}"
        args += f" quantum.engine_name={self.ENGINES[self.engine.get()]}"
        args += f" refine.mode={self.MODES[self.mode.get()]}"
        args += f" refine.number_of_macro_cycles={self.macroCycles.get()}"
        args += f" refine.number_of_refine_cycles={self.refineCycles.get()}"
        if not self.refineSites.get():
            args += " refine.refine_sites=False"
        if self.refineADP.get():
            args += " refine.refine_adp=True"

        args += (f" refine.number_of_weight_search_cycles="f"{self.weightCycles.get()}")
        if self.skipWeightSearch.get():
            args += " refine.skip_weight_search=True"
        if self.dataWeight.hasValue():
            args += f" refine.data_weight={self.dataWeight.get()}"
        args += f" refine.minimizer={self.MINIMIZERS[self.minimizer.get()]}"

        if self.exclude.get().strip():
            args += f' refine.exclude="{self.exclude.get()}"'
        args += f" parallel.nproc={self.numberOfThreads.get()}"

        return args
