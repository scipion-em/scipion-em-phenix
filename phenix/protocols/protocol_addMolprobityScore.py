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
from phenix.constants import MOLPROBITY, PHENIXVERSION
from phenix import Plugin
from .protocol_refinement_base import PhenixProtRunRefinementBase
from pwem.convert.atom_struct import retry

from pwem.objects import AtomStruct, SetOfAtomStructs, Volume
from pwem.convert.headers import Ccp4Header
from pyworkflow.protocol.params import (PointerParam, FloatParam, StringParam, LEVEL_ADVANCED, STEPS_PARALLEL)

class PhenixProtAddMolprobity(PhenixProtRunRefinementBase):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure inferred from an electron density map. This protocol adds its score to
the attribute of the AtomStructures in the set
"""
    _label = 'Add molprobity'
    _program = ""
    #_version = VERSION_1_2
    MOLPROBITYFILE = 'molprobity.mrc'
    TMPCIFFILENAME="inMolprobity.cif"
    TMPPDBFILENAME="inMolprobity.pdb"

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

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------

    def _insertAllSteps(self):
        cIds = []
        if self._getInputVolume() is not None:
            cIds.append(self._insertFunctionStep('convertInputStep', self.MOLPROBITYFILE))

        molProbIds = []
        for inStructFn in self.getAtomStructFileNames():
            molProbIds.append(self._insertFunctionStep('runMolprobityStep', inStructFn, prerequisites=cIds))

        self._insertFunctionStep('createOutputStep', prerequisites=molProbIds)

    # --------------------------- STEPS functions --------------------------
    def convertInputStep(self, tmpMapFileName):
        """ convert 3D maps to MRC '.mrc' format
        """
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        newFn = self._getExtraPath(tmpMapFileName)
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        Ccp4Header.fixFile(inVolName, newFn, origin, sampling, Ccp4Header.START)  # ORIGIN


    def runMolprobityStep(self, inStructFn):
        version = Plugin.getPhenixVersion()
        if version == '1.13':
            print("PHENIX version: 1.13")
        else:
            print(("PHENIX version: ", version))
        # PDBx/mmCIF

        outDir = self._getExtraPath(self.getAtomStructBasename(inStructFn))
        os.mkdir(outDir)
        if self._getInputVolume() is not None:
            self.volFn = os.path.abspath(self._getExtraPath(self.MOLPROBITYFILE))
            args = self._writeArgsMolProbityExpand(inStructFn, self.volFn)
        else:
            args = self._writeArgsMolProbityExpand(inStructFn, vol=None)
        # script with auxiliary files
        retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY),
              # args, cwd=os.path.abspath(self._getExtraPath()),
              args, cwd=outDir,
              listAtomStruct=[inStructFn], log=self._log)

    def createOutputStep(self):
        outStructs = SetOfAtomStructs.create(self._getPath())
        outVol = Volume(self._getExtraPath(self.MOLPROBITYFILE))
        outVol.copyInfo(self._getInputVolume())
        for inStructFn in self.getAtomStructFileNames():
            outDir = self._getExtraPath(self.getAtomStructBasename(inStructFn))
            MOLPROBITYOUTFILENAME = os.path.join(outDir, self.MOLPROBITYOUTFILENAME)
            try:
                self._parseFile(MOLPROBITYOUTFILENAME)
            except:
                if self.MOLPROBITYFILE is not None:
                    # self.vol = os.path.abspath(self._getExtraPath(self.MOLPROBITYFILE))
                    self.volFn = self._getExtraPath(self.MOLPROBITYFILE)
                    args = self._writeArgsMolProbityExpand(inStructFn, self.volFn)
                else:
                    args = self._writeArgsMolProbityExpand(inStructFn, vol=None)
                args += " allow_polymer_cross_special_position=True "
                retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY),
                      # args, cwd=os.path.abspath(self._getExtraPath()),
                      args, cwd=outDir,
                      listAtomStruct=[inStructFn], log=self._log)
                self._parseFile(MOLPROBITYOUTFILENAME)

            newAS = AtomStruct(filename=inStructFn)
            newAS.setVolume(outVol)
            newAS = self.addMolProbScores(newAS)
            outStructs.append(newAS)

        self._defineOutputs(outputAtomStructs=outStructs)
        self._defineSourceRelation(self.inputStructureSet, outStructs)



    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = self.validateBase(MOLPROBITY, 'MOLPROBITY')
        return errors

    def _summary(self):
        summary = PhenixProtRunRefinementBase._summary(self)
        summary.append("MolProbity: http://molprobity.biochem.duke.edu/")
        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")
        return methodsMsgs

    def _citations(self):
        return ['Chen_2010']

    # --------------------------- UTILS functions --------------------------

    def _writeArgsMolProbityExpand(self, atomStruct, vol=None):
        args = self._writeArgsMolProbity(atomStruct, vol)
        if Plugin.getPhenixVersion() != PHENIXVERSION:
            args += " pickle=True"
        args += " pdb_interpretation.clash_guard.nonbonded_distance_threshold=None"
        args += " %s " % self.extraParams.get()
        # args += " wxplots=True" # TODO: Avoid the direct opening of plots
        return args

    def getAtomStructFileNames(self):
        fns = []
        for inStruct in self.inputStructureSet.get():
            fns.append(os.path.abspath(inStruct.getFileName()))
        return fns

    def getAtomStructBasename(self, atomStructFn):
        return os.path.splitext(os.path.basename(atomStructFn))[0]

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
