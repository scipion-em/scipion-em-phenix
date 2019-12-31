# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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

class PhenixProtRunMolprobity(PhenixProtRunRefinementBase):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure inferred from an electron density map.
"""
    _label = 'molprobity'
    _program = ""
    #_version = VERSION_1_2
    MOLPROBITYFILE = 'molprobity.mrc'
    TMPCIFFILENAME="inMolprobity.cif"
    TMPPDBFILENAME="inMolprobity.pdb"

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtRunMolprobity, self)._defineParams(form)
        param = form.getParam('inputVolume')
        param.help.set("\nSet the starting volume.\nOnly with version 1.13, "
                       "Phenix will calculate real-space correlation.\n"
                       "If the volume and atomic structure are not correctly "
                       "fitted, values of real-space correlation will indicate "
                       "not correlation at all.\n")

    # --------------------------- INSERT steps functions --------------------

    def _insertAllSteps(self):
        if (self.inputVolume.get() or self.inputStructure.get().getVolume()) \
                is not None:
            self._insertFunctionStep('convertInputStep', self.MOLPROBITYFILE)
        self._insertFunctionStep('runMolprobityStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def runMolprobityStep(self):
        version = Plugin.getPhenixVersion()
        if version == '1.13':
            print("PHENIX version: 1.13")
        else:
            print(("PHENIX version: ", version))
        # PDBx/mmCIF
        fileName = self.inputStructure.get().getFileName()
        # self.atomStruct = os.path.abspath(fileName)
        self.atomStruct = os.getcwd() + "/" + fileName
        # starting volume (.mrc)
        if (self.inputVolume.get() or self.inputStructure.get().getVolume()) \
                is not None:
            tmpMapFile = self.MOLPROBITYFILE
            # self.vol = os.path.abspath(self._getExtraPath(tmpMapFile))
            self.vol = os.getcwd() + "/" + self._getExtraPath(tmpMapFile)
            args = self._writeArgsMolProbityExpand(self.atomStruct, self.vol)
        else:
            args = self._writeArgsMolProbityExpand(self.atomStruct, vol=None)
        # script with auxiliary files
        retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY),
              # args, cwd=os.path.abspath(self._getExtraPath()),
              args, cwd=self._getExtraPath(),
              listAtomStruct=[self.atomStruct], log=self._log)

    def createOutputStep(self):
        MOLPROBITYOUTFILENAME = self._getExtraPath(
            self.MOLPROBITYOUTFILENAME)
        try:
            self._parseFile(MOLPROBITYOUTFILENAME)
        except:
            if self.MOLPROBITYFILE is not None:
                # self.vol = os.path.abspath(self._getExtraPath(self.MOLPROBITYFILE))
                self.vol = self._getExtraPath(self.MOLPROBITYFILE)
                args = self._writeArgsMolProbityExpand(self.atomStruct, self.vol)
            else:
                args = self._writeArgsMolProbityExpand(self.atomStruct, vol=None)
            args += " allow_polymer_cross_special_position=True "
            retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY),
                  # args, cwd=os.path.abspath(self._getExtraPath()),
                  args, cwd=self._getExtraPath(),
                  listAtomStruct=[self.atomStruct], log=self._log)
            self._parseFile(MOLPROBITYOUTFILENAME)
        self._store()

    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = self.validateBase(MOLPROBITY,'MOLPROBITY')
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
