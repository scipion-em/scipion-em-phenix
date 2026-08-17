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
from phenix.constants import MOLPROBITY
from phenix import Plugin
from .protocol_refinement_base import PhenixProtRunRefinementBase
from pwem.convert.atom_struct import retry

class PhenixProtRunMolprobity(PhenixProtRunRefinementBase):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure inferred from an electron density map.
"""
    """
        MolProbity (PhenixProtRunMolprobity) — User Manual

            Overview

            The MolProbity protocol validates the stereochemical quality and
            structural consistency of an atomic model obtained from cryo-EM data.
            Using the Phenix implementation of MolProbity, the protocol evaluates
            whether the fitted structure is geometrically reliable and biologically
            plausible.

            Inputs and Workflow

            The protocol requires an atomic structure in PDBx/mmCIF format and can
            optionally use a cryo-EM density map. If a map is provided, the
            protocol also evaluates the agreement between the model and the
            experimental density through real-space correlation analysis.

            During execution, the input map is converted into MRC format when
            needed, and the Phenix MolProbity validation program is launched
            automatically. The workflow includes stereochemical validation,
            geometry analysis, and optional density-based validation depending on
            the installed Phenix version.

            Validation and Biological Interpretation

            MolProbity helps detect structural problems such as steric clashes,
            unrealistic bond geometries, backbone outliers, and incorrect side
            chain conformations. These metrics are essential for assessing model
            quality before publication, deposition, or downstream biological
            analysis.

            When density information is available, the protocol also evaluates how
            well the model fits the experimental map. Poor correlation values may
            indicate domain misplacement, incorrect residue assignment, or local
            overfitting.

            Outputs

            The protocol produces validation statistics and structural quality
            metrics that help determine the reliability of the atomic model.
            Results can be used to identify problematic regions requiring further
            refinement or rebuilding.

            Final Perspective

            MolProbity is an essential validation step in cryo-EM workflows because
            accurate biological interpretation depends not only on map fitting, but
            also on maintaining chemically and stereochemically correct atomic
            models.
        """
    _label = 'molprobity'
    _program = ""
    MOLPROBITYFILE = 'molprobity.mrc'
    TMPCIFFILENAME = "inMolprobity.cif"
    TMPPDBFILENAME = "inMolprobity.pdb"

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

        fileName = self.inputStructure.get().getFileName()
        self.atomStruct = os.path.join(os.getcwd(), fileName)
        # self.atomStruct = os.path.abspath(self.atomStruct)
        # self.atomStruct = os.getcwd() + "/" + fileName
        # starting volume (.mrc)
        if (self.inputVolume.get() or self.inputStructure.get().getVolume()) \
                is not None:
            tmpMapFile = self.MOLPROBITYFILE
            # self.vol = os.path.abspath(self._getExtraPath(tmpMapFile))
            self.vol = os.path.join(os.getcwd(), self._getExtraPath(tmpMapFile))
            # self.vol = os.getcwd() + "/" + self._getExtraPath(tmpMapFile)
            args = self._writeArgsMolProbityExpand(self.atomStruct, self.vol)
        else:
            args = self._writeArgsMolProbityExpand(self.atomStruct, vol=None)
        # script with auxiliary files
        retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY),
              # args, cwd=os.path.abspath(self._getExtraPath()),
              args, cwd=self._getExtraPath(),
              listAtomStruct=[self.atomStruct], log=self._log,
              sdterrLog = self.getLogsLastLines)

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
                  listAtomStruct=[self.atomStruct], log=self._log,
                  sdterrLog = self.getLogsLastLines)

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
        args += " pickle=True"
        args += " pdb_interpretation.clash_guard.nonbonded_distance_threshold=None"
        args += " %s " % self.extraParams.get()
        # args += " wxplots=True" # TODO: Avoid the direct opening of plots
        return args
