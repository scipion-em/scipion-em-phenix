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
from pyworkflow.object import Float, Integer
from pyworkflow.utils import importFromPlugin
from phenix.constants import MOLPROBITY
from phenix import Plugin


from protocol_refinement_base import PhenixProtRunRefinementBase

class PhenixProtRunMolprobity(PhenixProtRunRefinementBase):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure derived from a cryo-EM density map.
"""
    _label = 'molprobity'
    _program = ""
    #_version = VERSION_1_2
    MOLPROBITYFILE = 'molprobity.mrc'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtRunMolprobity, self)._defineParams(form)

    # --------------------------- INSERT steps functions --------------------

    def _insertAllSteps(self):
        if (self.inputVolume.get() or self.inputStructure.get().getVolume()) \
                is not None:
            self._insertFunctionStep('convertInputStep', self.MOLPROBITYFILE)
        if self.inputStructure.get().getFileName().endswith(".cif"):
            self._insertFunctionStep('sanitizeAtomStruct', self.inputStructure.get().getFileName())
        self._insertFunctionStep('runMolprobityStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def sanitizeAtomStruct(self, fileName):
        """molprobity does not like CIF files produced by biopython"""
        import re
        content = open(fileName,"r").read()
        f = open(self._sanitizedStructureFileName(fileName), "w")
        f.write(re.sub(' +', ' ', content))
        f.close()

    def runMolprobityStep(self):
        # PDBx/mmCIF
        fileName = self.inputStructure.get().getFileName()
        if fileName.endswith(".cif"):
            pdb = self._sanitizedStructureFileName(fileName)
        else:
            pdb = os.path.abspath(fileName)
        args = ""
        args += pdb
        # starting volume (.mrc)
        if (self.inputVolume.get() or self.inputStructure.get().getVolume()) \
                is not None:
            tmpMapFile = self.MOLPROBITYFILE
            volume = os.path.abspath(self._getExtraPath(tmpMapFile))
            args += " "
            args += "map_file_name=%s" % volume
            args += " "
            args += "d_min=%f" % self.resolution.get()
        args += " "
        args += "pickle=True"
        args += " pdb_interpretation.clash_guard.nonbonded_distance_threshold=None"
        args += " %s " % self.extraParams.get()

        numberOfThreads = self.numberOfThreads.get()
        if numberOfThreads > 1:
            args += " nproc=%d" % numberOfThreads
        # args += " wxplots=True" # TODO: Avoid the direct opening of plots
        # script with auxiliary files
        try:
            Plugin.runPhenixProgram(Plugin.getProgram(MOLPROBITY), args,
                         cwd=self._getExtraPath())
        except:
            print "WARNING!!!\nPHENIX error:\n pdb_interpretation.clash_guard" \
                  " failure: High number of nonbonded interaction distances " \
                  "< 0.5. This error has been disable by running the same " \
                  "command with the same following additional " \
                  "argument:\npdb_interpretation.clash_guard." \
                  "nonbonded_distance_threshold=None "
            args += " "
            args += "pdb_interpretation.clash_guard." \
                    "nonbonded_distance_threshold=None"
            Plugin.runPhenixProgram(Plugin.getProgram(MOLPROBITY), args,
                             cwd=self._getExtraPath())


    def createOutputStep(self):
        MOLPROBITYOUTFILENAME = self._getExtraPath(
            self.MOLPROBITYOUTFILENAME)
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

    def _sanitizedStructureFileName(self, fileName):
        return os.path.abspath(self._getExtraPath(os.path.basename(fileName)))




