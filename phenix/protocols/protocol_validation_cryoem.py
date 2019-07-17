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
from phenix.constants import MOLPROBITY, VALIDATION_CRYOEM, PHENIXVERSION
from phenix import Plugin
from pyworkflow.em.convert.atom_struct import toCIF
from protocol_refinement_base import PhenixProtRunRefinementBase

class PhenixProtRunValidationCryoEM(PhenixProtRunRefinementBase):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure inferred from an electron density map.
"""
    _label = 'validation_cryoem'
    _program = ""
    #_version = VERSION_1_2
    VALIDATIONCRYOEMFILE = 'validation_cryoem.mrc'
    VALIDATIONCRYOEMPKLFILE = 'validation_cryoem.pkl'

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtRunValidationCryoEM, self)._defineParams(form)
        param = form.getParam('inputVolume')
        param.help.set("\nSet the starting volume.\nPhenix will calculate "
                       "real-space correlation between map and atomic structure.\n"
                       "If the volume and atomic structure are not correctly "
                       "fitted, values of real-space correlation will indicate "
                       "not correlation at all.\n")

    # --------------------------- INSERT steps functions --------------------

    def _insertAllSteps(self):
        if self.inputVolume.get() is not None:
            self.vol = self.inputVolume.get()
        elif self.inputStructure.get().getVolume() is not None:
            self.vol = self.inputStructure.get().getVolume()
        if self.vol is not None:
            self._insertFunctionStep('convertInputStep', self.VALIDATIONCRYOEMFILE)
        if self.inputStructure.get().getFileName().endswith(".cif"):
            self._insertFunctionStep('sanitizeAtomStruct', self.inputStructure.get().getFileName())
        self._insertFunctionStep('runValidationCryoEMStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def sanitizeAtomStruct(self, fileName):

        # transform to CIF
        localCIFFileName = self._sanitizedStructureFileName(fileName)
        fileName = toCIF(fileName, localCIFFileName)

        # read file
        with open(fileName,'r') as f:
            content = f.read()

        import re
        f = open(localCIFFileName, "w")
        f.write(re.sub(' +', ' ', content)) # remove double spaces

        f.close()

    def runValidationCryoEMStep(self):
        version = Plugin.getPhenixVersion()
        if version == '1.13':
            print "PHENIX version: 1.13"
        else:
            print "PHENIX version: ", version

        # PDBx/mmCIF
        fileName = self.inputStructure.get().getFileName()
        if fileName.endswith(".cif"):
            pdb = self._sanitizedStructureFileName(fileName)
        else:
            pdb = os.path.abspath(fileName)

        # starting volume (.mrc)
        if self.vol is not None:
            tmpMapFile = self.VALIDATIONCRYOEMFILE
            volume = os.path.abspath(self._getExtraPath(tmpMapFile))
        else:
            volume = None

        # MolProbity is run to get the file molprobity.out
        # (necessary to get geometry outliers)
        numberOfThreads = self.numberOfThreads.get()

        args = " " + pdb
        if Plugin.getPhenixVersion() == PHENIXVERSION and volume is not None:
            args += (" map_file_name=%s" % volume) + \
                        (" d_min=%f" % self.resolution.get())
            args += " pickle=True"
        args += " pdb_interpretation.clash_guard.nonbonded_distance_threshold=None" \
                + (" %s " % self.extraParams.get())

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

        if Plugin.getPhenixVersion() != PHENIXVERSION and self.vol is not None:

            if self.vol.getHalfMaps():
                halves = []
                for halfMap in self.vol.getHalfMaps().split(','):
                    if not os.path.abspath(halfMap).endswith(".mrc"):
                        half = os.path.abspath(halfMap).split(".")[0] + ".mrc"
                    else:
                        half = os.path.abspath(halfMap)
                    halves.append(half)
                args = " " + pdb + (" " + volume) + (" " + halves[0]) + (" " + halves[1]) \
                       + " " + ("resolution=%f" % self.resolution.get()) + " pickle=True" + \
                       " slim=False" + \
                       " pdb_interpretation.clash_guard.nonbonded_distance_threshold=None" \
                       + (" %s " % self.extraParams.get())
            else:
                args = " " + pdb + (" " + volume) \
                   + " "+ ("resolution=%f" % self.resolution.get()) + " pickle=True" + \
                   " slim=False" + \
                   " pdb_interpretation.clash_guard.nonbonded_distance_threshold=None" \
                   + (" %s " % self.extraParams.get())
            Plugin.runPhenixProgram(Plugin.getProgram(VALIDATION_CRYOEM), args,
                                    cwd=self._getExtraPath())

    def createOutputStep(self):
        VALIDATIONCRYOEMPKLFILENAME = self._getExtraPath(
            self.VALIDATIONCRYOEMPKLFILE)
        self._readValidationPklFile(VALIDATIONCRYOEMPKLFILENAME)
        self._store()

    # --------------------------- INFO functions ---------------------------
    def _validate(self):
        errors = self.validateBase(VALIDATION_CRYOEM, 'VALIDATION_CRYOEM')

        # Check that the input volume exist
        if (self.inputVolume.get() or self.inputStructure.get().getVolume()) \
                is None:
            errors.append("Error: You should provide a volume.\n")

        return errors


    def _summary(self):
        summary = PhenixProtRunRefinementBase._summary(self)
        summary.append("https://www.phenix-online.org/documentation/"
                       "reference/validation_cryo_em.html")
        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")
        return methodsMsgs

    def _citations(self):
        return ['Afonine_2018']

    def _sanitizedStructureFileName(self, inFileName):
        if inFileName.endswith(".pdb") or inFileName.endswith(".ent"):
            inFileName = inFileName.replace(".pdb", ".cif").replace(".ent", ".cif")
        return os.path.abspath(self._getExtraPath(os.path.basename(inFileName)))




