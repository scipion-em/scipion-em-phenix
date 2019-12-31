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
from phenix.constants import MOLPROBITY, VALIDATION_CRYOEM, PHENIXVERSION
from phenix import Plugin
from pwem.convert.atom_struct import retry
from .protocol_refinement_base import PhenixProtRunRefinementBase

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
        self._insertFunctionStep('runValidationCryoEMStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def runValidationCryoEMStep(self):
        version = Plugin.getPhenixVersion()
        if version == '1.13':
            print("PHENIX version: 1.13")
        else:
            print("PHENIX version: ", version)

        fileName = self.inputStructure.get().getFileName()
        atomStruct = os.path.abspath(fileName)

        # starting volume (.mrc)
        if self.vol is not None:
            tmpMapFile = self.VALIDATIONCRYOEMFILE
            volume = os.path.abspath(self._getExtraPath(tmpMapFile))
        else:
            volume = None

        # MolProbity is run to get the file molprobity.out
        # (necessary to get geometry outliers)
        numberOfThreads = self.numberOfThreads.get()

        args = self._writeArgsMolProbity(atomStruct, vol=volume)
        cwd = os.getcwd() + "/" + self._getExtraPath()
        # script with auxiliary files
        retry(Plugin.runPhenixProgram, Plugin.getProgram(MOLPROBITY),
              args, cwd=cwd, listAtomStruct=[atomStruct], log=self._log)

        args = self._writeArgsValCryoEM(atomStruct, volume, self.vol)

        if Plugin.getPhenixVersion() != PHENIXVERSION and self.vol is not None:
            retry(Plugin.runPhenixProgram, Plugin.getProgram(VALIDATION_CRYOEM),
                  args, cwd=cwd, listAtomStruct=[atomStruct], log=self._log)

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
