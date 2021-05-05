# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import pyworkflow.wizard as pwizard
from pwem.wizards import GetStructureChainsWizard, pwobj, emconv
from .protocols import PhenixProtSearchFit
from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog


class GetChainResiduesWizardPhenix(pwizard.Wizard):
    _targets = [(PhenixProtSearchFit, ['firstaa'])]

    def getSequenceStep(self, form):
        residueList = []
        if form.protocol.inputSequence.get() is not None:
            sequence = list(form.protocol.inputSequence.get().getSequence())
            for seq_index, seq_residue in enumerate(sequence):
                residueList.append('{"residue": %d, "%s"}' %
                                   (seq_index + 1, seq_residue))
        finalResiduesList = []
        for i in residueList:
            finalResiduesList.append(pwobj.String(i))
        return finalResiduesList

    def show(self, form, *params):
        finalResiduesList = self.getSequenceStep(form)
        provider = ListTreeProviderString(finalResiduesList)
        dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                                "Select one residue (residue number, "
                                "residue name)")
        form.setVar('firstaa', dlg.values[0].get())

class GetChainResidues2WizardPhenix(GetChainResiduesWizardPhenix):
    _targets = [(PhenixProtSearchFit, ['lastaa'])]

    def show(self, form, *params):
        finalResiduesList = self.getSequenceStep(form)
        provider = ListTreeProviderString(finalResiduesList)
        dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                                "Select one residue (residue number, "
                                "residue name)")
        form.setVar('lastaa', dlg.values[0].get())

