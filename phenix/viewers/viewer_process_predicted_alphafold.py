# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
# *           Marta Martinez (mmmtnez@cnb.csic.es)
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

import os

from phenix.protocols import PhenixProtProcessPredictedAlphaFold2Model
from pyworkflow.viewer import Viewer
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pwem.viewers import Chimera
from pwem import Domain
from phenix import Plugin
from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.gui.text import openTextFileEditor

class PhenixProtRunProcessPredictedAlphaFoldViewer(ProtocolViewer):
    """ Visualize the output of protocol process predicted Alphafold2 model """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'Process predicted alphafold2 model viewer'
    _targets = [PhenixProtProcessPredictedAlphaFold2Model]
    REMAINDER = False

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        for fileName in os.listdir(self.protocol._getExtraPath()):
            if fileName.endswith(self.protocol.SEQREMAINDER):
                self.fileName = fileName
                self.REMAINDER = True

    def _defineParams(self, form):
        form.addSection(label='Visualize results')
        form.addParam('visualizeStructures', LabelParam,
                      label="Structures in ChimeraX",
                      help="Clik the eye to visualize the input predicted"
                           " model by Alphafold2, as well as the processed one.\n")
        if self.REMAINDER:
            form.addParam('RemainSequences', LabelParam,
                          label="Remaining Sequences",
                          help="Open a text file with the sequence fragments "
                               "removed from the processed one.\n")

    def _getVisualizeDict(self):
        return {
            'visualizeStructures': self._displayStructures,
            'RemainSequences': self._visualizeRemainSequences
        }

    def _displayStructures(self, e=None):
        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")
        # To show pdbs only
        dim = 150.
        sampling = 1.

        bildFileName = self.protocol._getExtraPath("axis_output.bild")
        Chimera.createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)

        with open(fnCmd, 'w') as f:
            # change to workingDir
            # If we do not use cd and the project name has an space
            # the protocol fails even if we pass absolute paths
            f.write('cd %s\n' % os.getcwd())
            f.write("open %s\n" % bildFileName) #1
            f.write("cofr 0,0,0\n")  # set center of coordinates

            pdbList = self._getPdbs()
            for filename in pdbList:
                f.write("open %s\n" % filename)
            f.write("color bfactor  #%d  palette alphafold" % 2)

        # run in the background
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&",
                                        cwd=os.getcwd())
        return []

    def _visualizeRemainSequences(self, e=None):
        # Visualization of extra file
        openTextFileEditor(os.path.abspath(
            self.protocol._getExtraPath(self.fileName)))

    def _getPdbs(self):
        pdbList = []
        pdbInput = self.protocol.inputPredictedModel.get().getFileName()
        pdbList.append(pdbInput)
        pdbProcessed = self.protocol.outputPdb.getFileName()
        pdbList.append(pdbProcessed)
        return pdbList





