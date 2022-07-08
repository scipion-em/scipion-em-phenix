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

from phenix.protocols.protocol_dock_in_map import PhenixProtRunDockInMap
from phenix.protocols.protocol_dock_predicted_alphafold2_model import \
    PhenixProtDockPredictedAlphaFold2Model
from phenix.protocols.protocol_rebuild_docked_predicted_alphafold2_model import \
    PhenixProtRebuildDockPredictedAlphaFold2Model
from phenix.protocols.protocol_dock_and_rebuild_alphafold_model import \
    PhenixProtDockAndRebuildAlphaFold2Model
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from pwem.viewers import Chimera
from pwem import Domain
from phenix import Plugin

class PhenixProtRunDockInMapViewer(Viewer):
    """ Visualize the output of protocol dock in map """
    _environments = [DESKTOP_TKINTER]
    _label = 'Dock in map viewer'
    _targets = [PhenixProtRunDockInMap, PhenixProtDockPredictedAlphaFold2Model,
                PhenixProtRebuildDockPredictedAlphaFold2Model,
                PhenixProtDockAndRebuildAlphaFold2Model]

    def _visualize(self, obj, **args):
        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")

        self._getVols()
        self._getPdbs()
        dim = float()
        sampling = float()
        if self.vols[0] is not None:
            dim, sampling = self._getDimSamplingFromVol(self.vols[0])
        else:
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
            if len(self.vols) > 0:
                for vol in self.vols:
                    sampling, volFileName, x, y, z = self._getXYZFromVol(vol)
                    f.write("open %s\n" % volFileName)
                    f.write("volume #2 style surface voxelSize %f\n"
                            "volume #2 origin %0.2f,%0.2f,%0.2f\n"
                            % (sampling, x, y, z))
            for filename in self.pdbList:
                f.write("open %s\n" % filename)
            if self.protocol.hasAttribute("inputVolume"):
                f.write("color bfactor  #%d  palette alphafold" % 3)
            elif not self.protocol.hasAttribute("inputVolume") and \
                self.protocol.hasAttribute("inputPredictedModel"):
                f.write("color bfactor  #%d  palette alphafold" % 2)
                # f.write("color bfactor palette alphafold\n")

        # run in the background
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&",
                                        cwd=os.getcwd())
        return []

    def _getVols(self):
        self.vols = []
        if self.protocol.hasAttribute("inputVolume1"):
            vol1 = self.protocol.inputVolume1.get()
            if vol1 is not None:
                self.vols.append(vol1)
        if self.protocol.hasAttribute("inputVolume"):
            vol2 = self.protocol.inputVolume.get()
            if vol2 is not None:
                self.vols.append(vol2)

    def _getPdbs(self):
        self.pdbList = []
        if self.protocol.hasAttribute("inputVolume1"):
            inputPdb1 = self.protocol.inputStructure.get().getFileName()
            self.pdbList.append(inputPdb1)
        if self.protocol.hasAttribute("inputVolume"):
            if self.protocol.hasAttribute("inputPredictedModel"):
                inputPdb2 = self.protocol.inputPredictedModel.get().getFileName()
                self.pdbList.append(inputPdb2)
            if self.protocol.hasAttribute("inputProcessedPredictedModel"):
                inputPdb3 = self.protocol.inputProcessedPredictedModel.get().getFileName()
                self.pdbList.append(inputPdb3)
            if self.protocol.hasAttribute("inputDockedPredictedModel"):
                inputPdb4 = self.protocol.inputDockedPredictedModel.get().getFileName()
                self.pdbList.append(inputPdb4)
        outputPdb = self.protocol.outputPdb.getFileName()
        self.pdbList.append(outputPdb)

    def _getDimSamplingFromVol(self, vol):
        dim = vol.getDim()[0]
        sampling = vol.getSamplingRate()

        return dim, sampling

    def _getXYZFromVol(self, vol):
        sampling = vol.getSamplingRate()
        volFileName = vol.getFileName()
        if vol.hasOrigin():
            x, y, z = vol.getOrigin().getShifts()
        else:
            x, y, z = vol.getOrigin(force=True).getShifts()
        return sampling, volFileName, x, y, z



