# **************************************************************************
# *
# * Authors:     Roberto Marabini
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
import math
import os
import sqlite3

import matplotlib.pyplot as plt
from tkinter import messagebox
from phenix.protocols.protocol_search_fit import PhenixProtSearchFit
from pyworkflow.protocol.params import LabelParam, EnumParam, IntParam, FloatParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pwem.viewers import TableView, Chimera
from pwem import Domain
from phenix.protocols.protocol_search_fit import (DATAFILE,
                                                  TABLE)
from phenix import Plugin
from phenix.constants import PHENIXVERSION19

def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        messagebox.showerror("Error",  # bar title
                              msg,  # message
                              parent=tkParent)
    except:
        print(("Error:", msg))


class PhenixProtRuSearchFitViewer(ProtocolViewer):
    """ Viewer for Phenix program DEARCHFIT
    """
    _label = 'Search Fit Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [PhenixProtSearchFit]
    SEARCHFITSUBPLOTSFILENAME = 'model_to_map_fit.py'

    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)

    def _defineParams(self, form):
        form.addSection(label="Volume and models")
        form.addParam('showMapModel', LabelParam,
                      label="Show volume and atomic structures in ChimeraX",
                      help="Display input volume and atom struct"
                           " plus 5 better atom struct candidates.")
        form.addParam("numAtomStruct", IntParam, label="Max. Number Atom Structs.",
                      default=1000,
                      help="Number of atom structs to show ordered by model_to_map_fit\n")
        form.addParam("zone", FloatParam, label="Show Area around input atomic struct (A)",
                      default=3,
                      help="Limit the display to a zone around the input atomic structure.\n"
                           "Units = A.")
        form.addParam('showPlot', LabelParam,
                      label="Summary Plot",
                      help="Plot showing 'model_to_map_fit' values. The X axis indicates \n"
                           "the XXXX number of each refined atomic structure that appears \n"
                           "in files named \n"
                           "coot_000000_Imol_0000_version_XXXX_real_space_refined_000.cif")

    def _getVisualizeDict(self):
        return {
            'showMapModel': self._showMapModel,
            'showPlot': self._showPlot
        }

    def _showMapModel(self, e=None):
        bildFileName = self.protocol._getExtraPath("axis_output.bild")

        _inputVol = self.protocol.inputVolume.get()
        dim = _inputVol.getDim()[0]
        sampling = _inputVol.getSamplingRate()
        counter = 1

        Chimera.createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)
        fnCmd = os.path.abspath(self.protocol._getExtraPath("chimera_output.cxc"))
        f = open(fnCmd, 'w')
        # change to workingDir
        # If we do not use cd and the project name has an space
        # the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())
        # reference axis model = 0
        f.write("open %s\n" % os.path.abspath(bildFileName)) # 1
        f.write("cofr 0,0,0\n")  # set center of coordinates

        # input 3D map
        counter += 1  # 2
        vol = self._getInputVolume()
        fnVol = os.path.abspath(vol.getFileName())
        f.write("open %s\n" % fnVol)
        x, y, z = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        f.write("volume #%d style surface voxelSize %f\nvolume #%d origin "
                "%0.2f,%0.2f,%0.2f\n"
                % (counter, sampling, counter, x, y, z))

        # input PDB (usually from coot)
        counter += 1  # 3
        pdbFileName = os.path.abspath(self.protocol.inputStructure.get().getFileName())
        f.write("open %s\n" % pdbFileName)
        # zone command
        f.write("volume zone #%d near #%d range %f newMap false\n" %
                (counter -1, counter, self.zone.get()))
        f.write("cofr #%d\n" % counter)
        # open database and retrieve all files
        conn = sqlite3.connect(os.path.abspath(self.protocol._getExtraPath(DATAFILE)))
        c = conn.cursor()
        sqlCommand = """SELECT filename, phenix_id
                        FROM   %s
                        WHERE model_to_map_fit != -1
                        ORDER BY model_to_map_fit DESC
                        LIMIT %d""" % (TABLE, self.numAtomStruct)
        c.execute(sqlCommand)
        rows = c.fetchall()

        for row in rows:
            if Plugin.getPhenixVersion() >= PHENIXVERSION19:
                atomStructFn = row[0][:-4] + "_real_space_refined_000.cif"
            else:
                atomStructFn = row[0][:-4] + "_real_space_refined.cif"
            f.write("open %s\n" % atomStructFn)
        c.close()
        conn.close()
        f.close()
        # run in the background
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&",
                                        cwd=os.getcwd())
        return []

    def _getInputVolume(self):
        if self.protocol.inputVolume.get() is None:
            fnVol = self.protocol.inputStructure.get().getVolume()
        else:
            fnVol = self.protocol.inputVolume.get()
        return fnVol


    def _showPlot(self, e=None):

        xList = []
        yList = []
        conn = sqlite3.connect(os.path.abspath(self.protocol._getExtraPath(DATAFILE)))
        c = conn.cursor()
        sqlCommand = """
                        SELECT id, model_to_map_fit
                            FROM (
                                  SELECT id, model_to_map_fit
                                  FROM   %s
                                  WHERE model_to_map_fit != -1
                                  ORDER BY model_to_map_fit desc
                                  LIMIT %d) AS a
                           ORDER BY id
                          """ % (TABLE, self.numAtomStruct)

        c.execute(sqlCommand)
        rows = c.fetchall()

        for row in rows:
            xList.append(float(row[0]) - 1.0)
            yList.append(float(row[1]))
        # compute avg
        sqlCommand = """SELECT AVG(model_to_map_fit) 
                        FROM (SELECT model_to_map_fit
                             FROM %s
                             WHERE model_to_map_fit != -1
                             ORDER BY model_to_map_fit desc
                             )""" % (TABLE)
        c.execute(sqlCommand)
        rows = c.fetchone(); avg = rows[0]
        print("avg", avg)
        fromRelation = """(SELECT model_to_map_fit
                             FROM %s
                             WHERE model_to_map_fit != -1
                             ORDER BY model_to_map_fit desc
                             ) AS mainTable""" % (TABLE)

        sqlCommand = """SELECT AVG((mainTable.model_to_map_fit - sub.a) * (mainTable.model_to_map_fit - sub.a)) as var
                        FROM %s,
                        (SELECT AVG(model_to_map_fit) AS a
                         FROM %s
                         WHERE model_to_map_fit != -1
                        ) AS sub
                        WHERE model_to_map_fit != -1
                        """ % (fromRelation, fromRelation)
        c.execute(sqlCommand)
        rows = c.fetchone(); std_2 = rows[0]
        std = math.sqrt(std_2)
        c.close()
        conn.close()

        if not (xList or yList):
            errorWindow(self.getTkRoot(), "No data available")
            return

        title = 'avg (all data) = %f, std (all data) = %f'%(avg, std)
        plt.plot(xList, yList, 'x')
        plt.axis([-1.0, max(xList) + 1.0, 0.0, max(yList)+0.1])
        plt.title(title)
        plt.xlabel('#Atom Structs')
        plt.ylabel('Map Model Fit Score')
        plt.show()
        print(xList, yList)

        # SELECT FLOOR(model_to_map_fit/5.00)*5 As Grade,
        #        COUNT(*) AS [Grade Count]
        # FROM TableName
        # GROUP BY FLOOR(model_to_map_fit/5.00)*5
        # ORDER BY 1
