# ***************************************************************************
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
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
# ***************************************************************************/

# protocol to test the phenix protocol dock in map
import os
from phenix.protocols import PhenixProtRunDockInMap
from pwem.protocols.protocol_import import (ProtImportPdb,
                                            ProtImportVolumes)
from pyworkflow.tests import *
from chimera.protocols import ChimeraProtOperate
from xmipp3.protocols import XmippProtExtractUnit


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import atomic structures(PDBx/mmCIF files)
    """
    pdbID = '5ni1'

    def _importVolume(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/emd_3488.map'),
                'samplingRate': 1.05,
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_3488.map\n')
        self.launchProtocol(protImportVol)
        volume1 = protImportVol.outputVolume
        return volume1

    def _importStructure(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': self.pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1')
        self.launchProtocol(protImportPDB)
        return protImportPDB


class TestProtDockInMap(TestImportData):
    """ Test the protocol Dock In Map
    """

    def testDockInMap1(self):
        """ This test checks that phenix dock in map protocol runs with
        an atomic structure and a map"""
        print("Run phenix dock_in_map protocol from an imported atomic "
              "structure file and an imported map")

        # import PDB
        protImportPDB = self._importStructure()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        # import map
        map = self._importVolume()
        self.assertTrue(map.getFileName())

        args = {
                'inputVolume1': map,
                'resolution': 3.2,
                'inputStructure': structure1
               }

        protDockInMap = self.newProtocol(PhenixProtRunDockInMap, **args)
        protDockInMap.setObjLabel('DockInMap\n'
                                  'haemoglobin\n')
        self.launchProtocol(protDockInMap)
        self.assertTrue(os.path.exists(
            protDockInMap.outputPdb.getFileName()))

    def testDockInMap2(self):
        """ This test checks that phenix dock in map protocol runs with
        the unit cell of a map and a chain of the atomic structure"""
        print("Run phenix dock_in_map protocol from the unit cell "
              "of a map and a chain of the atomic structure")
        # import PDB from Database
        protImportPDB = self._importStructure()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        # import map
        map = self._importVolume()
        self.assertTrue(map.getFileName())

        # extract chain using chimera
        # create auxiliary CMD file for chimera operate
        extraCommands = ""
        extraCommands += "sel #2/A\n"
        extraCommands += "save /tmp/chainA.cif format mmcif models #2 relModel #1 selectedOnly true\n"
        extraCommands += "open /tmp/chainA.cif\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_A_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure1
                }
        protChimera = self.newProtocol(ChimeraProtOperate,
                                       **args)
        protChimera.setObjLabel('chimera operate\n extract chain A')
        self.launchProtocol(protChimera)
        chainA = eval("protChimera.DONOTSAVESESSION_A_Atom_struct__3_%06d" % \
                      protChimera.getObjId())
        self.assertIsNotNone(chainA.getFileName(),
                             "There was a problem with the alignment")

        # extract unit cell from map
        args = {'inputVolumes': map,
                'symmetryGroup': 0,
                'symmetryOrder': 2,
                'offset': -45.0,
                'outerRadius': 44,
                'expandFactor': .2
               }
        protExtractUnitCell = self.newProtocol(XmippProtExtractUnit,
                                               **args)
        self.launchProtocol(protExtractUnitCell)
        unit_cell = protExtractUnitCell.outputVolume
        self.assertTrue(unit_cell.getFileName())

        #dock chainA in map unit cell
        args = {
                'inputVolume1': unit_cell,
                'resolution': 3.2,
                'inputStructure': chainA
                }

        protDockInMap = self.newProtocol(PhenixProtRunDockInMap, **args)
        protDockInMap.setObjLabel('DockInMap\n'
                                  'haemoglobin unit cell\n')
        self.launchProtocol(protDockInMap)
        self.assertTrue(os.path.exists(
            protDockInMap.outputPdb.getFileName()))

    def testDockInMap3(self):
        """ This test checks that phenix dock in map protocol runs with
        a chain of the atomic structure and a whole map"""
        print("Run phenix dock_in_map protocol with a chain of the "
              "imported atomic structure file and an imported whole map")

        # import PDB from Database
        protImportPDB = self._importStructure()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        # import map
        map = self._importVolume()
        self.assertTrue(map.getFileName())

        # extract chain using chimera
        # create auxiliary CMD file for chimera operate
        extraCommands = ""
        extraCommands += "sel #2/A\n"
        extraCommands += "save /tmp/chainA.cif format mmcif models #2 relModel #1 selectedOnly true\n"
        extraCommands += "open /tmp/chainA.cif\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_A_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure1
                }
        protChimera = self.newProtocol(ChimeraProtOperate,
                                       **args)
        protChimera.setObjLabel('chimera operate\n extract chain A')
        self.launchProtocol(protChimera)
        chainA = eval("protChimera.DONOTSAVESESSION_A_Atom_struct__3_%06d" % \
                      protChimera.getObjId())
        self.assertIsNotNone(chainA.getFileName(),
                              "There was a problem with the alignment")

        # dock chainA in the full map
        args = {
                'inputVolume1': map,
                'resolution': 3.2,
                'inputStructure': chainA
               }

        protDockInMap = self.newProtocol(PhenixProtRunDockInMap, **args)
        protDockInMap.setObjLabel('DockInMap\n'
                                  'haemoglobin\n'
                                  'and chain A')
        self.launchProtocol(protDockInMap)
        self.assertTrue(os.path.exists(
            protDockInMap.outputPdb.getFileName()))

    def testDockInMap4(self):
        """ This test checks that phenix dock in map protocol runs with
        two identical chains of the atomic structure and a whole map"""
        print("Run phenix dock_in_map protocol with two identical chains of the "
              "imported atomic structure file and an imported whole map")

        # import PDB from Database
        protImportPDB = self._importStructure()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        # import map
        map = self._importVolume()
        self.assertTrue(map.getFileName())

        # extract chain using chimera
        # create auxiliary CMD file for chimera operate
        extraCommands = ""
        extraCommands += "sel #2/A\n"
        extraCommands += "save /tmp/chainA.cif format mmcif models #2 relModel #1 selectedOnly true\n"
        extraCommands += "open /tmp/chainA.cif\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_A_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure1
                }
        protChimera = self.newProtocol(ChimeraProtOperate,
                                       **args)
        protChimera.setObjLabel('chimera operate\n extract chain A')
        self.launchProtocol(protChimera)
        chainA = eval("protChimera.DONOTSAVESESSION_A_Atom_struct__3_%06d" % \
                      protChimera.getObjId())
        self.assertIsNotNone(chainA.getFileName(),
                             "There was a problem with the alignment")

        # dock two chainA in the full map
        args = {
                'inputVolume1': map,
                'resolution': 3.2,
                'inputStructure': chainA,
                'modelCopies': 2
               }

        protDockInMap = self.newProtocol(PhenixProtRunDockInMap, **args)
        protDockInMap.setObjLabel('DockInMap\n'
                                  'haemoglobin\n'
                                  'and 2 chains A_C')
        self.launchProtocol(protDockInMap)
        self.assertTrue(os.path.exists(
            protDockInMap.outputPdb.getFileName()))