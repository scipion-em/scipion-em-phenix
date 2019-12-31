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

# protocol to test the phenix protocol superpose_pdbs
import os

from phenix.protocols import PhenixProtRunDockInMap
from pwem.protocols.protocol_import import (ProtImportPdb,
                                            ProtImportVolumes)
from pyworkflow.tests import *


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

    def _importUnitCell(self):
        args = {'filesPath': self.dsModBuild.getFile(
                'volumes/hemoglobin_unit_cell.mrc'),
                'samplingRate': 1.05,
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_3488 unit cell\n')
        self.launchProtocol(protImportVol)
        volume2 = protImportVol.outputVolume
        return volume2

    def _importStructure(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': self.pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1')
        self.launchProtocol(protImportPDB)
        structure1_PDB = protImportPDB.outputPdb
        return structure1_PDB

    def _importStructureChainA(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/hemoglobin_chainA.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1 chainA')
        self.launchProtocol(protImportPDB)
        chainA = protImportPDB.outputPdb
        return chainA


class TestProtDockInMap(TestImportData):
    """ Test the protocol Dock In Map
    """

    def testDockInMap1(self):
        """ This test checks that phenix dock in map protocol runs with
        an atomic structure and a map"""
        print("Run phenix dock_in_map protocol from an imported atomic "
              "structure file and an imported map")

        # import PDB
        structure1 = self._importStructure()
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
        # import PDB
        chainA = self._importStructureChainA
        self.assertTrue(chainA.getFileName())
        self.assertFalse(chainA.getVolume())

        # import map
        unit_cell = self._importUnitCell
        self.assertTrue(unit_cell.getFileName())

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

        # import PDB
        chainA = self._importStructureChainA
        self.assertTrue(chainA.getFileName())
        self.assertFalse(chainA.getVolume())

        # import map
        map = self._importVolume()
        self.assertTrue(map.getFileName())

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

        # import PDB
        chainA = self._importStructureChainA
        self.assertTrue(chainA.getFileName())
        self.assertFalse(chainA.getVolume())

        # import map
        map = self._importVolume()
        self.assertTrue(map.getFileName())

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