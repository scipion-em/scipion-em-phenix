# ***************************************************************************
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

from pwem.protocols.protocol_import import (ProtImportSetOfAtomStructs, ProtImportVolumes)
from phenix.protocols import PhenixProtRunRSRefineSet
from pyworkflow.tests import *


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')

    @classmethod
    def _importVolume(cls):
        args = {'filesPath': cls.dsModBuild.getFile('volumes/1ake_4-5A.mrc'),
                'samplingRate': 1.5,
                'setOrigCoord': False
                }
        protImportVol = cls.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume 1ake_4-5A\n with default '
                                  'origin\n')
        cls.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        return volume

    @classmethod
    def _importStructurePDBSet(cls):
      args = {'inputPdbData': ProtImportSetOfAtomStructs.IMPORT_FROM_FILES,
              'filesPath': cls.dsModBuild.getFile('PDBx_mmCIF'),
              'filesPattern': cls.dsModBuild.getFile('PDBx_mmCIF/1ake_*.pdb')
              }
      protImportPDBs = cls.newProtocol(ProtImportSetOfAtomStructs, **args)
      protImportPDBs.setObjLabel('import 1akes structures\n')
      cls.launchProtocol(protImportPDBs)

      return getattr(protImportPDBs, protImportPDBs._OUTNAME)


class TestRSRefineSet(TestImportBase):

  def testRSRefineSet(self):
    """"""
    args = {'inputVolume': self._importVolume(),
            'resolution': 3.5,
            'inputStructureSet': self._importStructurePDBSet(),
            }

    protMolProb = self.newProtocol(PhenixProtRunRSRefineSet, **args)
    protMolProb.setObjLabel('Real space refine on SetOfAtomStruct')
    self.launchProtocol(protMolProb)
    self.assertEqual(len(protMolProb.outputAtomStructs), len(protMolProb.inputStructureSet.get()))