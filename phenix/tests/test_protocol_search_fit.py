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


# protocol to test the phenix search fit between a sequence and a map density TODO
from phenix.protocols import PhenixProtSearchFit
from pwem.protocols.protocol_import import (ProtImportPdb,
                                                    ProtImportVolumes)
from chimera.protocols import ChimeraProtOperate
from pyworkflow.tests import *
import pwem.protocols as emprot
import os.path


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import map volumes and atomic structures(PDBx/mmCIF files)
    """
    pdbID1 = "5ni1"  # Haemoglobin atomic structure
    NAME1 = '5ni1_A-seq'
    CHAIN1 = '{"model": 0, "chain": "A", "residues": 141}'
    FirstResidue = '{"residue": 90, "K"}'
    LastResidue = '{"residue": 120, "A"}'
    def _importVolume(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/emd_3488.map'),
                'samplingRate': 1.05,
                'setOrigCoord': True,
                'x': 0.0,
                'y': 0.0,
                'z': 0.0
                }
        protImportVol = self.newProtocol(emprot.ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume haemoglobin\n with default '
                                  'origin\n')
        self.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        return volume

    def _importAtomStruct(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': self.pdbID1
                }
        protImportPDB = self.newProtocol(emprot.ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1')
        self.launchProtocol(protImportPDB)
        structure = protImportPDB.outputPdb
        return structure

    def _importSequence(self):
        args = {'inputSequenceName': self.NAME1,
                'inputProteinSequence':
                    emprot.ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    emprot.ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': self.pdbID1,
                'inputStructureChain': self.CHAIN1
                }
        protImportSequence = self.newProtocol(emprot.ProtImportSequence, **args)
        protImportSequence.setObjLabel('import sequence\n 5ni1_A_seq')
        self.launchProtocol(protImportSequence)
        sequence = protImportSequence.outputSequence
        return sequence


class TestPhenixProtSearchFit(TestImportData):
    """ Test the chimera subtraction map protocol
    """

    def testPhenixSearchFit1(self):
        """ This test checks  """
        print("Run Chimera subtraction of a derived-model map "
              "from the imported volume\n")

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure = self._importAtomStruct()

        # import sequence
        sequence = self._importSequence()

        # create auxiliary CMD file for chimera operate to select only
        # a small fragment of the structure (helix between residues 94 and 118)
        extraCommands = ""
        extraCommands += "select all & ~ #3/A:94-118\n"
        extraCommands += "del sel\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_5ni1_chainA_94_118_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate,
                                       **args)
        protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera1)
        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        try:
            result = eval(
                "protChimera1.DONOTSAVESESSION_5ni1_chainA_94_118_Atom_struct__3_%06d.getFileName()"
                          % protChimera1.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        args = {'inputVolume': volume,
                'resolution' : 3.2,
                'inputStructure': result,
                'inputSequence': sequence,
                'firstaa': self.FirstResidue,
                'lastaa' : self.LastResidue
                }
        protSearchFit1 = self.newProtocol(PhenixProtSearchFit,
                                        **args)
        protSearchFit1.setObjLabel('search fit\n volume 3844\npdb 5ni1_A_94_118\nseq 5ni1_A')
        self.launchProtocol(protSearchFit1)

    def testPhenixSearchFit2(self):
        """ This test checks  """
        print("Run Chimera subtraction of a derived-model map "
              "from the imported volume\n")

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure = self._importAtomStruct()

        # import sequence
        sequence = self._importSequence()

        # create auxiliary CMD file for chimera operate to select only
        # a small fragment of the structure (helix between residues 94 and 118)
        extraCommands = ""
        extraCommands += "select all & ~ #3/A:94-118\n"
        extraCommands += "del sel\n"
        extraCommands += "swapaa #3/A:94-118 ALA\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_5ni1_chainA_94_118_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera2 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera2.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera2)
        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        try:
            result = eval(
                "protChimera2.DONOTSAVESESSION_5ni1_chainA_94_118_MutALA_Atom_struct__3_%06d.getFileName()"
                % protChimera2.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        args = {'inputVolume': volume,
                'resolution': 3.2,
                'inputStructure': result,
                'inputSequence': sequence,
                'firstaa': self.FirstResidue,
                'lastaa': self.LastResidue
                }
        protSearchFit2 = self.newProtocol(PhenixProtSearchFit,
                                          **args)
        protSearchFit2.setObjLabel('search fit\n volume 3844\npdb 5ni1_A_94_118_MutALA\nseq 5ni1_A')
        self.launchProtocol(protSearchFit2)


        # # TODO: These steps of protChimera2 can not be performed in protChimera1 because
        # # when the map is saved keep an inappropriate origin
        # extraCommands = ""
        # # extraCommands += "molmap #2 2.1 gridSpacing 1.05 modelId 3\n"
        # extraCommands += "molmap #2 2.1 gridSpacing 1.05 replace false\n"
        # extraCommands += "scipionwrite #3 " \
        #                  "prefix DONOTSAVESESSION_\n"
        # extraCommands += "exit\n"
        # result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
        #               % protChimera1.getObjId())
        # args = {'extraCommands': extraCommands,
        #         'pdbFileToBeRefined': result,
        #         }
        # protChimera2 = self.newProtocol(ChimeraProtOperate,
        #                                 **args)
        # protChimera2.setObjLabel('chimera operate\n pdb\n save '
        #                          'model-derived map')
        # self.launchProtocol(protChimera2)
        # try:
        #     result = eval("protChimera2.DONOTSAVESESSION_Map__3_%06d.getFileName()"
        #                   % protChimera2.getObjId())
        # except:
        #     self.assertTrue(False, "There was a problem with the alignment")
        #
        # self.assertTrue(os.path.exists(result))
        #
        # # protocol chimera map subtraction
        # extraCommands = "run(session, 'select all')\n"
        # extraCommands += "run(session, 'exit')\n"
        # result = eval("protChimera2.DONOTSAVESESSION_Map__3_%06d"
        #               % protChimera2.getObjId())
        # args = {'extraCommands': extraCommands,
        #         'inputVolume': volume,
        #         'mapOrModel': 0,
        #         'inputVolume2': result,
        #         }
        # protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
        # protChimera3.setObjLabel('chimera subtract\n map -\n'
        #                          'model-derived map\n')
        # self.launchProtocol(protChimera3)
        #
        # # Dynamically defined name of the variable because it does depend on
        # # the protocol ID
        # try:
        #     result = eval("protChimera3.difference_Map__8_%06d.getFileName()"
        #          % protChimera3.getObjId())
        # except:
        #     self.assertTrue(False,  "There was a problem with the alignment")
        #
        # self.assertTrue(os.path.exists(result))
        #
        # try:
        #     result = eval("protChimera3.filtered_Map__9_%06d.getFileName()"
        #          % protChimera3.getObjId())
        # except:
        #     self.assertTrue(False,  "There was a problem with the alignment")
        #
        # self.assertTrue(os.path.exists(result))