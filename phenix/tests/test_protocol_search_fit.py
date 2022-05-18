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


# protocol to test the phenix search fit between a sequence and a map density
from phenix.protocols import PhenixProtSearchFit
from pwem.protocols.protocol_import import (ProtImportPdb,
                                                    ProtImportVolumes)
from chimera.protocols import ChimeraProtOperate
from xmipp3.protocols.protocol_extract_asymmetric_unit import XmippProtExtractUnit
from pwem.constants import SCIPION_SYM_NAME
from xmipp3.constants import XMIPP_SYM_NAME, XMIPP_TO_SCIPION, XMIPP_I222r
from pyworkflow.tests import *
import pwem.protocols as emprot
import os.path


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import map volumes, atomic structures(PDBx/mmCIF files)
    and sequences of specific chains of the atomic structures.
    """
    pdbID1 = "5ni1"  # Haemoglobin atomic structure
    CHAIN1 = '{"model": 0, "chain": "A", "residues": 141}'
    pdbID2 = "6qi5"  # Atadenovirus atomic structure
    CHAIN2 = '{"model": 0, "chain": "M", "residues": 451}'
    NAME1 = '5ni1_A-seq' # Haemoglobin sequence (chain A)
    NAME2 = '6qi5_M-seq' # Atadenovirus sequence (chain M)
    removeResidues1 = '{"index": "90-120", "residues": "KLRVDPVNFKLLSHCLLVTLAAHLPAEFTPA"}'
    removeResidues1_2 = '{"index": "80-130", "residues": "LSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLA"}'
    removeResidues2 = '{"index": "42-55", "residues": "YFPHFDLSHGSAQV"}'
    removeResidues3 = '{"index": "128-157", "residues": "ELSIPEGDYTVGSLIDMLNNAVVENYLEVG"}'

    def importVolume(self, args, label):
        protImportVol = self.newProtocol(emprot.ProtImportVolumes, **args)
        protImportVol.setObjLabel(label)
        self.launchProtocol(protImportVol)
        return protImportVol.outputVolume

    def _importVolumeHEMOGLOBINA(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/emd_3488.map'),
                'samplingRate': 1.05,
                'setOrigCoord': True,
                'x': 0.0,
                'y': 0.0,
                'z': 0.0
                }
        return self.importVolume(args, 'import volume haemoglobin\n'
                                       ' with default origin\n')

    def _importVolumeATADENOVIRUS(self):
        args = {'importFrom': ProtImportVolumes.IMPORT_FROM_EMDB,
                'emdbId': 4551
                }
        return self.importVolume(args, 'import volume 4551\n'
                                       'atadenovirus\n')

    def importAtomStruct(self, args, label):
        protImportPDB = self.newProtocol(emprot.ProtImportPdb, **args)
        protImportPDB.setObjLabel(label)
        self.launchProtocol(protImportPDB)
        return protImportPDB.outputPdb

    def _importAtomStructHEMOGLOBINA(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': self.pdbID1
                }
        return self.importAtomStruct(args, 'import pdb\n 5ni1')

    def _importAtomStructATADENOVIRUS(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': self.pdbID2
                }
        return self.importAtomStruct(args,'import pdb\n 6qi5')

    def importSequence(self, args, label):
        protImportSequence = self.newProtocol(emprot.ProtImportSequence, **args)
        protImportSequence.setObjLabel(label)
        self.launchProtocol(protImportSequence)
        return protImportSequence.outputSequence

    def _importSequenceHEMOGLOBINA(self):
        args = {'inputSequenceName': self.NAME1,
                'inputProteinSequence':
                    emprot.ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    emprot.ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': self.pdbID1,
                'inputStructureChain': self.CHAIN1
                }
        return self.importSequence(args, 'import sequence\n 5ni1_A_seq')

    def _importSequenceATADENOVIRUS(self):
        args = {'inputSequenceName': self.NAME2,
                'inputProteinSequence':
                    emprot.ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    emprot.ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': self.pdbID2,
                'inputStructureChain': self.CHAIN2
                }
        return self.importSequence(args,'import sequence\n 6qi5_M_seq')

class TestPhenixProtSearchFit(TestImportData):
    """ Test the chimera subtraction map protocol
    """

    def testPhenixSearchFit1(self):
        """ This test checks the fitting of the structure of a traced fragment (residues)
          of the haemoglobin chain A in the whole density map"""
        print("Run Phenix Search Fit to fit the alpha-helix structure "
              "(traced fragment) of the haemoglobin chain A in"
              " the whole density map\n")

        # Import Volume
        volume = self._importVolumeHEMOGLOBINA()

        # import PDB
        structure = self._importAtomStructHEMOGLOBINA()

        # import sequence
        sequence = self._importSequenceHEMOGLOBINA()

        # create auxiliary CMD file for chimera operate to select only
        # a small fragment of the structure (helix between residues 94 and 118)
        extraCommands = ""
        extraCommands += "select all & ~ #2/A:94-118\n"
        extraCommands += "del sel\n"
        extraCommands += "scipionwrite #2 " \
                         "prefix DONOTSAVESESSION_5ni1_chainA_94_118_\n"
        extraCommands += "exit\n"

        args1 = {'extraCommands': extraCommands,
                 'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate, **args1)
        protChimera1.setObjLabel('chimera operate\n fragment_5ni1_A')
        self.launchProtocol(protChimera1)

        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        result = ""
        try:
            result = eval(
                "protChimera1.DONOTSAVESESSION_5ni1_chainA_94_118_Atom_struct__2_%06d" %
                protChimera1.getObjId())
        except (NameError, SyntaxError) as e:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result.getFileName()))

        args2 = {'inputVolume': volume,
                 'resolution': 3.2,
                 'inputStructure': result,
                 'inputSequence': sequence,
                 'residues': self.removeResidues1,
                 'numberOfMpi' : 8
                }
        protSearchFit1 = self.newProtocol(PhenixProtSearchFit, **args2)
        protSearchFit1.setObjLabel('search fit\n volume 3844\n5ni1_A_94_118')
        self.launchProtocol(protSearchFit1)

        self.assertTrue(os.path.exists(protSearchFit1.outputAtomStruct_4.getFileName()))

    def testPhenixSearchFit2(self):
        """ This test checks the fitting of the structure of the traced fragment
            carbon skeleton (ALA) of the haemoglobin chain A in the whole density map"""
        print("Run Phenix Search Fit to fit the alpha-helix structure "
              "(carbon skeleton fragment) of the haemoglobin chain A in "
              "the whole density map to retrieve the residues\n")

        # Import Volume
        volume = self._importVolumeHEMOGLOBINA()

        # import PDB
        structure = self._importAtomStructHEMOGLOBINA()

        # import sequence
        sequence = self._importSequenceHEMOGLOBINA()

        # create auxiliary CMD file for chimera operate to select only
        # a small fragment of the structure (helix between residues 94 and 118)
        # and mutate it to ALA
        extraCommands = ""
        extraCommands += "select all & ~ #2/A:94-118\n"
        extraCommands += "del sel\n"
        extraCommands += "swapaa #2/A:94-118 ALA\n"
        extraCommands += "scipionwrite #2 " \
                         "prefix DONOTSAVESESSION_5ni1_chainA_94_118_MutALA_\n"
        extraCommands += "exit\n"

        args1 = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure
                }
        protChimera2 = self.newProtocol(ChimeraProtOperate, **args1)
        protChimera2.setObjLabel('chimera operate\n fragment_5ni1_A_MutALA')
        self.launchProtocol(protChimera2)

        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        result = ""
        try:
            result = eval(
                "protChimera2.DONOTSAVESESSION_5ni1_chainA_94_118_MutALA_Atom_struct__2_%06d" %
                protChimera2.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result.getFileName()))

        args2 = {'inputVolume': volume,
                 'resolution': 3.2,
                 'inputStructure': result,
                 'inputSequence': sequence,
                 'residues': self.removeResidues1,
                 'numberOfMpi' : 8
                }
        protSearchFit2 = self.newProtocol(PhenixProtSearchFit, **args2)
        protSearchFit2.setObjLabel(
            'search fit\n volume 3844\nfragment_5ni1_A_94_118_MutALA\nseq 5ni1_A')
        self.launchProtocol(protSearchFit2)

        self.assertTrue(os.path.exists(protSearchFit2.outputAtomStruct_4.getFileName()))

        # Branch to observe the effect of modifying the residue overlapping between
        # the sequence and the ALA chain (firstaa and lastaa)
        args3 = {'inputVolume': volume,
                 'resolution': 3.2,
                 'inputStructure': result,
                 'inputSequence': sequence,
                 'residues': self.removeResidues1_2,
                 'numberOfMpi' : 8
                 }
        protSearchFit3 = self.newProtocol(PhenixProtSearchFit, **args3)
        protSearchFit3.setObjLabel(
            'search fit\n volume 3844\nfragment_5ni1_A_94_118_MutALA\noverlap_80_130')
        self.launchProtocol(protSearchFit3)

        self.assertTrue(os.path.exists(protSearchFit3.outputAtomStruct_4.getFileName()))

    def testPhenixSearchFit3(self):
        """ This test checks the fitting of the structure of a traced fragment (residues)
          of the haemoglobin chain A in the whole density map"""
        print("Run Phenix Search Fit to fit the loop like structure "
              "(traced fragment) of the haemoglobin chain A in"
              " the whole density map\n")

        # Import Volume
        volume = self._importVolumeHEMOGLOBINA()

        # import PDB
        structure = self._importAtomStructHEMOGLOBINA()

        # import sequence
        sequence = self._importSequenceHEMOGLOBINA()

        # create auxiliary CMD file for chimera operate to select only
        # a small fragment of the structure (loop between residues 42 and 55)
        extraCommands = ""
        extraCommands += "select all & ~ #2/A:42-55\n"
        extraCommands += "del sel\n"
        extraCommands += "scipionwrite #2 " \
                         "prefix DONOTSAVESESSION_5ni1_chainA_42_55_\n"
        extraCommands += "exit\n"

        args1 = {'extraCommands': extraCommands,
                 'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate, **args1)
        protChimera1.setObjLabel('chimera operate\n fragment_5ni1_A')
        self.launchProtocol(protChimera1)

        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        result = ""
        try:
            result = eval(
                "protChimera1.DONOTSAVESESSION_5ni1_chainA_42_55_Atom_struct__2_%06d" %
                protChimera1.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result.getFileName()))

        args2 = {'inputVolume': volume,
                 'resolution': 3.2,
                 'inputStructure': result,
                 'inputSequence': sequence,
                 'residues': self.removeResidues2,
                 'numberOfMpi' : 8
                }
        protSearchFit3 = self.newProtocol(PhenixProtSearchFit, **args2)
        protSearchFit3.setObjLabel('search fit\n volume 3844\n5ni1_A_42_55')
        self.launchProtocol(protSearchFit3)

        self.assertTrue(os.path.exists(protSearchFit3.outputAtomStruct_4.getFileName()))

    def testPhenixSearchFit4(self):
        """ This test checks the fitting of the structure of the traced fragment
            carbon skeleton (ALA) of the haemoglobin chain A in the whole density map"""
        print("Run Phenix Search Fit to fit the loop like "
              "(carbon skeleton fragment) of the haemoglobin chain A in "
              "the whole density map to retrieve the residues\n")

        # Import Volume
        volume = self._importVolumeHEMOGLOBINA()

        # import PDB
        structure = self._importAtomStructHEMOGLOBINA()

        # import sequence
        sequence = self._importSequenceHEMOGLOBINA()

        # create auxiliary CMD file for chimera operate to select only
        # a small fragment of the structure (loop between residues 42 and 55)
        extraCommands = ""
        extraCommands += "select all & ~ #2/A:42-55\n"
        extraCommands += "del sel\n"
        extraCommands += "swapaa #2/A:42-55 ALA\n"
        extraCommands += "scipionwrite #2 " \
                         "prefix DONOTSAVESESSION_5ni1_chainA_42_55_MutALA_\n"
        extraCommands += "exit\n"

        args1 = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure
                }
        protChimera2 = self.newProtocol(ChimeraProtOperate, **args1)
        protChimera2.setObjLabel('chimera operate\n fragment_5ni1_A_MutALA')
        self.launchProtocol(protChimera2)

        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        result = ""
        try:
            result = eval(
                "protChimera2.DONOTSAVESESSION_5ni1_chainA_42_55_MutALA_Atom_struct__2_%06d" %
                protChimera2.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result.getFileName()))

        args2 = {'inputVolume': volume,
                 'resolution': 3.2,
                 'inputStructure': result,
                 'inputSequence': sequence,
                 'residues': self.removeResidues2,
                 'numberOfMpi' : 8
                }
        protSearchFit4 = self.newProtocol(PhenixProtSearchFit, **args2)
        protSearchFit4.setObjLabel(
            'search fit\n volume 3844\nfragment_5ni1_A_42_55_MutALA\nseq 5ni1_A')
        self.launchProtocol(protSearchFit4)

        self.assertTrue(os.path.exists(protSearchFit4.outputAtomStruct_4.getFileName()))

    def testPhenixSearchFit5(self):
        """ This test checks the fitting of the structure of a traced fragment (residues)
          of the adenovirus chain M in the whole density map"""
        print("Run Phenix Search Fit to fit the loop like structure "
              "(traced fragment) of the adenovirus chain M in"
              " the whole density map\n")

        # Import Volume
        volume = self._importVolumeATADENOVIRUS()

        # import PDB
        structure = self._importAtomStructATADENOVIRUS()

        # import sequence
        sequence = self._importSequenceATADENOVIRUS()

        # run extract asymmetric unit
        # sym = XMIPP_SYM_NAME[XMIPP_I222r] + " " \
        #       "(" + SCIPION_SYM_NAME[XMIPP_TO_SCIPION[XMIPP_I222r]] + ")"
        sym = XMIPP_I222r
        args = {'inputVolumes': volume,
                'symmetryGroup': sym,
                'innerRadius': 202.0,
                'outerRadius': 389.0,
                'expandFactor': .2}

        protExtractAsymUnit = self.newProtocol(XmippProtExtractUnit, **args)
        protExtractAsymUnit.setObjLabel('extract asym unit')
        self.launchProtocol(protExtractAsymUnit)

        # create auxiliary CMD file for chimera operate to select only
        # a small fragment of the structure (helix between residues 130 and 156)
        extraCommands = ""
        extraCommands += "select #3 & ~ #3/M:130-156\n"
        extraCommands += "del sel\n"
        extraCommands += "sym #3 i,222r copies true\n"
        extraCommands += "move -462.88,16.37,26.51 coordinateSystem #1 models #3\n"
        extraCommands += "mmaker #3 to #4.37\n"
        extraCommands += "fitmap #3 inMap #2\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_6qi5_chainM_130_156_\n"
        extraCommands += "exit\n"

        args1 = {'extraCommands': extraCommands,
                 'inputVolume': protExtractAsymUnit.outputVolume,
                 'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate, **args1)
        protChimera1.setObjLabel('chimera operate\n fragment_6qi5_M')
        self.launchProtocol(protChimera1)

        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        result = ""
        try:
            result = eval(
                "protChimera1.DONOTSAVESESSION_6qi5_chainM_130_156_Atom_struct__3_%06d" %
                protChimera1.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result.getFileName()))

        args2 = {'inputVolume': protExtractAsymUnit.outputVolume,
                 'resolution': 3.4,
                 'inputStructure': result,
                 'inputSequence': sequence,
                 'residues': self.removeResidues3,
                 'numberOfMpi' : 8
                }
        protSearchFit5 = self.newProtocol(PhenixProtSearchFit, **args2)
        protSearchFit5.setObjLabel('search fit\n volume 4551\n5ni1_M_130_156')
        self.launchProtocol(protSearchFit5)

        self.assertTrue(os.path.exists(protSearchFit5.outputAtomStruct_4.getFileName()))


