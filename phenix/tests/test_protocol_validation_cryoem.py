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

# protocol to test the phenix method real_spacerefine

from pyworkflow.em.protocol.protocol_import import (ProtImportPdb,
                                                    ProtImportVolumes)
from phenix.protocols.protocol_real_space_refine import (PhenixProtRunRSRefine,
                                                         mmCIF)
from phenix.protocols.protocol_molprobity import PhenixProtRunMolprobity
from phenix.protocols.protocol_validation_cryoem import PhenixProtRunValidationCryoEM
from pyworkflow.tests import *


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import map volumes and atomic structures(PDBx/mmCIF files)
    """

    def _importVolRefmac3(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/refmac3.mrc'),
            'samplingRate': 1.5,
            'setOrigCoord': True,
            'x': 37.5,
            'y': 37.5,
            'z': 37.5
        }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume refmac3.mrc\nset '
                                  'origin 37.5 37.5 37.5\n')
        self.launchProtocol(protImportVol)
        volume_refmac3 = protImportVol.outputVolume
        return volume_refmac3

    def _importVolHemoOrg(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/emd_3488.map'),
            'samplingRate': 1.05,
            'setOrigCoord': True,
            'x': 0.0,
            'y': 0.0,
            'z': 0.0
        }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_3488.map\nset '
                                  'origin 0.0 0.0 0.0\n')
        self.launchProtocol(protImportVol)
        volume_hemo_orig = protImportVol.outputVolume
        return volume_hemo_orig

    def _importVolNucleosome(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/emd_6838.mrc'),
            'samplingRate': 1.4,
            'setOrigCoord': True,
            'x': 0.0,
            'y': 0.0,
            'z': 0.0
        }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_6838.mrc\nset '
                                  'origin 0.0 0.0 0.0\n')
        self.launchProtocol(protImportVol)
        volume_nucleosome = protImportVol.outputVolume
        return volume_nucleosome

    def _importHemoHalf1Half2(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/emd_3488.map'),
                'samplingRate': 1.05,
                'setOrigCoord': True,
                'x': 0.0,
                'y': 0.0,
                'z': 0.0,
                'setHalfMaps': True,
                'half1map': self.dsModBuild.getFile(
                            'volumes/emd_3488_Noisy_half1.map'),
                'half2map': self.dsModBuild.getFile(
                            'volumes/emd_3488_Noisy_half2.map'),
        }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume hemoglobin\nset '
                                  'origin 0.0 0.0 0.0\nhalf1 and half2')
        self.launchProtocol(protImportVol)
        volume_hemo_half1_hal2 = protImportVol.outputVolume
        return volume_hemo_half1_hal2

    def _importStructRefmac3(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/refmac3.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n refmac3.pdb')
        self.launchProtocol(protImportPDB)
        structure_refmac3 = protImportPDB.outputPdb
        return structure_refmac3

    def _importStructHemoPDB(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/5ni1.pdb'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1.pdb')
        self.launchProtocol(protImportPDB)
        structure_hemo_pdb = protImportPDB.outputPdb
        return structure_hemo_pdb

    def _importStructNucleosomePDB(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/chimera_3lz0.pdb'),
                'inputVolume': self._importVolNucleosome()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import cif\n nucleosome\nassociated volume '
                                  'emd_6838.mrc')
        self.launchProtocol(protImportPDB)
        structure_nucleosome_pdb = protImportPDB.outputPdb
        self.assertTrue(structure_nucleosome_pdb.getFileName())
        return structure_nucleosome_pdb

class TestValCryoEM(TestImportData):
    """ Test the protocol of Real Space Refine
    """
    def checkMPResults(self, ramOutliers, ramFavored, rotOutliers,
                       cbetaOutliers, clashScore, overallScore,
                       protMolProbity, places=2):
        # method to check MolProbity statistic results of the Final Results
        # Table
        self.assertAlmostEqual(protMolProbity.ramachandranOutliers,
                               ramOutliers, places)
        self.assertAlmostEqual(protMolProbity.ramachandranFavored,
                               ramFavored, places)
        self.assertAlmostEqual(protMolProbity.rotamerOutliers,
                               rotOutliers, places)
        self.assertAlmostEqual(protMolProbity.cbetaOutliers,
                               cbetaOutliers, places)
        self.assertAlmostEqual(protMolProbity.clashscore,
                               clashScore, places)
        self.assertAlmostEqual(protMolProbity.overallScore,
                               overallScore, places)

    def checkRSRefineResults(self, ramOutliers, ramFavored, rotOutliers,
                             cbetaOutliers, clashScore, overallScore,
                             protRSRefine, places=2):
        # method to check Real Space Refine statistic results of the Final Results
        # Table
        self.assertAlmostEqual(protRSRefine.ramachandranOutliers.get(),
                               ramOutliers, places)
        self.assertAlmostEqual(protRSRefine.ramachandranFavored.get(),
                               ramFavored, delta=1)
        self.assertAlmostEqual(protRSRefine.rotamerOutliers.get(),
                               rotOutliers, delta=3)
        self.assertAlmostEqual(protRSRefine.cbetaOutliers.get(),
                               cbetaOutliers, places)
        self.assertAlmostEqual(protRSRefine.clashscore.get(),
                               clashScore, delta=0.5)
        self.assertAlmostEqual(protRSRefine.overallScore.get(),
                               overallScore, delta=0.5)

    def checkValCryoEMResults(self, ramOutliers, ramFavored, rotOutliers,
                             cbetaOutliers, clashScore, overallScore,
                             protValCryoEM, places=2):
        # method to check Validation CryoEM statistic results of the Final Results
        # Table
        self.assertAlmostEqual(protValCryoEM.ramachandranOutliers.get(),
                               ramOutliers, places)
        self.assertAlmostEqual(protValCryoEM.ramachandranFavored.get(),
                               ramFavored, delta=1)
        self.assertAlmostEqual(protValCryoEM.rotamerOutliers.get(),
                               rotOutliers, delta=3)
        self.assertAlmostEqual(protValCryoEM.cbetaOutliers.get(),
                               cbetaOutliers, places)
        self.assertAlmostEqual(protValCryoEM.clashscore.get(),
                               clashScore, delta=0.5)
        self.assertAlmostEqual(protValCryoEM.overallScore.get(),
                               overallScore, delta=0.5)

    def testValCryoEMFromPDB(self):
        """ This test checks that phenix real space refine protocol runs
        with an atomic structure; No Volume was provided and an error message
        is expected
        """
        print ("Run phenix real space refine protocol from imported pdb file "
               "without imported or pdb-associated volume")

        # import PDB
        structure_refmac3 = self._importStructRefmac3()
        self.assertTrue(structure_refmac3.getFileName())
        self.assertFalse(structure_refmac3.getVolume())

        # validation_cryoem
        args = {
                'resolution': 3.5,
                'inputStructure': structure_refmac3
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM without\nvolume, should NOT work')
        try:
            self.launchProtocol(protValCryoEM)
        except Exception as e:
            self.assertTrue(True, "This test should return a error message as '"
                  " Input volume cannot be EMPTY.'\n")

            return
        self.assertTrue(False)

    def testValCryoEMFromVolumeAndPDB1(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol, the input PDB was fitted
        to the volume and refined previously by coot and refmac withouth mask
        in another project; (MolProbity has been run to compare values before
        and after refinement); default refine strategy;
        MolProbity and ValCryoEM after RSRefine should show identical values
        """
        print ("Run phenix real_space_refine from imported volume and pdb file "
               "previously fitted and refined by Coot and Refmac without mask "
               "(MolProbity has been run to compare values before and after "
               "refinement); default refine strategy; MolProbity and ValCryoEM "
               "after RSRefine should show identical values")

        # Import Volume
        volume_refmac3 = self._importVolRefmac3()

        # import PDB
        structure_refmac3 = self._importStructRefmac3()

        #MolProbity
        args = {
                'inputVolume': volume_refmac3,
                'resolution': 3.5,
                'inputStructure': structure_refmac3,
                'numberOfThreads': 4
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'volume and pdb\n')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkMPResults(ramOutliers=0.47,
                            ramFavored=83.96,
                            rotOutliers=5.68,
                            cbetaOutliers=1,
                            clashScore=4.77,
                            overallScore=2.50,
                            protMolProbity=protMolProbity)

        # real_space_refine
        args = {'inputVolume': volume_refmac3,
                'resolution': 3.5,
                'inputStructure': structure_refmac3
                # default parameters in Optimization strategy options
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine\n refmac3.mrc and '
                                   'refmac3.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        # self.checkRSRefineResults(ramOutliers=0.00,
        #                           ramFavored=95.75,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=2.09,
        #                           overallScore=1.27,
        #                           protRSRefine=protRSRefine)
        self.checkRSRefineResults(ramOutliers=0.00,
                                  ramFavored=96.70,
                                  rotOutliers=3.98,
                                  cbetaOutliers=0,
                                  clashScore=4.47,
                                  overallScore=1.89,
                                  protRSRefine=protRSRefine)
        # MolProbity2
        args = {
            'inputVolume': volume_refmac3,
            'resolution': 3.5,
            'inputStructure': protRSRefine.outputPdb,
            'numberOfThreads': 4
        }
        protMolProbity2 = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity2.setObjLabel('MolProbity\n'
                                   'after RSRefine\n')
        self.launchProtocol(protMolProbity2)

        # check MolProbity results
        self.checkMPResults(ramOutliers=0.00,
                            ramFavored=96.70,
                            rotOutliers=3.98,
                            cbetaOutliers=0,
                            clashScore=4.47,
                            overallScore=1.89,
                            protMolProbity=protMolProbity2)

        # validation_cryoEM
        args = {'inputVolume': volume_refmac3,
                'resolution': 3.5,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nrefmac3.mrc and '
                                 'protRSRefine.outputPdb\n')
        self.launchProtocol(protValCryoEM)

        # check validation cryoem results
        # self.checkValCryoEMResults(ramOutliers=0.00,
        #                           ramFavored=95.75,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=2.09,
        #                           overallScore=1.27,
        #                           protValCryoEM=protValCryoEM)
        self.checkValCryoEMResults(ramOutliers=0.00,
                                  ramFavored=96.70,
                                  rotOutliers=3.98,
                                  cbetaOutliers=0,
                                  clashScore=4.47,
                                  overallScore=1.89,
                                  protValCryoEM=protValCryoEM)

    def testValCryoEMFFromVolumeAndPDB2(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks; default refine strategy; (MolProbity has been run to
        compare values before and after refinement). MolProbity and
        ValCryoEM after RSRefine should show identical values.
        """
        print("Run phenix real_space_refine from imported volume and pdb file "
              "from data banks (vol origin 0.0, 0.0, 0.0); default refine "
              "strategy; (MolProbity has been run to compare values before "
              "and after refinement). MolProbity and ValCryoEM after RSRefine "
              "should show identical values.")

        # Import Volume
        volume_hemo_org = self._importVolHemoOrg()

        # import PDB
        structure_hemo_pdb = self._importStructHemoPDB()

        # MolProbity
        args = {
                'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': structure_hemo_pdb,
                'numberOfThreads': 4
        }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'volume and pdb\nhemoglobin')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkMPResults(ramOutliers=0.00,
                            ramFavored=95.23,
                            rotOutliers=0.43,
                            cbetaOutliers=0,
                            clashScore=3.53,
                            overallScore=1.48,
                            protMolProbity=protMolProbity)

        # real_space_refine
        args = {'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': structure_hemo_pdb,
                'numberOfThreads': 4
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        # self.checkRSRefineResults(ramOutliers=0.00,
        #                           ramFavored=98.23,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=1.99,
        #                           overallScore=0.97,
        #                           protRSRefine=protRSRefine)

        self.checkRSRefineResults(ramOutliers=0.00,
                                  ramFavored=97.70,
                                  rotOutliers=2.82,
                                  cbetaOutliers=0,
                                  clashScore=4.86,
                                  overallScore=1.66,
                                  protRSRefine=protRSRefine)

        # MolProbity2
        args = {
            'inputVolume': volume_hemo_org,
            'resolution': 3.2,
            'inputStructure': protRSRefine.outputPdb,
            'numberOfThreads': 4
        }
        protMolProbity2 = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity2.setObjLabel('MolProbity validation\n'
                                   'volume and pdb\nhemoglobin')
        self.launchProtocol(protMolProbity2)

        # check MolProbity results
        self.checkMPResults(ramOutliers=0.00,
                            ramFavored=97.70,
                            rotOutliers=2.82,
                            cbetaOutliers=0,
                            clashScore=4.97,
                            overallScore=1.67,
                            protMolProbity=protMolProbity2)

        # validation_cryoEM
        args = {'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nvolume_hemo_org and '
                                  'protRSRefine.outputPdb\n')
        self.launchProtocol(protValCryoEM)

        # check validation cryoem results
        # self.checkValCryoEMResults(ramOutliers=0.00,
        #                           ramFavored=98.23,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=1.99,
        #                           overallScore=0.97,
        #                           protValCryoEM=protValCryoEM)

        self.checkValCryoEMResults(ramOutliers=0.00,
                                  ramFavored=97.70,
                                  rotOutliers=2.82,
                                  cbetaOutliers=0,
                                  clashScore=4.86,
                                  overallScore=1.66,
                                  protValCryoEM=protValCryoEM)

    def testValCryoEMFFromVolumeAssociatedToPDB3(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks; default refine strategy; (MolProbity has been run to
        compare values before and after refinement). MolProbity and
        ValCryoEM after RSRefine should show identical values.
        """
        print("Run phenix real_space_refine from imported volume and pdb file "
              "from data banks (vol origin 0.0, 0.0, 0.0); default refine "
              "strategy; (MolProbity has been run to compare values before "
              "and after refinement). MolProbity and ValCryoEM after RSRefine "
              "should show identical values.")

        # import PDB
        structure_nucleosome_pdb = self._importStructNucleosomePDB()

        # MolProbity
        args = {
                'inputStructure': structure_nucleosome_pdb,
                'numberOfThreads': 4
        }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'volume and pdb\nnucleosome')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkMPResults(ramOutliers=1.77,
                            ramFavored=87.62,
                            rotOutliers=10.40,
                            cbetaOutliers=2,
                            clashScore=12.17,
                            overallScore=2.98,
                            protMolProbity=protMolProbity)

        # real_space_refine
        args = {
                'resolution': 4.0,
                'inputStructure': structure_nucleosome_pdb,
                'numberOfThreads': 4
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine nucleosome\nvolume and '
                                 'pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        # self.checkRSRefineResults(ramOutliers=0.00,
        #                           ramFavored=98.23,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=1.99,
        #                           overallScore=0.97,
        #                           protRSRefine=protRSRefine)

        self.checkRSRefineResults(ramOutliers=0.00,
                                  ramFavored=94.15,
                                  rotOutliers=10.08,
                                  cbetaOutliers=0,
                                  clashScore=14.97,
                                  overallScore=2.84,
                                  protRSRefine=protRSRefine)

        # MolProbity2
        args = {
            'inputStructure': protRSRefine.outputPdb,
            'numberOfThreads': 4
        }
        protMolProbity2 = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity2.setObjLabel('MolProbity validation\n'
                                   'volume and pdb\nnucleosome')
        self.launchProtocol(protMolProbity2)

        # check MolProbity results
        self.checkMPResults(ramOutliers=0.00,
                            ramFavored=94.15,
                            rotOutliers=10.08,
                            cbetaOutliers=0,
                            clashScore=15.02,
                            overallScore=2.84,
                            protMolProbity=protMolProbity2)

        # validation_cryoEM
        args = {
                'resolution': 4.0,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nnucleosome and '
                                  'protRSRefine.outputPdb\n')
        self.launchProtocol(protValCryoEM)

        # check validation cryoem results
        # self.checkValCryoEMResults(ramOutliers=0.00,
        #                           ramFavored=98.23,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=1.99,
        #                           overallScore=0.97,
        #                           protValCryoEM=protValCryoEM)

        self.checkValCryoEMResults(ramOutliers=0.00,
                                  ramFavored=94.15,
                                  rotOutliers=10.08,
                                  cbetaOutliers=0,
                                  clashScore=14.97,
                                  overallScore=2.84,
                                  protValCryoEM=protValCryoEM)

    def testValCryoEMFFromVolumeAndPDB4(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol (with half1 and half2)
        and the input PDB from data banks; default refine strategy;
        (MolProbity has been run to
        compare values before and after refinement). MolProbity and
        ValCryoEM after RSRefine should show identical values.
        """
        print("Run phenix real_space_refine from imported volume and pdb file "
              "from data banks (vol origin 0.0, 0.0, 0.0); default refine "
              "strategy; (MolProbity has been run to compare values before "
              "and after refinement). MolProbity and ValCryoEM after RSRefine "
              "should show identical values.")

        # Import Volume
        volume_hemo_half1_hal2 = self._importHemoHalf1Half2()

        # import PDB
        structure_hemo_pdb = self._importStructHemoPDB()

        # MolProbity
        args = {
                'inputVolume': volume_hemo_half1_hal2,
                'resolution': 3.2,
                'inputStructure': structure_hemo_pdb,
                'numberOfThreads': 4
        }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'full and half volumes\nand pdb\nhemoglobin')
        self.launchProtocol(protMolProbity)

        # check MolProbity results
        self.checkMPResults(ramOutliers=0.00,
                            ramFavored=95.23,
                            rotOutliers=0.43,
                            cbetaOutliers=0,
                            clashScore=3.53,
                            overallScore=1.48,
                            protMolProbity=protMolProbity)

        # real_space_refine
        args = {'inputVolume': volume_hemo_half1_hal2,
                'resolution': 3.2,
                'inputStructure': structure_hemo_pdb,
                'numberOfThreads': 4
                }
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map (full, half1, half2)\nand '
                                 '5ni1.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        # self.checkRSRefineResults(ramOutliers=0.00,
        #                           ramFavored=98.23,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=1.99,
        #                           overallScore=0.97,
        #                           protRSRefine=protRSRefine)

        self.checkRSRefineResults(ramOutliers=0.00,
                                  ramFavored=97.70,
                                  rotOutliers=2.82,
                                  cbetaOutliers=0,
                                  clashScore=4.86,
                                  overallScore=1.66,
                                  protRSRefine=protRSRefine)

        # MolProbity2
        args = {
            'inputVolume': volume_hemo_half1_hal2,
            'resolution': 3.2,
            'inputStructure': protRSRefine.outputPdb,
            'numberOfThreads': 4
        }
        protMolProbity2 = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity2.setObjLabel('MolProbity validation\n'
                                   'full and half volumes\nand pdb\nhemoglobin')
        self.launchProtocol(protMolProbity2)

        # check MolProbity results
        self.checkMPResults(ramOutliers=0.00,
                            ramFavored=97.70,
                            rotOutliers=2.82,
                            cbetaOutliers=0,
                            clashScore=4.97,
                            overallScore=1.67,
                            protMolProbity=protMolProbity2)

        # validation_cryoEM
        args = {'inputVolume': volume_hemo_half1_hal2,
                'resolution': 3.2,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nvolume_hemo_org and '
                                  'protRSRefine.outputPdb\n')
        self.launchProtocol(protValCryoEM)

        # check validation cryoem results
        # self.checkValCryoEMResults(ramOutliers=0.00,
        #                           ramFavored=98.23,
        #                           rotOutliers=0.00,
        #                           cbetaOutliers=0,
        #                           clashScore=1.99,
        #                           overallScore=0.97,
        #                           protValCryoEM=protValCryoEM)

        self.checkValCryoEMResults(ramOutliers=0.00,
                                  ramFavored=97.70,
                                  rotOutliers=2.82,
                                  cbetaOutliers=0,
                                  clashScore=4.86,
                                  overallScore=1.66,
                                  protValCryoEM=protValCryoEM)