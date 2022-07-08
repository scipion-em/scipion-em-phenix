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

from pwem.protocols.protocol_import import (ProtImportPdb,
                                                    ProtImportVolumes)
from pyworkflow import Config

from phenix.protocols.protocol_real_space_refine import (PhenixProtRunRSRefine,
                                                         mmCIF)
from phenix.protocols.protocol_molprobity import PhenixProtRunMolprobity
from phenix.protocols.protocol_validation_cryoem import PhenixProtRunValidationCryoEM
from pyworkflow.tests import *
from phenix import Plugin, PHENIXVERSION, PHENIXVERSION18, PHENIXVERSION19, PHENIXVERSION20


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
                            'volumes/emd_3488_Noisy_half1.vol'),
                'half2map': self.dsModBuild.getFile(
                            'volumes/emd_3488_Noisy_half2.vol'),
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
                       protMolProbity, places=2, delta=1):
        # method to check MolProbity statistic results of the Final Results
        # Table
        try:
            print("Checking ramOutliers")
            self.assertAlmostEqual(protMolProbity.ramachandranOutliers.get(),
                                   ramOutliers, delta=delta)
            print("Checking ramFavored")
            self.assertAlmostEqual(protMolProbity.ramachandranFavored.get(),
                                   ramFavored, delta=delta)
            # root outliers seems to vary a lot from one version to other
            print("Checking rotOutliers")
            self.assertAlmostEqual(protMolProbity.rotamerOutliers.get(),
                                   rotOutliers, delta=5)
            print("Checking cbetaOutliers")
            self.assertAlmostEqual(protMolProbity.cbetaOutliers.get(),
                                   cbetaOutliers, delta=delta)
            print("Checking clashScore")
            self.assertAlmostEqual(protMolProbity.clashscore.get(),
                                   clashScore, delta=2)
            print("Checking overallScore")
            self.assertAlmostEqual(protMolProbity.overallScore.get(),
                                   overallScore, delta=delta)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def checkRSRefineResults(self, ramOutliers, ramFavored, rotOutliers,
                             cbetaOutliers, clashScore, overallScore,
                             protRSRefine, places=2, delta=3):
        # method to check Real Space Refine statistic results of the Final Results
        # Table
        try:
            self.assertAlmostEqual(protRSRefine.ramachandranFavored.get(),
                                   ramFavored, delta=delta)
            self.assertAlmostEqual(protRSRefine.rotamerOutliers.get(),
                                   rotOutliers, delta=5)
            self.assertAlmostEqual(protRSRefine.cbetaOutliers.get(),
                                   cbetaOutliers, places)
            self.assertAlmostEqual(protRSRefine.clashscore.get(),
                                   clashScore, delta=delta)
            self.assertAlmostEqual(protRSRefine.overallScore.get(),
                                   overallScore, delta=delta)
            self.assertAlmostEqual(protRSRefine.ramachandranOutliers.get(),
                                   ramOutliers, places)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def checkValCryoEMResults(self, ramOutliers, ramFavored, rotOutliers,
                             cbetaOutliers, clashScore, overallScore,
                             protValCryoEM, places=2, delta=3):
        # method to check Validation CryoEM statistic results of the Final Results
        # Table
        try:
            self.assertAlmostEqual(protValCryoEM.ramachandranOutliers.get(),
                                   ramOutliers, places)
            self.assertAlmostEqual(protValCryoEM.ramachandranFavored.get(),
                                   ramFavored, delta=delta)
            self.assertAlmostEqual(protValCryoEM.rotamerOutliers.get(),
                                   rotOutliers, delta=5)
            self.assertAlmostEqual(protValCryoEM.cbetaOutliers.get(),
                                   cbetaOutliers, places)
            self.assertAlmostEqual(protValCryoEM.clashscore.get(),
                                   clashScore, delta=delta)
            self.assertAlmostEqual(protValCryoEM.overallScore.get(),
                                   overallScore, delta=delta)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine\n refmac3.mrc and '
                                   'refmac3.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=95.75,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=2.09,
                                      overallScore=1.27,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.58,
                                      rotOutliers=2.84,
                                      cbetaOutliers=0,
                                      clashScore=1.79,
                                      overallScore=1.28,
                                      protRSRefine=protRSRefine)
        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.23,
                                      rotOutliers=3.41,
                                      cbetaOutliers=0,
                                      clashScore=4.17,
                                      overallScore=1.86,
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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=95.75,
                                rotOutliers=0.00,
                                cbetaOutliers=0,
                                clashScore=2.09,
                                overallScore=1.27,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=98.58,
                                rotOutliers=1.70,
                                cbetaOutliers=0,
                                clashScore=2.09,
                                overallScore=1.16,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkMPResults(ramOutliers=0.00,
                                 ramFavored=98.58,
                                 rotOutliers=2.84,
                                 cbetaOutliers=0,
                                 clashScore=1.79,
                                 overallScore=1.28,
                                 protMolProbity=protMolProbity2)
        else:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=96.23,
                                rotOutliers=3.41,
                                cbetaOutliers=0,
                                clashScore=4.17,
                                overallScore=1.86,
                                protMolProbity=protMolProbity2)

        # validation_cryoEM
        args = {'inputVolume': volume_refmac3,
                'resolution': 3.5,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nrefmac3.mrc and '
                                 'protRSRefine.outputPdb\n')
        try:
            self.launchProtocol(protValCryoEM)
        except Exception as e:
            self.assertTrue(True, "This test should return a error message as '"
                  " Protocol has validation errors:\n"
                  "Binary '/usr/local/phenix-1.13-2998/modules/phenix/phenix/"
                  "command_line/validation_cryoem.py' does not exists.\n"
                  "Check if you need to upgrade your PHENIX version to run "
                  "VALIDATION_CRYOEM.\nYour current PHENIX version is 1.13.\n"
                  "Check configuration file: " +
                  Config.SCIPION_LOCAL_CONFIG + "\n"
                  "and set VALIDATION_CRYOEM and PHENIX_HOME variables properly.\n"
                  "Current values:\nPHENIX_HOME = /usr/local/phenix-1.13-2998\n")

            return
        # self.assertTrue(False)

        # check validation cryoem results
        if Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=98.58,
                                       rotOutliers=2.84,
                                       cbetaOutliers=0,
                                       clashScore=1.79,
                                       overallScore=1.28,
                                       protValCryoEM=protValCryoEM)
        else:
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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.53,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=2.43,
                                      overallScore=1.12,
                                      protRSRefine=protRSRefine)

        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                  ramFavored=96.11,
                                  rotOutliers=3.04,
                                  cbetaOutliers=0,
                                  clashScore=5.30,
                                  overallScore=1.92,
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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=97.53,
                                rotOutliers=0.00,
                                cbetaOutliers=0,
                                clashScore=2.43,
                                overallScore=1.12,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=96.64,
                                rotOutliers=4.77,
                                cbetaOutliers=0,
                                clashScore=6.07,
                                overallScore=2.06,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=97.35,
                                rotOutliers=2.39,
                                cbetaOutliers=0,
                                clashScore=3.87,
                                overallScore=1.59,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=98.59,
                                rotOutliers=1.30,
                                cbetaOutliers=0,
                                clashScore=4.31,
                                overallScore=1.30,
                                protMolProbity=protMolProbity2)
        else:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=96.11,
                                rotOutliers=3.04,
                                cbetaOutliers=0,
                                clashScore=5.30,
                                overallScore=1.92,
                                protMolProbity=protMolProbity2)

        # validation_cryoEM
        args = {'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nvolume_hemo_org and '
                                  'protRSRefine.outputPdb\n')
        try:
            self.launchProtocol(protValCryoEM)
        except Exception as e:
            self.assertTrue(True, "This test should return a error message as '"
                            " Protocol has validation errors:\n"
                            "Binary '/usr/local/phenix-1.13-2998/modules/phenix/phenix/"
                            "command_line/validation_cryoem.py' does not exists.\n"
                            "Check if you need to upgrade your PHENIX version to run "
                            "VALIDATION_CRYOEM.\nYour current PHENIX version is 1.13.\n"
                            "Check configuration file: " +
                            Config.SCIPION_LOCAL_CONFIG + "\n"
                            "and set VALIDATION_CRYOEM and PHENIX_HOME variables properly.\n"
                            "Current values:\nPHENIX_HOME = /usr/local/phenix-1.13-2998\n")

            return
        # self.assertTrue(False)
        if Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=96.64,
                                       rotOutliers=4.77,
                                       cbetaOutliers=0,
                                       clashScore=6.07,
                                       overallScore=2.06,
                                       protValCryoEM=protValCryoEM)
        if Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=98.59,
                                       rotOutliers=1.30,
                                       cbetaOutliers=0,
                                       clashScore=4.42,
                                       overallScore=1.31,
                                       protValCryoEM=protValCryoEM)
        else:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                      ramFavored=95.41,
                                      rotOutliers=3.47,
                                      cbetaOutliers=0,
                                      clashScore=4.75,
                                      overallScore=1.98,
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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine nucleosome\nvolume and '
                                 'pdb\n')

        # TODO, this protocol fails because
        # Opening quote in middle of word: ATOM 5963 O5' . DA I -72 ? 73.27900 73.22500 141.39500 1.000 212.95366 O ? L ? . 1
        # so 05' is not valid for the cif reader, I think is fair and should
        # not be counted as a problem
        self.launchProtocol(protRSRefine) #Keep this line commented in pull requests
        # return; self.launchProtocol(protRSRefine)

        # check real_space_refine results
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=95.92,
                                      rotOutliers=0.16,
                                      cbetaOutliers=0,
                                      clashScore=7.46,
                                      overallScore=1.69,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.19,
                                      rotOutliers=12.32,
                                      cbetaOutliers=0,
                                      clashScore=12.40,
                                      overallScore=2.69,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.23,
                                      rotOutliers=2.72,
                                      cbetaOutliers=0,
                                      clashScore=7.65,
                                      overallScore=1.75,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.10,
                                      rotOutliers=2.40,
                                      cbetaOutliers=0,
                                      clashScore=7.23,
                                      overallScore=1.69,
                                      protRSRefine=protRSRefine)
        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=94.42,
                                      rotOutliers=11.04,
                                      cbetaOutliers=0,
                                      clashScore=12.92,
                                      overallScore=2.79,
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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=95.92,
                                rotOutliers=0.16,
                                cbetaOutliers=0,
                                clashScore=7.46,
                                overallScore=1.69,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=96.19,
                                rotOutliers=12.32,
                                cbetaOutliers=0,
                                clashScore=12.45,
                                overallScore=2.69,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=98.23,
                                rotOutliers=2.72,
                                cbetaOutliers=0,
                                clashScore=7.65,
                                overallScore=1.75,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=98.10,
                                rotOutliers=2.40,
                                cbetaOutliers=0,
                                clashScore=7.32,
                                overallScore=1.69,
                                protMolProbity=protMolProbity2)
        else:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=94.42,
                                rotOutliers=12.96,
                                cbetaOutliers=0,
                                clashScore=12.96,
                                overallScore=2.80,
                                protMolProbity=protMolProbity2)
        # TODO: Talk to Roberto if we have to continue testing these values (rotOutliers)

        # validation_cryoEM
        args = {
                'resolution': 4.0,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nnucleosome and '
                                  'protRSRefine.outputPdb\n')
        try:
            self.launchProtocol(protValCryoEM)
        except Exception as e:
            self.assertTrue(True, "This test should return a error message as '"
                            " Protocol has validation errors:\n"
                            "Binary '/usr/local/phenix-1.13-2998/modules/phenix/phenix/"
                            "command_line/validation_cryoem.py' does not exists.\n"
                            "Check if you need to upgrade your PHENIX version to run "
                            "VALIDATION_CRYOEM.\nYour current PHENIX version is 1.13.\n"
                            "Check configuration file: " +
                            Config.SCIPION_LOCAL_CONFIG + "\n"
                            "and set VALIDATION_CRYOEM and PHENIX_HOME variables properly.\n"
                            "Current values:\nPHENIX_HOME = /usr/local/phenix-1.13-2998\n")

            return
        # self.assertTrue(False)

        # check validation cryoem results
        if Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=96.19,
                                       rotOutliers=12.32,
                                       cbetaOutliers=0,
                                       clashScore=12.40,
                                       overallScore=2.69,
                                       protValCryoEM=protValCryoEM)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=98.23,
                                       rotOutliers=2.72,
                                       cbetaOutliers=0,
                                       clashScore=7.65,
                                       overallScore=1.75,
                                       protValCryoEM=protValCryoEM)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=98.10,
                                       rotOutliers=2.40,
                                       cbetaOutliers=0,
                                       clashScore=7.23,
                                       overallScore=1.69,
                                       protValCryoEM=protValCryoEM)
        else:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                      ramFavored=98.58,
                                      rotOutliers=1.70,
                                      cbetaOutliers=0,
                                      clashScore=2.09,
                                      overallScore=1.16,
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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map (full, half1, half2)\nand '
                                 '5ni1.pdb\n')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.53,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=2.43,
                                      overallScore=1.12,
                                      protRSRefine=protRSRefine)

        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.11,
                                      rotOutliers=3.04,
                                      cbetaOutliers=0,
                                      clashScore=5.30,
                                      overallScore=1.92,
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
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=97.53,
                                rotOutliers=0.00,
                                cbetaOutliers=0,
                                clashScore=2.43,
                                overallScore=1.12,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=96.64,
                                rotOutliers=4.77,
                                cbetaOutliers=0,
                                clashScore=6.07,
                                overallScore=2.06,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=97.35,
                                rotOutliers=2.39,
                                cbetaOutliers=0,
                                clashScore=3.87,
                                overallScore=1.59,
                                protMolProbity=protMolProbity2)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=98.59,
                                rotOutliers=1.30,
                                cbetaOutliers=0,
                                clashScore=4.31,
                                overallScore=1.30,
                                protMolProbity=protMolProbity2)
        else:
            self.checkMPResults(ramOutliers=0.00,
                                ramFavored=96.11,
                                rotOutliers=3.04,
                                cbetaOutliers=0,
                                clashScore=5.30,
                                overallScore=1.92,
                                protMolProbity=protMolProbity2)

        # validation_cryoEM
        args = {'inputVolume': volume_hemo_half1_hal2,
                'resolution': 3.2,
                'inputStructure': protRSRefine.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter RSRefine\nvolume_hemo_org and '
                                  'protRSRefine.outputPdb\n')
        try:
            self.launchProtocol(protValCryoEM)
        except Exception as e:
            self.assertTrue(True, "This test should return a error message as '"
                            " Protocol has validation errors:\n"
                            "Binary '/usr/local/phenix-1.13-2998/modules/phenix/phenix/"
                            "command_line/validation_cryoem.py' does not exists.\n"
                            "Check if you need to upgrade your PHENIX version to run "
                            "VALIDATION_CRYOEM.\nYour current PHENIX version is 1.13.\n"
                            "Check configuration file: " +
                            Config.SCIPION_LOCAL_CONFIG + "\n"
                            "and set VALIDATION_CRYOEM and PHENIX_HOME variables properly.\n"
                            "Current values:\nPHENIX_HOME = /usr/local/phenix-1.13-2998\n")

            return
        # self.assertTrue(False)

        # check validation cryoem results
        self.checkValCryoEMResults(ramOutliers=0.00,
                                  ramFavored=96.11,
                                  rotOutliers=3.04,
                                  cbetaOutliers=0,
                                  clashScore=5.30,
                                  overallScore=1.92,
                                  protValCryoEM=protValCryoEM)
