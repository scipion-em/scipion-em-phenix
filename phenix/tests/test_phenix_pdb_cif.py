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

# protocol to test the different phenix methods with .pdb and .cif files
from chimera.protocols import ChimeraProtOperate
from pwem.protocols.protocol_import import (ProtImportPdb,
                                                    ProtImportVolumes)
from pyworkflow import Config

from phenix.protocols.protocol_real_space_refine import (PhenixProtRunRSRefine,
                                                         mmCIF)
from phenix.protocols.protocol_molprobity import PhenixProtRunMolprobity
from phenix.protocols.protocol_validation_cryoem import PhenixProtRunValidationCryoEM
from pyworkflow.tests import *
from phenix import Plugin, PHENIXVERSION, PHENIXVERSION18, PHENIXVERSION19, PHENIXVERSION20
import json
from phenix.protocols.protocol_emringer import PhenixProtRunEMRinger
from phenix.protocols.protocol_superpose_pdbs import PhenixProtRunSuperposePDBs


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import map volumes and atomic structures(PDBx/mmCIF files)
    """
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

    def _importStructHemoCIF(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/5ni1.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import cif\n 5ni1.cif')
        self.launchProtocol(protImportPDB)
        structure_hemo_cif = protImportPDB.outputPdb
        return structure_hemo_cif

    def _importStructHemoFromDB(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': '5ni1'
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import cif from PDB\n 5ni1.cif')
        self.launchProtocol(protImportPDB)
        structure_hemo_cif_PDB = protImportPDB.outputPdb
        return structure_hemo_cif_PDB

class TestPhenixPdbCif(TestImportData):
    """ Test every phenix protocol
    """
    def checkEMResults(self, optThresh, rotRatio, maxZscore, modLength,
                     EMScore, protEMRinger, places=4):
        # method to check EMRinger statistic results of the Final Results Table
        textFileName = protEMRinger._getExtraPath(
            protEMRinger.EMRINGERTRANSFERFILENAME.replace('py', 'txt'))
        with open(textFileName, "r") as f:
            self.resultsDict = json.loads(str(f.read()))
        try:
            self.assertAlmostEqual(self.resultsDict[
                                'Optimal Threshold'], optThresh, delta=0.5)
            self.assertAlmostEqual(self.resultsDict[
                                'Rotamer-Ratio'], rotRatio, delta=0.5)
            self.assertAlmostEqual(self.resultsDict[
                                'Max Zscore'], maxZscore, delta=0.5)
            self.assertAlmostEqual(self.resultsDict[
                                'Model Length'], modLength, delta=0.5)
            self.assertAlmostEqual(self.resultsDict[
                                'EMRinger Score'], EMScore, delta=0.5)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def checkMPResults(self, ramOutliers, ramFavored, rotOutliers,
                       cbetaOutliers, clashScore, overallScore,
                       protMolProbity, places=2, delta=1):
        # method to check MolProbity statistic results of the Final Results
        # Table
        try:
            self.assertAlmostEqual(protMolProbity.ramachandranOutliers.get(),
                                   ramOutliers, delta=delta)
            self.assertAlmostEqual(protMolProbity.ramachandranFavored.get(),
                                   ramFavored, delta=delta)
            self.assertAlmostEqual(protMolProbity.rotamerOutliers.get(),
                                   rotOutliers, delta=3)
            self.assertAlmostEqual(protMolProbity.cbetaOutliers.get(),
                                   cbetaOutliers, delta=delta)
            self.assertAlmostEqual(protMolProbity.clashscore.get(),
                                   clashScore, delta=delta)
            self.assertAlmostEqual(protMolProbity.overallScore.get(),
                                   overallScore, delta=delta)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def checkRSRefineResults(self, ramOutliers, ramFavored, rotOutliers,
                             cbetaOutliers, clashScore, overallScore,
                             protRSRefine, places=2):
        # method to check Real Space Refine statistic results of the Final Results
        # Table
        try:
            self.assertAlmostEqual(protRSRefine.ramachandranOutliers.get(),
                                   ramOutliers, places)
            self.assertAlmostEqual(protRSRefine.ramachandranFavored.get(),
                                   ramFavored, delta=1)
            self.assertAlmostEqual(protRSRefine.rotamerOutliers.get(),
                                   rotOutliers, delta=3)
            self.assertAlmostEqual(protRSRefine.cbetaOutliers.get(),
                                   cbetaOutliers, places)
            self.assertAlmostEqual(protRSRefine.clashscore.get(),
                                   clashScore, delta=0.75)
            self.assertAlmostEqual(protRSRefine.overallScore.get(),
                                   overallScore, delta=0.75)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def checkValCryoEMResults(self, ramOutliers, ramFavored, rotOutliers,
                             cbetaOutliers, clashScore, overallScore,
                             protValCryoEM, places=2):
        # method to check Validation CryoEM statistic results of the Final Results
        # Table
        try:
            print("Checking ramOutliers")
            self.assertAlmostEqual(protValCryoEM.ramachandranOutliers.get(),
                                   ramOutliers, places)
            print("Checking ramFavored")
            self.assertAlmostEqual(protValCryoEM.ramachandranFavored.get(),
                                   ramFavored, delta=1)
            print("Checking rotOutliers")
            self.assertAlmostEqual(protValCryoEM.rotamerOutliers.get(),
                                   rotOutliers, delta=4)
            print("Checking cbetaOutliers")
            self.assertAlmostEqual(protValCryoEM.cbetaOutliers.get(),
                                   cbetaOutliers, places)
            print("Checking clashScore")
            self.assertAlmostEqual(protValCryoEM.clashscore.get(),
                                   clashScore, delta=2)
            print("Checking overallScore")
            self.assertAlmostEqual(protValCryoEM.overallScore.get(),
                                   overallScore, delta=0.75)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def checkSuperposeResults(self, startRMSD, finalRMSD, protSuperposePdbs, places=2):
        # method to check start and final RMSD values showed in Summary
        logFile = os.path.abspath(protSuperposePdbs._getLogsPath()) +\
                  "/run.stdout"
        protSuperposePdbs._parseLogFile(logFile)
        try:
            self.assertAlmostEqual(protSuperposePdbs.startRMSD.get(),
                                   startRMSD, places)
            self.assertAlmostEqual(protSuperposePdbs.finalRMSD.get(),
                                   finalRMSD, places)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def testEMRingerValidationFromVolumeAndPDB(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol and the input PDB fitted to
        the volume
         """
        print("Run EMRinger validation from imported volume and pdb file")

        # Import Volume
        volume_hemo_orig = self._importVolHemoOrg()

        # import PDB
        structure_hemo_pdb = self._importStructHemoPDB()

        # EMRinger
        args = {'inputVolume': volume_hemo_orig,
                'inputStructure': structure_hemo_pdb,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n hemoglobin 5ni1.pdb\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        # Values obtained fron phenix GUI v.1.16
        self.checkEMResults(optThresh=0.07,
                            rotRatio=0.77,
                            maxZscore=7.65,
                            modLength=339,
                            EMScore=4.16,
                            protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolumeAndCIF(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol and input PDB
         """
        print("Run EMRinger validation from imported volume and cif file")

        # Import Volume
        volume_hemo_orig = self._importVolHemoOrg()

        # import CIF
        structure_hemo_cif = self._importStructHemoCIF()

        # chimera operate to repair cif file
        extraCommands = ""
        extraCommands += "scipionwrite #2 " \
                         "prefix repaired_CIF_ChimeraX_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure_hemo_cif
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        protChimera.setObjLabel('chimera operate\n repairing CIF\n')
        self.launchProtocol(protChimera)
        result = eval("protChimera.repaired_CIF_ChimeraX_Atom_struct__2_%06d" % \
                 protChimera.getObjId())

        # EMRinger
        args = {'inputVolume': volume_hemo_orig,
                'inputStructure': result,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n hemoglobin 5ni1.cif\n')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        # Values obtained fron phenix GUI v.1.16
        self.checkEMResults(optThresh=0.07,
                            rotRatio=0.77,
                            maxZscore=7.65,
                            modLength=339,
                            EMScore=4.16,
                            protEMRinger=protEMRinger)

    def testEMRingerValidationFromVolumeAndCIFFromPDB(self):
        """ This test checks that EMRinger validation protocol runs with a
        volume provided directly as inputVol and input PDB
         """
        print("Run EMRinger validation from imported volume and cif file from PDB")

        # Import Volume
        volume_hemo_orig = self._importVolHemoOrg()

        # import CIF
        structure_hemo_cif_PDB = self._importStructHemoFromDB()

        # EMRinger
        args = {'inputVolume': volume_hemo_orig,
                'inputStructure': structure_hemo_cif_PDB,
                'doTest': True
                }
        protEMRinger = self.newProtocol(PhenixProtRunEMRinger, **args)
        protEMRinger.setObjLabel('EMRinger validation\n hemoglobin 5ni1.cif\nfrom PDB')
        self.launchProtocol(protEMRinger)

        # check EMRinger results
        self.checkEMResults(optThresh=0.07,
                            rotRatio=0.77,
                            maxZscore=7.65,
                            modLength=339,
                            EMScore=4.16,
                            protEMRinger=protEMRinger)

    def testValCryoEMFFromVolumeAndPDB(self):
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
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.64,
                                      rotOutliers=4.77,
                                      cbetaOutliers=0,
                                      clashScore=6.07,
                                      overallScore=2.06,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.35,
                                      rotOutliers=2.39,
                                      cbetaOutliers=0,
                                      clashScore=3.87,
                                      overallScore=1.59,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.59,
                                      rotOutliers=1.30,
                                      cbetaOutliers=0,
                                      clashScore=4.31,
                                      overallScore=1.30,
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
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.59,
                                      rotOutliers=1.30,
                                      cbetaOutliers=0,
                                      clashScore=4.31,
                                      overallScore=1.30,
                                      protRSRefine=protRSRefine)
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
        if Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=97.35,
                                       rotOutliers=2.39,
                                       cbetaOutliers=0,
                                       clashScore=3.87,
                                       overallScore=1.59,
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
                                      ramFavored=96.11,
                                      rotOutliers=3.04,
                                      cbetaOutliers=0,
                                      clashScore=5.30,
                                      overallScore=1.92,
                                      protValCryoEM=protValCryoEM)

    def testValCryoEMFFromVolumeAndCIF(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol (with half1 and half2)
        and the input CIF from data banks; default refine strategy;
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
        volume_hemo_org = self._importVolHemoOrg()

        # import PDB
        structure_hemo_cif = self._importStructHemoCIF()

        # chimera operate to repair cif file
        extraCommands = ""
        extraCommands += "scipionwrite #2 " \
                         "prefix repaired_CIF_ChimeraX_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure_hemo_cif
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        protChimera.setObjLabel('chimera operate\n repairing CIF\n')
        self.launchProtocol(protChimera)
        result = eval("protChimera.repaired_CIF_ChimeraX_Atom_struct__2_%06d" % \
                 protChimera.getObjId())
        # MolProbity
        args = {
                'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': result,
                'numberOfThreads': 4
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'volume\nand cif\nhemoglobin')
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
                'inputStructure': structure_hemo_cif,
                'numberOfThreads': 4
                }
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map\nand '
                                 '5ni1.cif\n')
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
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.64,
                                      rotOutliers=4.77,
                                      cbetaOutliers=0,
                                      clashScore=6.07,
                                      overallScore=2.06,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.35,
                                      rotOutliers=2.39,
                                      cbetaOutliers=0,
                                      clashScore=3.87,
                                      overallScore=1.59,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.59,
                                      rotOutliers=1.30,
                                      cbetaOutliers=0,
                                      clashScore=4.31,
                                      overallScore=1.30,
                                      protRSRefine=protRSRefine)
        else:
            # values obtained from phenix GUI v. 1.16
            # (minimization_global + adp)
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
                                   'volume\nand cif\nhemoglobin')
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
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.59,
                                      rotOutliers=1.30,
                                      cbetaOutliers=0,
                                      clashScore=4.31,
                                      overallScore=1.30,
                                      protRSRefine=protRSRefine)
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

        # check validation cryoem results
        if Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                    ramFavored=97.35,
                                    rotOutliers=2.39,
                                    cbetaOutliers=0,
                                    clashScore=3.87,
                                    overallScore=1.59,
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
                                      ramFavored=96.11,
                                      rotOutliers=3.04,
                                      cbetaOutliers=0,
                                      clashScore=5.30,
                                      overallScore=1.92,
                                      protValCryoEM=protValCryoEM)

    def testValCryoEMFFromVolumeAndCIFFromPDB(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol (with half1 and half2)
        and the input CIF from data banks; default refine strategy;
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
        volume_hemo_org = self._importVolHemoOrg()

        # import cif
        structure_hemo_cif_PDB = self._importStructHemoFromDB()

        # chimera operate to repair cif file
        extraCommands = ""
        extraCommands += "scipionwrite #2 " \
                         "prefix repaired_CIF_ChimeraX_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure_hemo_cif_PDB
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        protChimera.setObjLabel('chimera operate\n repairing CIF\n')
        self.launchProtocol(protChimera)
        result = eval("protChimera.repaired_CIF_ChimeraX_Atom_struct__2_%06d" % \
                 protChimera.getObjId())

        # MolProbity
        args = {
                'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': result,
                'numberOfThreads': 4
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'volume\nand cif from PDB\nhemoglobin')
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
                'inputStructure': structure_hemo_cif_PDB,
                'numberOfThreads': 4
                }
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map\nand '
                                 '5ni1.cif from PDB\n')
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

        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.64,
                                      rotOutliers=4.77,
                                      cbetaOutliers=0,
                                      clashScore=6.07,
                                      overallScore=2.06,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.35,
                                      rotOutliers=2.39,
                                      cbetaOutliers=0,
                                      clashScore=3.87,
                                      overallScore=1.59,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.59,
                                      rotOutliers=1.30,
                                      cbetaOutliers=0,
                                      clashScore=4.31,
                                      overallScore=1.30,
                                      protRSRefine=protRSRefine)
        else:
            # values obtained from phenix GUI v. 1.16
            # (minimization_global + adp)
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
                                   'volume\nand cif\nhemoglobin')
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
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.59,
                                      rotOutliers=1.30,
                                      cbetaOutliers=0,
                                      clashScore=4.31,
                                      overallScore=1.30,
                                      protRSRefine=protRSRefine)
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

        # check validation cryoem results
        if Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkValCryoEMResults(ramOutliers=0.00,
                                       ramFavored=97.35,
                                       rotOutliers=2.39,
                                       cbetaOutliers=0,
                                       clashScore=3.87,
                                       overallScore=1.59,
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
                                      ramFavored=96.11,
                                      rotOutliers=3.04,
                                      cbetaOutliers=0,
                                      clashScore=5.30,
                                      overallScore=1.92,
                                      protValCryoEM=protValCryoEM)

    def testSuperposePdbsFromPDBAndCIF(self):
        """ This test checks that phenix superpose_pdbs protocol runs with
        two atomic structures (pdb and cif)"""
        print("Run phenix superpose_pdbs protocol from an imported pdb file " \
              "and an imported cif file without volumes associated")

        # Import Volume
        volume_hemo_org = self._importVolHemoOrg()

        # import PDBs
        structure_PDB = self._importStructHemoPDB()
        self.assertTrue(structure_PDB.getFileName())
        self.assertFalse(structure_PDB.getVolume())

        structure_CIF = self._importStructHemoCIF()
        self.assertTrue(structure_CIF.getFileName())
        self.assertFalse(structure_CIF.getVolume())

        args = {
            'inputStructureFixed': structure_PDB,
            'inputStructureMoving': structure_CIF
        }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs, **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                      'pdb fixed\n')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))

        # check RMSD results
        self.checkSuperposeResults(startRMSD=0.000,
                                   finalRMSD=0.000,
                                   protSuperposePdbs=protSuperposePdbs)

        # validation_cryoEM
        args = {'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': protSuperposePdbs.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter Superpose\nvolume_hemo_org and\n'
                                  'protSuperposePdbs.outputPdb\n')
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
                                   ramFavored=95.23,
                                   rotOutliers=0.43,
                                   cbetaOutliers=0,
                                   clashScore=3.53,
                                   overallScore=1.48,
                                   protValCryoEM=protValCryoEM)
        args = {
            'inputStructureFixed': structure_CIF,
            'inputStructureMoving': structure_PDB
        }

        protSuperposePdbs = self.newProtocol(PhenixProtRunSuperposePDBs, **args)
        protSuperposePdbs.setObjLabel('SuperposePDBs\n'
                                      'cif fixed\n')
        self.launchProtocol(protSuperposePdbs)
        self.assertTrue(os.path.exists(
            protSuperposePdbs.outputPdb.getFileName()))

        # check RMSD results
        self.checkSuperposeResults(startRMSD=0.000,
                                   finalRMSD=0.000,
                                   protSuperposePdbs=protSuperposePdbs)

        # validation_cryoEM
        args = {'inputVolume': volume_hemo_org,
                'resolution': 3.2,
                'inputStructure': protSuperposePdbs.outputPdb
                }
        protValCryoEM = self.newProtocol(PhenixProtRunValidationCryoEM, **args)
        protValCryoEM.setObjLabel('ValCryoEM\nafter Superpose\nvolume_hemo_org and\n'
                                  'protSuperposePdbs.outputPdb\n')
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
                                   ramFavored=95.23,
                                   rotOutliers=0.43,
                                   cbetaOutliers=0,
                                   clashScore=3.53,
                                   overallScore=1.48,
                                   protValCryoEM=protValCryoEM)

