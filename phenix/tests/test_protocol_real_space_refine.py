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
from chimera.protocols import ChimeraProtOperate
from pwem.protocols.protocol_import import (ProtImportPdb,
                                                    ProtImportVolumes)
from phenix.protocols.protocol_real_space_refine import (PhenixProtRunRSRefine,
                                                         mmCIF)
from phenix.protocols.protocol_molprobity import PhenixProtRunMolprobity
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

    def _importVolume3(self):
        args = {'filesPath': self.dsModBuild.getFile(
            'volumes/emd_4116.map'),
                'samplingRate': 0.637,
                'setOrigCoord': False
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_4116\nwith default '
                                  'origin\n')
        self.launchProtocol(protImportVol)
        volume3 = protImportVol.outputVolume
        return volume3

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

    def _importStructurePDBWoVol2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/3i3e_fitted.pdb')
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 3i3e_fitted')
        self.launchProtocol(protImportPDB)
        structure5_PDB = protImportPDB.outputPdb
        return structure5_PDB

    def _importStructurePDBWithVol2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/3i3e_fitted.pdb'),
                'inputVolume': self._importVolume3()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 3i3e_fitted with\n associated '
                                  'volume')
        self.launchProtocol(protImportPDB)
        structure6_PDB = protImportPDB.outputPdb
        return structure6_PDB

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

    def _importStructHemoCIF(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': '5ni1'
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import cif\n 5ni1.cif')
        self.launchProtocol(protImportPDB)
        structure_hemo_cif = protImportPDB.outputPdb
        return structure_hemo_cif

    def _importStructHemoCIF2(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': '5ni1',
                'inputVolume': self._importVolHemoOrg()
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import cif\n 5ni1.cif\nassociated volume '
                                  'emd_3488.map')
        self.launchProtocol(protImportPDB)
        structure_hemo_cif2 = protImportPDB.outputPdb
        self.assertTrue(structure_hemo_cif2.getFileName())
        return structure_hemo_cif2

class TestPhenixRSRefine(TestImportData):
    """ Test the protocol of Real Space Refine
    """
    def checkMPResults(self, ramOutliers, ramFavored, rotOutliers,
                       cbetaOutliers, clashScore, overallScore,
                       protMolProbity, places=2):
        # method to check MolProbity statistic results of the Final Results
        # Table
        try:
            self.assertAlmostEqual(protMolProbity.ramachandranOutliers.get(),
                                   ramOutliers, places)
            self.assertAlmostEqual(protMolProbity.ramachandranFavored.get(),
                                   ramFavored, places)
            self.assertAlmostEqual(protMolProbity.rotamerOutliers.get(),
                                   rotOutliers, places)
            self.assertAlmostEqual(protMolProbity.cbetaOutliers.get(),
                                   cbetaOutliers, places)
            self.assertAlmostEqual(protMolProbity.clashscore.get(),
                                   clashScore, places)
            self.assertAlmostEqual(protMolProbity.overallScore.get(),
                                   overallScore, places)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def checkRSRefineResults(self, ramOutliers, ramFavored, rotOutliers,
                             cbetaOutliers, clashScore, overallScore,
                             protRSRefine, places=3, delta=1):
        # method to check MolProbity statistic results of the Final Results
        # Table
        try:
            self.assertAlmostEqual(protRSRefine.ramachandranOutliers.get(),
                                   ramOutliers, places)
            self.assertAlmostEqual(protRSRefine.ramachandranFavored.get(),
                                   ramFavored, delta=3)
            self.assertAlmostEqual(protRSRefine.rotamerOutliers.get(),
                                   rotOutliers, delta=5)
            self.assertAlmostEqual(protRSRefine.cbetaOutliers.get(),
                                   cbetaOutliers, places)
            self.assertAlmostEqual(protRSRefine.clashscore.get(),
                                   clashScore, delta=delta)
            self.assertAlmostEqual(protRSRefine.overallScore.get(),
                                   overallScore, delta=delta)
        except Exception as e:
            # print error since test does not print it
            print(("Exception error:", str(e)))
            raise self.failureException(str(e))

    def testPhenixRSRefineFromPDB(self):
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

        # real_space_refine
        args = {
                'resolution': 3.5,
                'inputStructure': structure_refmac3,
                'numberOfThreads': 4
                }
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine without\nvolume, should NOT work')
        try:
            self.launchProtocol(protRSRefine)
        except Exception as e:
            self.assertTrue(True, "This test should return a error message as '"
                            " Input volume cannot be EMPTY.'\n")

            return
        self.assertTrue(False)

    def testPhenixRSRefineFromVolumeAndPDB1(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol, the input PDB was fitted
        to the volume and refined previously by coot and refmac withouth mask
        in another project; (MolProbity has been run to compare values before
        and after refinement); default refine strategy
        """
        print ("Run phenix real_space_refine from imported volume and pdb file "
               "previously fitted and refined by Coot and Refmac without mask "
               "(MolProbity has been run to compare values before and after "
               "refinement); default refine strategy")

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
                'inputStructure': structure_refmac3,
                'numberOfThreads': 4
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
        elif Plugin.getPhenixVersion() == PHENIXVERSION19:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.56,
                                      rotOutliers=1.70,
                                      cbetaOutliers=0,
                                      clashScore=2.09,
                                      overallScore=1.16,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.0,
                                      ramFavored=98.58,
                                      rotOutliers=2.84,
                                      cbetaOutliers=0,
                                      clashScore=1.79,
                                      overallScore=1.28,
                                      protRSRefine=protRSRefine)
        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.70,
                                      rotOutliers=3.98,
                                      cbetaOutliers=0,
                                      clashScore=4.47,
                                      overallScore=1.89,
                                      protRSRefine=protRSRefine)

    def testPhenixRSRefineFromVolumeAndPDB4(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks; default refine strategy; (MolProbity has been run to
        compare values before and after refinement).
        """
        print("Run phenix real_space_refine from imported volume and pdb file "
              "from data banks (vol origin 0.0, 0.0, 0.0); default refine "
              "strategy; (MolProbity has been run to compare values before "
              "and after refinement).")

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
                                   'volume and pdb\n')
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

    def testPhenixRSRefineFromVolumeAndCIF5(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks; default refine strategy; (MolProbity has been run to
        compare values before and after refinement).
        """
        print ("Run phenix real_space_refine from imported volume and cif file "
               "from data banks (vol origin 0.0, 0.0, 0.0); default refine "
               "strategy; (MolProbity has been run to compare values before "
               "and after refinement).")

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
                                   'volume and cif\n')
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
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.cif\ndefault refine strategy')
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
            # self.checkRSRefineResults(ramOutliers=0.00,
            #                           ramFavored=95.41,
            #                           rotOutliers=3.47,
            #                           cbetaOutliers=0,
            #                           clashScore=4.75,
            #                           overallScore=1.98,
            #                           protRSRefine=protRSRefine)

    def testPhenixRSRefineFromVolumeAndPDB6(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks; alternative refine strategy; (MolProbity has been run to
        compare values before and after refinement).
        """
        print ("Run phenix real_space_refine from imported volume and pdb "
               "file from data banks (vol origin 0.0, 0.0, 0.0); alternative "
               "refine strategy; (MolProbity has been run to compare values "
               "before and after refinement).")

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
                                   'volume and pdb\n')
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
                'localGridSearch': True,
                'morphing': True,
                'simulatedAnnealing': True,
                'adp': False,
                'numberOfThreads': 4
                }
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.pdb\nalternative refine strategy')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.35,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=2.10,
                                      overallScore=1.11,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.89,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=3.76,
                                      overallScore=1.19,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=98.06,
                                      rotOutliers=0.22,
                                      cbetaOutliers=0,
                                      clashScore=5.52,
                                      overallScore=1.30,
                                      protRSRefine=protRSRefine)
        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.29,
                                      rotOutliers=3.47,
                                      cbetaOutliers=0,
                                      clashScore=4.86,
                                      overallScore=1.50,
                                      protRSRefine=protRSRefine)

    # TODO: Talk to Roberto if we have to continue testing these values (rotOutliers)

    def testPhenixRSRefineFromVolumeAndCIF7(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks; alternative refine strategy; (MolProbity has been run
        to compare values before and after refinement).
        """
        print("Run phenix real_space_refine from imported volume and cif "
              "file from data banks (vol origin 0.0, 0.0, 0.0); "
              "alternative refine strategy; (MolProbity has been run to "
              "compare values before and after refinement).")

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
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'volume and cif\n')
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
                'localGridSearch': True,
                'morphing': True,
                'simulatedAnnealing': True,
                'adp': False
                }
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.cif\nalternative refine strategy')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.35,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=2.10,
                                      overallScore=1.11,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.52,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=3.64,
                                      overallScore=1.25,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.70,
                                      rotOutliers=0.22,
                                      cbetaOutliers=0,
                                      clashScore=5.63,
                                      overallScore=1.37,
                                      protRSRefine=protRSRefine)
        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.29,
                                      rotOutliers=3.47,
                                      cbetaOutliers=0,
                                      clashScore=4.86,
                                      overallScore=1.50,
                                      protRSRefine=protRSRefine)
            # self.checkRSRefineResults(ramOutliers=0.00,
            #                           ramFavored=96.47,
            #                           rotOutliers=0.00,
            #                           cbetaOutliers=0,
            #                           clashScore=4.20,
            #                           overallScore=1.43,
            #                           protRSRefine=protRSRefine)
    # TODO: Talk to Roberto if we have to continue testing these values (rotOutliers)

    def testPhenixRSRefineFromVolumeAndCIF8(self):
        """ This test checks that phenix real_space_refine protocol runs
        with a volume provided directly as inputVol and the input PDB from
        data banks; default refine strategy; (MolProbity has been run to
        compare values before and after refinement); outputFormat = mmCIF.
        """
        print ("Run phenix real_space_refine from imported volume and cif file "
               "from data banks (vol origin 0.0, 0.0, 0.0); default refine "
               "strategy; (MolProbity has been run to compare values before "
               "and after refinement); outputFormat = mmCIF.")

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
               }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                       'volume and cif\n')
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
                'outputFormat': mmCIF
                }
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n emd_3488.map and '
                                 '5ni1.cif\ndefault refine strategy')
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
            # self.checkRSRefineResults(ramOutliers=0.00,
            #                           ramFavored=95.41,
            #                           rotOutliers=3.47,
            #                           cbetaOutliers=0,
            #                           clashScore=4.75,
            #                           overallScore=1.98,
            #                           protRSRefine=protRSRefine)

    def testPhenixRSRefineFromVolumeAndCIF9(self):
        """ This test checks that phenix real_space_refine protocol runs
            with a volume associated to the input mmCIF from the
            data bank; alternative refine strategy; (MolProbity has been run
            to compare values before and after refinement); outputFormat =
            mmCIF.
        """
        print("Run phenix real_space_refine from imported cif "
              "file from data banks (vol origin 0.0, 0.0, 0.0) and "
              "associated volume; "
              "alternative refine strategy; (MolProbity has been run to "
              "compare values before and after refinement); outputFormat = "
              "mmCIF.")

        # import PDB
        structure_hemo_cif2 = self._importStructHemoCIF2()

        # chimera operate to repair cif file
        extraCommands = ""
        extraCommands += "scipionwrite #3 " \
                         "prefix repaired_CIF_ChimeraX_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': structure_hemo_cif2
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        protChimera.setObjLabel('chimera operate\n repairing CIF\n')
        self.launchProtocol(protChimera)
        result = eval("protChimera.repaired_CIF_ChimeraX_Atom_struct__3_%06d" % \
                      protChimera.getObjId())

        # MolProbity
        args = {
                'resolution': 3.2,
                'inputStructure': result,
                }
        protMolProbity = self.newProtocol(PhenixProtRunMolprobity, **args)
        protMolProbity.setObjLabel('MolProbity validation\n'
                                   'cif and associated map\n')
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
        args = {
                'resolution': 3.2,
                'inputStructure': structure_hemo_cif2,
                'outputFormat': mmCIF,
                'localGridSearch': True,
                'morphing': True,
                'simulatedAnnealing': True,
                'adp': False
                }
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args['doSecondary'] = False
        protRSRefine = self.newProtocol(PhenixProtRunRSRefine, **args)
        protRSRefine.setObjLabel('RSRefine hemo\n'
                                 '5ni1.cif and associated map\n'
                                 'emd_3488.map\nalternative '
                                 'refine strategy')
        self.launchProtocol(protRSRefine)

        # check real_space_refine results
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.35,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=2.10,
                                      overallScore=1.11,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION18:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.52,
                                      rotOutliers=0.00,
                                      cbetaOutliers=0,
                                      clashScore=3.64,
                                      overallScore=1.25,
                                      protRSRefine=protRSRefine)
        elif Plugin.getPhenixVersion() == PHENIXVERSION20:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=97.70,
                                      rotOutliers=0.22,
                                      cbetaOutliers=0,
                                      clashScore=5.63,
                                      overallScore=1.37,
                                      protRSRefine=protRSRefine)
        else:
            self.checkRSRefineResults(ramOutliers=0.00,
                                      ramFavored=96.29,
                                      rotOutliers=3.47,
                                      cbetaOutliers=0,
                                      clashScore=4.86,
                                      overallScore=1.50,
                                      protRSRefine=protRSRefine)
            # self.checkRSRefineResults(ramOutliers=0.00,
            #                           ramFavored=96.47,
            #                           rotOutliers=0.00,
            #                           cbetaOutliers=0,
            #                           clashScore=4.20,
            #                           overallScore=1.43,
            #                           protRSRefine=protRSRefine)

    # TODO: Talk to Roberto if we have to continue testing these values (rotOutliers)
