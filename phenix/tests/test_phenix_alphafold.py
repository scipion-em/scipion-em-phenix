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
from phenix.protocols import PhenixProtProcessPredictedAlphaFold2Model,  \
    PhenixProtDockPredictedAlphaFold2Model, PhenixProtRebuildDockPredictedAlphaFold2Model, \
    PhenixProtDockAndRebuildAlphaFold2Model
from pwem.protocols.protocol_import import (ProtImportPdb,
                                            ProtImportVolumes)
from pyworkflow.tests import *
from chimera.protocols import ChimeraProtOperate
from xmipp3.protocols import XmippProtExtractUnit
import pwem.protocols as emprot
import requests
import xml.etree.ElementTree as ET


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import maps and structure predictions from Alphafold2
    """
    def _importVolume1(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/emd_3488.map'),
                'samplingRate': 1.05,
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_3488.map\n')
        self.launchProtocol(protImportVol)
        volume1 = protImportVol.outputVolume
        return volume1

    def _importVolume2(self):
        args = {'importFrom': ProtImportVolumes.IMPORT_FROM_EMDB,
                'emdbId': 24230
                }
        protImportVol = self.newProtocol(emprot.ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume\n emd_24230.map\n from EMDB\n')
        self.launchProtocol(protImportVol)
        volume2 = protImportVol.outputVolume
        # volume.setOrigin(None)
        # The volume has no origin
        t = volume2.getOrigin(force=True)
        x, y, z = t.getShifts()
        # self.assertTrue(os.path.exists(volume2._getExtraPath('emd_24230.map')))
        return volume2

    def fetch_file(self, url, retry=3, json=False, outFilename=None):
        """ fetch file from url, retry 'retries' number of times
        :param str url: full url to file
        :param int reties: number of attemps
        :param bool json: parse json response
        :param outFilename
        """
        import requests
        for r in range(retry):
            try:
                n = os.path.basename(url)
                response = requests.get(url)  # Downloading the file and saving it at app/test with the file name n
            except Exception as e:
                if r < 2:
                    print(f'Failed. Attempt # {r + 1}')
                else:
                    print('Error encoutered at third attempt downloading {n}')
                    print(e)
            else:
                print(f"Success: {n}")
                break

        if json:
            return response.json()

        if outFilename is not None:
            # check if we have a valid answer
            if response.text.find('<?xml version=') != -1:
                # error case, parse return xml string
                root = ET.fromstring(response.text)
                errorMessage = root.find("Message")
                return False, errorMessage.text
            with open(outFilename, 'wb') as f:
                f.write(response.content)
                return True, ''
        else:  # return content
            return True, response.content

    def _get_alphafold_database_settings(self):
        """ get alphafold database settings from
             https://www.rbvi.ucsf.edu/chimerax/data/status/alphafold_database.json"""
        url = 'https://www.rbvi.ucsf.edu/chimerax/data/status/alphafold_database.json'
        self.settings = {}
        try:
            self.settings = self.fetch_file(url, json=True, outFilename=None)
        except Exception:
            print("Could not reach update site")

        if not self.settings:
            raise Exception('No alphafold database settings found')

    def _getAlphaFoldModelFromEBI(self, uniProtID):
        '''Copied from Roberto's protocol chimerax-protocol alphafold:
        Fetch structures from EBI AlphaFold database using UniProt sequence ID.
                   Example for UniProt P29474.
                   https://alphafold.ebi.ac.uk/files/AF-P29474-F1-model_v1.cif'''

        # get alphafold EBI database url
        self._get_alphafold_database_settings()
        data = {'uniprot_id': uniProtID, 'version': self.settings['database_version']}
        # get model
        extraDir = os.getcwd()
        model_url = self.settings['database_url'].format(**data)
        outFileName = os.path.join(extraDir, uniProtID + ".cif")
        status, msg = self.fetch_file(model_url, retry=3, json=False, outFilename=outFileName)
        if not status:
            error_message = f"ERROR: Can not retieve {uniProtID} from EBI Alphafold database. Is {uniProtID} a valid UNIPROT ID?"
            error_message += f"SYSTEM report error {msg}"
            raise Exception(error_message)
        else:
            return outFileName

    def _importStructure1(self):
        """Import atom structure prdiction by alphafold from EBI using uniprotid"""
        uniProtID = 'Q9BXJ8'
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self._getAlphaFoldModelFromEBI(uniProtID)
                }
        protImportAlphaFoldPDB = self.newProtocol(ProtImportPdb, **args)
        protImportAlphaFoldPDB.setObjLabel('import from EBI with \n uniProtID %s' % uniProtID)
        self.launchProtocol(protImportAlphaFoldPDB)
        return protImportAlphaFoldPDB

    def _importStructure2(self):
        """Import atom structure prdiction by alphafold from EBI using uniprotid"""
        uniProtID = 'P69905'
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self._getAlphaFoldModelFromEBI(uniProtID)
                }
        protImportAlphaFoldPDB = self.newProtocol(ProtImportPdb, **args)
        protImportAlphaFoldPDB.setObjLabel('import from EBI with \n uniProtID %s' % uniProtID)
        self.launchProtocol(protImportAlphaFoldPDB)
        return protImportAlphaFoldPDB

    def _importStructure3(self):
        """Import atom structure prdiction by alphafold from EBI using uniprotid"""
        uniProtID = 'A0A087WSY6'
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self._getAlphaFoldModelFromEBI(uniProtID)
                }
        protImportAlphaFoldPDB = self.newProtocol(ProtImportPdb, **args)
        protImportAlphaFoldPDB.setObjLabel('import from EBI with \n uniProtID %s' % uniProtID)
        self.launchProtocol(protImportAlphaFoldPDB)
        return protImportAlphaFoldPDB


class TestAProtProcessDockBuildPredictedAlphaFold(TestImportData):

    def testAProcessPrediction1(self):
        """ Test the protocol process alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure2()
        structure2 = protImportPDB.outputPdb
        self.assertTrue(structure2.getFileName())
        self.assertFalse(structure2.getVolume())

        args = {
                'inputPredictedModel': structure2,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nP69905\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

    def testBProcessDockPrediction1(self):
        """ Test the protocol process and dock alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure2()
        structure2 = protImportPDB.outputPdb
        self.assertTrue(structure2.getFileName())
        self.assertFalse(structure2.getVolume())

        args = {
                'inputPredictedModel': structure2,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nP69905\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

        print("Run phenix dock_predicted_model protocol from the previous"
              "predicted and processed atomic structure ")

        # import map
        map = self._importVolume1()
        self.assertTrue(map.getFileName())

        args = {
                'inputPredictedModel': structure2,
                'inputProcessedPredictedModel': protProcessPrediction.outputPdb,
                'inputVolume': map,
                'resolution': 3.2
                }
        protProcessDockPrediction = self.newProtocol(
            PhenixProtDockPredictedAlphaFold2Model, **args)
        protProcessDockPrediction.setObjLabel('Dock\nprediction processed\nP69905\n')
        self.launchProtocol(protProcessDockPrediction)
        self.assertTrue(os.path.exists(
            protProcessDockPrediction.outputPdb.getFileName()))

    def testCProcessDockBuildPrediction1(self):
        """ Test the protocol process, dock and rebuild alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure2()
        structure2 = protImportPDB.outputPdb
        self.assertTrue(structure2.getFileName())
        self.assertFalse(structure2.getVolume())

        args = {
                'inputPredictedModel': structure2,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nP69905\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

        print("Run phenix dock_predicted_model protocol from the previous"
              "predicted and processed atomic structure ")

        # import map
        map = self._importVolume1()
        self.assertTrue(map.getFileName())

        args = {
                'inputPredictedModel': structure2,
                'inputProcessedPredictedModel': protProcessPrediction.outputPdb,
                'inputVolume': map,
                'resolution': 3.2
                }
        protProcessDockPrediction = self.newProtocol(
            PhenixProtDockPredictedAlphaFold2Model, **args)
        protProcessDockPrediction.setObjLabel('Dock\nprediction processed\nP69905\n')
        self.launchProtocol(protProcessDockPrediction)
        self.assertTrue(os.path.exists(
            protProcessDockPrediction.outputPdb.getFileName()))

        print("Run phenix rebuild_predicted_model protocol from the previous"
              "predicted, processed and docked atomic structure ")

        args = {
            'inputPredictedModel': structure2,
            'inputDockedPredictedModel': protProcessDockPrediction.outputPdb,
            'inputVolume': map,
            'resolution': 3.2
        }
        protRebuildDockPrediction = self.newProtocol(
            PhenixProtRebuildDockPredictedAlphaFold2Model, **args)
        protRebuildDockPrediction.setObjLabel('Rebuild\ndock prediction\nP69905\n')
        self.launchProtocol(protRebuildDockPrediction)
        self.assertTrue(os.path.exists(
            protRebuildDockPrediction.outputPdb.getFileName()))

    def testDDockBuildPrediction1(self):
        """ Test the protocol dock and rebuild alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure2()
        structure2 = protImportPDB.outputPdb
        self.assertTrue(structure2.getFileName())
        self.assertFalse(structure2.getVolume())

        args = {
                'inputPredictedModel': structure2,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nP69905\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

        print("Run phenix dock_and_rebuild_predicted_model protocol "
              "from the previous atomic structure ")

        # import map
        map = self._importVolume1()
        self.assertTrue(map.getFileName())

        args = {
                'inputPredictedModel': structure2,
                'inputVolume': map,
                'resolution': 3.2
                }
        protDockBuildPrediction = self.newProtocol(
            PhenixProtDockAndRebuildAlphaFold2Model, **args)
        protDockBuildPrediction.setObjLabel('Dock and rebuild\nprediction\nP69905\n')
        self.launchProtocol(protDockBuildPrediction)
        self.assertTrue(os.path.exists(
            protDockBuildPrediction.outputPdb.getFileName()))


class TestBProtProcessDockBuildPredictedAlphaFold(TestImportData):

    def testAProcessPrediction2(self):
        """ Test the protocol process alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure1()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        args = {
                'inputPredictedModel': structure1,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nQ9BXJ8\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

    def testBProcessDockPrediction2(self):
        """ Test the protocol process and dock alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure1()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        args = {
                'inputPredictedModel': structure1,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nQ9BXJ8\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

        print("Run phenix dock_predicted_model protocol from the previous"
              "predicted and processed atomic structure ")

        # import map
        map = self._importVolume2()
        self.assertTrue(map.getFileName())

        args = {
                'inputPredictedModel': structure1,
                'inputProcessedPredictedModel': protProcessPrediction.outputPdb,
                'inputVolume': map,
                'resolution': 3.24
                }
        protProcessDockPrediction = self.newProtocol(
            PhenixProtDockPredictedAlphaFold2Model, **args)
        protProcessDockPrediction.setObjLabel('Dock\nprediction processed\nQ9BXJ8\n')
        self.launchProtocol(protProcessDockPrediction)
        self.assertTrue(os.path.exists(
            protProcessDockPrediction.outputPdb.getFileName()))

    def testCProcessDockBuildPrediction2(self):
        """ Test the protocol process, dock and rebuild alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure1()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        args = {
                'inputPredictedModel': structure1,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nQ9BXJ8\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

        print("Run phenix dock_predicted_model protocol from the previous"
              "predicted and processed atomic structure ")

        # import map
        map = self._importVolume2()
        self.assertTrue(map.getFileName())

        args = {
                'inputPredictedModel': structure1,
                'inputProcessedPredictedModel': protProcessPrediction.outputPdb,
                'inputVolume': map,
                'resolution': 3.24
                }
        protProcessDockPrediction = self.newProtocol(
            PhenixProtDockPredictedAlphaFold2Model, **args)
        protProcessDockPrediction.setObjLabel('Dock\nprediction processed\nQ9BXJ8\n')
        self.launchProtocol(protProcessDockPrediction)
        self.assertTrue(os.path.exists(
            protProcessDockPrediction.outputPdb.getFileName()))

        print("Run phenix rebuild_predicted_model protocol from the previous"
              "predicted, processed and docked atomic structure ")

        args = {
            'inputPredictedModel': structure1,
            'inputDockedPredictedModel': protProcessDockPrediction.outputPdb,
            'inputVolume': map,
            'resolution': 3.24
        }
        protRebuildDockPrediction = self.newProtocol(
            PhenixProtRebuildDockPredictedAlphaFold2Model, **args)
        protRebuildDockPrediction.setObjLabel('Rebuild\ndock prediction\nQ9BXJ8\n')
        self.launchProtocol(protRebuildDockPrediction)
        self.assertTrue(os.path.exists(
            protRebuildDockPrediction.outputPdb.getFileName()))

    def testDDockBuildPrediction2(self):
        """ Test the protocol dock and rebuild alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure1()
        structure1 = protImportPDB.outputPdb
        self.assertTrue(structure1.getFileName())
        self.assertFalse(structure1.getVolume())

        args = {
                'inputPredictedModel': structure1,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nQ9BXJ8\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))

        print("Run phenix dock_and_rebuild_predicted_model protocol "
              "from the previous atomic structure ")

        # import map
        map = self._importVolume2()
        self.assertTrue(map.getFileName())

        args = {
                'inputPredictedModel': structure1,
                'inputVolume': map,
                'resolution': 3.24
                }
        protDockBuildPrediction = self.newProtocol(
            PhenixProtDockAndRebuildAlphaFold2Model, **args)
        protDockBuildPrediction.setObjLabel('Dock and rebuild\nprediction\nQ9BXJ8\n')
        self.launchProtocol(protDockBuildPrediction)
        self.assertTrue(os.path.exists(
            protDockBuildPrediction.outputPdb.getFileName()))

class TestCProtProcessDockBuildPredictedAlphaFold(TestImportData):

    def testCProcessPrediction1(self):
        """ Test the protocol process alpahafold2 predicted model
        """
        print("Run phenix process_predicted_model protocol from a imported "
              "predicted atomic structure ")

        # import PDB
        protImportPDB = self._importStructure3()
        structure3 = protImportPDB.outputPdb
        self.assertTrue(structure3.getFileName())
        self.assertFalse(structure3.getVolume())

        args = {
                'inputPredictedModel': structure3,
               }

        protProcessPrediction = self.newProtocol(
            PhenixProtProcessPredictedAlphaFold2Model, **args)
        protProcessPrediction.setObjLabel('Process prediction\nA0A08WSY6\n')
        self.launchProtocol(protProcessPrediction)
        self.assertTrue(os.path.exists(
            protProcessPrediction.outputPdb.getFileName()))
