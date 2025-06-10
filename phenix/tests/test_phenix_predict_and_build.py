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
from phenix.protocols import PhenixPredictAndBuildCryoEM  \
    
from pwem.protocols.protocol_import import (ProtImportSequence,
                                            ProtImportVolumes)
from pwem.objects.data import Alphabet

from pyworkflow.tests import *
import pwem.protocols as emprot
import requests
import xml.etree.ElementTree as ET


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import maps and sequencesof chains alpha and beta from haemoglobin 5ni1
    """
    def _importVolume(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/emd_3488.map'),
                'samplingRate': 1.05,
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume emd_3488.map\n')
        self.launchProtocol(protImportVol)
        volume1 = protImportVol.outputVolume
        return volume1

    def _seqChainAlpha(self):
        """Download seq chain alpha from UniProt
        """
        args = {'inputSequenceName': '5ni1_alpha_chain',
                'inputProteinSequence':
                    emprot.ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': 'P69905'
                }
        protImportSequence = self.newProtocol(emprot.ProtImportSequence, **args)
        protImportSequence.setObjLabel('5ni1_alphaChain')
        self.launchProtocol(protImportSequence)
        sequence_5ni1_A = protImportSequence.outputSequence
        self.assertEqual("P69905", protImportSequence._getUniProtID())
        return sequence_5ni1_A
    
    def _seqChainBeta(self):
        """Download seq chain beta from UniProt
        """
        args = {'inputSequenceName': '5ni1_beta_chain',
                'inputProteinSequence':
                    emprot.ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': 'P68871'
                }
        protImportSequence = self.newProtocol(emprot.ProtImportSequence, **args)
        protImportSequence.setObjLabel('5ni1_betaChain')
        self.launchProtocol(protImportSequence)
        sequence_5ni1_B = protImportSequence.outputSequence
        self.assertEqual("P68871", protImportSequence._getUniProtID())
        return sequence_5ni1_B


class TestPredictAndBuildCryoEM(TestImportData):


    def testAPredictAndBuildCryoEM(self):
        """ Test the protocol run phenix.predict_and_build method with one full map 
        and one sequence
        """
        print("Run phenix predict_and_build protocol from an imported map"
              " and the sequence imported from an atom structure ")
        
        # import full map
        volume1 = self._importVolume()
        self.assertTrue(volume1.getFileName())

        # import sequence chain alpha
        sequence_5ni1_A = self._seqChainAlpha()
        sequence_5ni1_B = self._seqChainBeta()

        args = {
                'inputVolume': volume1,
                'inputSequenceS': [sequence_5ni1_A, sequence_5ni1_B,sequence_5ni1_A, sequence_5ni1_B],
                'resolution': 3.0
               }

        protPredictAndBuildCryoEM = self.newProtocol(PhenixPredictAndBuildCryoEM, **args)
        protPredictAndBuildCryoEM.setObjLabel('Predict and build\n5ni1_chainA_chainB\n')
        self.launchProtocol(protPredictAndBuildCryoEM)
        # self.assertTrue(os.path.exists(
        #     protPredictAndBuildCryoEM.outputPdb.getFileName()))

        # self.assertTrue(os.path.exists(
        #     protDockInMap.outputPdb.getFileName()))