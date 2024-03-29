# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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
# *
# **************************************************************************

import glob
import json
import os

from pyworkflow import Config

from phenix.constants import EMRINGER, DISPLAY
from pwem.convert.headers import Ccp4Header
from pwem.protocols import EMProtocol
from pyworkflow.object import String
from pyworkflow.protocol.params import BooleanParam, PointerParam
from phenix import Plugin
from phenix.constants import PHENIX_HOME
from pwem.convert.atom_struct import retry

class PhenixProtRunEMRinger(EMProtocol):
    """EMRinger is a Phenix application to validate the agreement between
the initial map and the derived low-resolution atomic structure. This program
samples the density around Chi1 angles of protein sidechains. Electronic
density and appropriate rotameric angles must coincide for each residue if
the atomic structure backbone has been perfectly fitted to the map.
"""
    _label = 'emringer'
    _program = ""
    # _version = VERSION_1_2
    EMRINGERFILE = 'emringer.map'

    EMRINGERTRANSFERFILENAME = 'emringer_transfer.txt'
    EMRINGERSCORESFILENAME = 'emringer_scores.pkl'
    EMRINGEROPTIMALTHRESHOLDFILENAME = 'thresholds.pkl'
    EMRINGERROTAMERRATIOSFILENAME = 'rotamer_ratios.pkl'
    EMRINGERMAXZSCOREFILENAME = 'zscores.pkl'

    # --------------------------- DEFINE param functions -------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Set the starting volume")
        form.addParam('inputStructure', PointerParam,
                      pointerClass="AtomStruct",
                      label='Input atomic structure',
                      help="PDBx/mmCIF to be validated against the volume. ")
        form.addParam('doTest', BooleanParam, default=False,
                      label='Test', condition='False',
                      help="""Saves the temporary file 
                      'emringer_transfer.py' in the extra folder. Use for 
                      testing""")

    # --------------------------- INSERT steps functions ---------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runEMRingerStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def convertInputStep(self):
        """ convert 3D maps to MRC '.mrc' format
        """
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        newFn = self._getExtraPath(self.EMRINGERFILE)
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        Ccp4Header.fixFile(inVolName, newFn, origin, sampling, Ccp4Header.START)  # ORIGIN

    def runEMRingerStep(self):
        inPDBAtomStructFile = self.inputStructure.get().getFileName()
        # atomStruct = os.path.abspath(inPDBAtomStructFile)
        if not inPDBAtomStructFile == os.path.abspath(inPDBAtomStructFile):
            atomStruct = os.getcwd() + "/" + inPDBAtomStructFile
        else:
            atomStruct = inPDBAtomStructFile
        # vol = os.path.abspath(self._getExtraPath(self.EMRINGERFILE))
        vol = os.getcwd() + "/" + self._getExtraPath(self.EMRINGERFILE)
        args = self._writeArgsEMRinger(atomStruct, vol)
        # script with auxiliary files
        import sys
        sys.stdout.flush()
        # import time
        # time.sleep(30)
        retry(Plugin.runPhenixProgram, Plugin.getProgram(EMRINGER),
              args, cwd=self._getExtraPath(),
              listAtomStruct=[atomStruct], log=self._log,
              messages=[("max_max = max(maxima)", 
                        "No features could be detected in the density around the model, so EMRinger can't proceed."), 
                        ("Sorry: No residues could be scanned by EMRinger, so scores cannot be generated.",
                        "Sorry: No residues could be scanned by EMRinger, so scores cannot be generated."),
                        ("Sorry: No features could be detected in the density around the model, so EMRinger can't",
                        "Sorry: No features could be detected in the density around the model, so EMRinger can't")], 
                        sdterrLog = self.getLogsLastLines)

    def createOutputStep(self):
        # get emringer information
        # values are stored in python pickle files
        # string to be run by phenix python

        # temporary file with emringer values
        test = self.doTest
        if test == True:
            EMRINGERTRANSFERFILENAME = self._getExtraPath(
                self.EMRINGERTRANSFERFILENAME)
        else:
            EMRINGERTRANSFERFILENAME = self._getTmpPath(
                self.EMRINGERTRANSFERFILENAME)
        # directory with files
        plots = glob.glob(self._getExtraPath("*_plots"))[0]

        # emringer main file
        mainDataFile = glob.glob(self._getExtraPath("*_emringer.pkl"))[0]

        # file with emringer score
        EMRINGERSCORESFILENAME = os.path.join(
            plots, self.EMRINGERSCORESFILENAME)

        EMRINGERFILES = []
        keyList = []
        # files with other scores
        EMRINGEROPTIMALTHRESHOLDFILENAME = \
            os.path.join(plots, self.EMRINGEROPTIMALTHRESHOLDFILENAME)
        EMRINGERFILES.append(EMRINGEROPTIMALTHRESHOLDFILENAME)
        keyList.append('Optimal Threshold')

        EMRINGERROTAMERRATIOSFILENAME = \
            os.path.join(plots, self.EMRINGERROTAMERRATIOSFILENAME)
        EMRINGERFILES.append(EMRINGERROTAMERRATIOSFILENAME)
        keyList.append('Rotamer-Ratio')

        EMRINGERMAXZSCOREFILENAME = \
            os.path.join(plots, self.EMRINGERMAXZSCOREFILENAME)
        EMRINGERFILES.append(EMRINGERMAXZSCOREFILENAME)
        keyList.append('Max Zscore')

        # I do not want to use the local pickle version
        # but the version used to create the files
        # that is the phenix.python version

        command = """import pickle
import json
import sys
import collections
from mmtbx.ringer.em_scoring import parse_pickle
from mmtbx.ringer.em_rolling import easy_pickle, RingerDict

def pickleData(file):
    with open(file,"r") as f:
        return pickle.load(f)

# process file %s"
data = pickleData('%s')
dataDict = collections.OrderedDict()
# index to maxscore
maxScore = max(data)
maxScoreIndex = data.index(maxScore)
""" % (EMRINGERSCORESFILENAME, EMRINGERSCORESFILENAME)

        for key, fileName in zip(keyList, EMRINGERFILES):
            command += "# process file %s\n" % fileName
            command += "data = pickleData('%s')\n" % fileName
            command += "dataDict['%s'] =  data[maxScoreIndex]\n" % key

        command += """file_name='%s'
waves, thresholds = parse_pickle(file_name, out=sys.stdout)
# waves: list of ringer_residue objects (Non-gamma-branched, non-proline
# aminoacids with a non-H gamma atom) used in global EMRinger score
# computation
dataDict['Model Length'] = len(waves)
dataDict['EMRinger Score'] = maxScore  # find max score
# dataDict is shown in a table
# except for those attributes starting with "_"
dataDict['_thresholds'] = thresholds
dataDict['_maxScoreIndex'] = maxScoreIndex
ringer_results = easy_pickle.load(file_name)
# ringer_results: list of ringer_residue objects (aminoacids with gamma
# carbon; at least one Chi angle)
hierarchy = RingerDict(ringer_results, 0)
chains = sorted(hierarchy.get_chains())
dataDict['_chains'] = chains

# Aminoacids that contain non-H gamma atom (at least one Chi angle; included in
# ringer_results)
dictResidue = {}
formatResidue = []
for residue in ringer_results:
    resFormatTmp = residue.format()
    formatResidue.append(resFormatTmp)
    dictResidue[resFormatTmp] = ringer_results.index(residue)
dataDict['_residues_format'] = formatResidue
dataDict['_residues_dict'] = dictResidue
""" % mainDataFile
        command += """with open('%s',"w") as f:
    f.write(json.dumps(dataDict))
""" % (EMRINGERTRANSFERFILENAME)

        pythonFileName = EMRINGERTRANSFERFILENAME.replace('.txt', '.py')
        # write script file
        with open(pythonFileName, "w") as f:
            f.write(command)

        # execute file with phenix.python
        Plugin.runPhenixProgram("", pythonFileName)

        # read file in scipion python
        with open(EMRINGERTRANSFERFILENAME, "r") as f:
            self.stringDataDict = String(f.read())
            # self.dataDict = json.loads(f.read())

        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram(EMRINGER)
        if program is None:
            errors.append("Missing variables EMRINGER and/or PHENIX_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)
        elif not self.is_tool(DISPLAY):
            errors.append("display program missing.\n "
                          "Please install imagemagick package")

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: " +
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set EMRINGER and PHENIX_HOME variables "
                          "properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % Plugin.getVar(PHENIX_HOME))
                errors.append("EMRINGER = %s" % EMRINGER)

        # Check that the input volume exist
        if self._getInputVolume() is None:
            errors.append("Error: You should provide a volume.\n")

        return errors

    def _summary(self):
        #  Think on how to update this summary with created PDB
        summary = []
        try:
            dataDict = json.loads(str(self.stringDataDict))
            summary.append("Optimal Threshold: %0.2f   Rotamer-Ratio: %0.2f"
                           % (dataDict['Optimal Threshold'],
                              dataDict['Rotamer-Ratio']))
            summary.append("Max Zscore:         %0.2f  Model Length: %d"
                           % (dataDict['Max Zscore'],
                              dataDict['Model Length']))
            summary.append("EMRinger Score:  %0.2f"
                           % dataDict['EMRinger Score'])
        except:
            summary = ["EMRinger Score not yet computed"]

        return summary

    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("TODO")

        return methodsMsgs

    def _citations(self):
        return ['Barad_2015']

    # --------------------------- UTILS functions --------------------------

    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputStructure.get().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol

    def is_tool(self, name):
        """Check whether `name` is on PATH."""
        from shutil import which
        return which(name) is not None

    def _writeArgsEMRinger(self, atomStruct, vol):
        args = " "
        args += atomStruct
        args += " "
        args += vol
        return args

