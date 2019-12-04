
import os

from pyworkflow import Config
from pyworkflow.object import Float, Integer
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import (PointerParam, FloatParam, \
    StringParam)
from phenix.constants import (PHENIXVERSION)
from pwem.convert.headers import Ccp4Header
from phenix import Plugin
from pyworkflow.protocol.constants import LEVEL_ADVANCED
import collections
import json



class PhenixProtRunRefinementBase(EMProtocol):
    """MolProbity is a Phenix application to validate the geometry of an
atomic structure derived from a cryo-EM density map.
"""
    _label = 'refinementBase'
    _program = ""

    MOLPROBITYOUTFILENAME = 'molprobity.out'
    MOLPROBITYCOOTFILENAME = 'molprobity_coot.py'
    MOLPROBITYPKLFILENAME = 'molprobity.pkl'
    SUMMARYFILENAME = 'validationSummary.txt'


    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addParallelSection(threads=1, mpi=0)
        form.addSection(label='Input')
        # TODO: Consider the possibility of including .mtz files in Scipion
        # TODO: input volume (read the next comment)
        # Real-space refine map volume: It is possible to use either CCP4 type
        # .map files (which are three-dimensional grid of voxels with electron
        # density values) or .mtz files (which contain Fourier coefficients
        # for an electron density map). If CCP4 type map files are used,
        # the resolution is needed as well. If mtz files are used,
        # the resolution is detected automatically from the map coefficients.
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="")
        form.addParam('resolution', FloatParam,
                      label='Resolution (A):',
                      default=3.0,
                      help='Set the resolution of the input volume.')
        form.addParam('inputStructure', PointerParam,
                      pointerClass="AtomStruct", allowsNull=False,
                      label='Input atomic structure.',
                      help="Set the atomic structure to be processed.\n"
                           "Supported formats are PDB or mmCIF; this last one"
                           " is especially useful for very large structures.")
        form.addParam('extraParams', StringParam,
                       label="Extra Params ",
                       default="",
                       expertLevel=LEVEL_ADVANCED,
                       help="This string will be added to the phenix command.\n"
                            "Syntax: paramName1=value1 paramName2=value2 ")

    # --------------------------- STEPS functions --------------------------

    def convertInputStep(self, tmpMapFileName):
        """ convert 3D maps to MRC '.mrc' format
        """
        vol = self._getInputVolume()
        inVolName = vol.getFileName()
        newFn = self._getExtraPath(tmpMapFileName)
        origin = vol.getOrigin(force=True).getShifts()
        sampling = vol.getSamplingRate()
        Ccp4Header.fixFile(inVolName, newFn, origin, sampling, Ccp4Header.START)  # ORIGIN


    def _summary(self):
        #  Think on how to update this summary with created PDB
        summary = []
        try:
            summary.append("MolProbity statistics:\n")
            summary.append("Ramachandran outliers: %0.2f %%   (Goal: < 0.2%%)  "\
                           "Ramachandran favored: %0.2f %%   (Goal: > 98%%) " %\
                           (self.ramachandranOutliers.get(),
                            self.ramachandranFavored.get())
                           )
            summary.append("Rotamer outliers:           %0.2f"\
                           " %%   (Goal: < 1%%)    "\
                           "C-beta outliers:              %d"\
                           "            (Goal: 0) " % \
                            (self.rotamerOutliers.get(),
                             self.cbetaOutliers.get()
                            ))
            summary.append("Clashscore:                   %0.2f"\
                           "                             Overall "\
                           "score:                %0.2f" %\
                           (self.clashscore.get(), self.overallScore.get())
                           )
        except:
            summary = ["Overall score not yet computed"]
        return summary

    def validateBase(self, program, label):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram(program)
        if program is None:
            errors.append("Missing variables %s and/or PHENIX_HOME" % label)
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)
            errors.append("Check if you need to upgrade your PHENIX version"
                          " to run %s.\nYour current PHENIX version is %s.\n"
                          % (label, PHENIXVERSION))

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: "+
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set %s and PHENIX_HOME variables "
                          "properly."% label)
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % Plugin.getHome())
                #errors.append("%s = %s" % label, program)

        return errors

    # --------------------------- UTILS functions --------------------------

    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputStructure.get().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol

    def _parseFile(self, fileName):
        with open(fileName, encoding = "ISO-8859-1") as f:
            line = f.readline()
            while line:
                words = line.strip().split()
                if len(words) > 1:
                    if (words[0] == 'Ramachandran' and words[1] == 'outliers'):
                        self.ramachandranOutliers = Float(words[3])
                    elif (words[0] == 'favored' and words[1] == '='):
                        self.ramachandranFavored = Float(words[2])
                    elif (words[0] == 'Rotamer' and words[1] == 'outliers'):
                        self.rotamerOutliers = Float(words[3])
                    elif (words[0] == 'C-beta' and words[1] == 'deviations'):
                        self.cbetaOutliers = Integer(words[3])
                    elif (words[0] == 'Clashscore' and words[1] == '='):
                        self.clashscore = Float(words[2])
                    elif (words[0] == 'MolProbity' and words[1] == 'score'):
                        self.overallScore = Float(words[3])
                line = f.readline()


    def _readValidationPklFile(self, fileName):
        self.SUMMARYFILENAME = self._getTmpPath(self.SUMMARYFILENAME)
        command = """import pickle
import collections
import json

def pickleData(file):
    with open(file,"r") as f:
        return pickle.load(f)

# process file {VALIDATIONCRYOEMPKLFILENAME}"
data = pickleData('{VALIDATIONCRYOEMPKLFILENAME}')
dictSummary = collections.OrderedDict()

dictSummary['Rhama_Outliers'] = data.model.geometry.ramachandran.outliers
dictSummary['Rhama_Favored'] = data.model.geometry.ramachandran.favored
dictSummary['Rota_Outliers'] = data.model.geometry.rotamer.outliers
dictSummary['Cbeta_Outliers_n'] = data.model.geometry.c_beta.cbetadev.n_outliers
dictSummary['Clash_score'] = data.model.geometry.clash.score
dictSummary['MolProbity_score'] = data.model.geometry.molprobity_score
""".format(VALIDATIONCRYOEMPKLFILENAME=fileName)

        command += """with open('%s',"w") as f:
    f.write(json.dumps(dictSummary))
""" % (self.SUMMARYFILENAME)

        pythonFileName = self.SUMMARYFILENAME.replace('.txt', '.py')
        # write script file
        with open(pythonFileName, "w") as f:
            f.write(command)

        # execute file with phenix.python
        Plugin.runPhenixProgram("", pythonFileName)

        # read file in scipion python
        with open(self.SUMMARYFILENAME, "r") as f:
            dictSummary = f.read()

        dictSummary = json.loads(
            dictSummary, object_pairs_hook=collections.OrderedDict)

        self.ramachandranOutliers = Float(dictSummary['Rhama_Outliers'])
        self.ramachandranFavored = Float(dictSummary['Rhama_Favored'])
        self.rotamerOutliers = Float(dictSummary['Rota_Outliers'])
        self.cbetaOutliers = Integer(dictSummary['Cbeta_Outliers_n'])
        self.clashscore = Float(dictSummary['Clash_score'])
        self.overallScore = Float(dictSummary['MolProbity_score'])

    def _writeArgsMolProbity(self, atomStruct, vol=None):
        args = ""
        args += atomStruct
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            args += " "
            args += " map_file_name=%s" % vol
            args += " pickle=True"
            args += " "
            args += " d_min=%f" % self.resolution.get()
        args += " "
        numberOfThreads = self.numberOfThreads.get()
        if numberOfThreads > 1:
            args += " nproc=%d" % numberOfThreads
        return args

    def _writeArgsValCryoEM(self, atomStruct, volume, vol):
        args = " " + atomStruct
        args += " " + volume
        if vol.getHalfMaps():
            halves = []
            for halfMap in vol.getHalfMaps().split(','):
                if not os.path.abspath(halfMap).endswith(".mrc"):
                    half = os.path.abspath(halfMap).split(".")[0] + ".mrc"
                else:
                    half = os.path.abspath(halfMap)
                halves.append(half)
            args += " " + halves[0] + " " + halves[1]
        args += " " + ("resolution=%f" % self.resolution.get())
        args += " pickle=True"
        args += " slim=False"
        args += " pdb_interpretation.clash_guard.nonbonded_distance_threshold=None"
        args += " %s " % self.extraParams.get()
        return args
