# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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

import os

from pyworkflow import Config
from pyworkflow.object import String, Float, Integer
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam
from phenix.constants import SUPERPOSE, PHENIX_HOME
from pwem.convert.atom_struct import fromCIFToPDB, \
    fromPDBToCIF, fromCIFTommCIF, AtomicStructHandler

try:
    from pwem.objects import AtomStruct
except:
    from pwem.objects import PdbFile as AtomStruct
from phenix import Plugin


class PhenixProtRunSuperposePDBs(EMProtocol):
    """Superpose two PDBs so that they optimally match """
    _label = 'superpose pdbs'
    _program = ""
    # _version = VERSION_1_2

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructureFixed', PointerParam,
                      pointerClass="AtomStruct",
                      label='Fixed atomic structure', important=True,
                      help="The moving PDB will be aligned to the fixed one")
        form.addParam('inputStructureMoving', PointerParam,
                      pointerClass="AtomStruct",
                      label='Moving atomic structure',
                      help="PDBx/mmCIF to be aligned")

    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runSuperposePDBsStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------

    def runSuperposePDBsStep(self):
        args = os.path.abspath(self.inputStructureFixed.get().getFileName())
        args += " "
        args += os.path.abspath(self.inputStructureMoving.get().getFileName())
        cwd = os.getcwd() + "/" + self._getExtraPath()
        try:
            Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args,
                                    extraEnvDict=None, cwd=cwd)
        except:
            # This exception will run when using phenix v. 1.16 after running
            # real space refine the .cif file generated can not be recognized by
            # superpose pdbs program and an error is produced
            list_args = args.split()
            self._runChangingCifFormatSuperpose(list_args)

    def createOutputStep(self):
        fnPdb = os.path.basename(self.inputStructureMoving.get().getFileName())
        fnPdb = fnPdb.split('.')[0]
        for file in os.listdir(self._getExtraPath()):
            if file.endswith('_fitted.pdb'):
                os.rename(self._getExtraPath() + "/" + file,
                          self._getExtraPath() + "/" + fnPdb + "_fitted.pdb")
        pdb = AtomStruct()
        pdb.setFileName(self._getExtraPath(fnPdb + "_fitted.pdb"))
        if self.inputStructureFixed.get().getVolume() is not None:
            pdb.setVolume(self.inputStructureFixed.get().getVolume())
        self._defineOutputs(outputPdb=pdb)
        self._defineSourceRelation(self.inputStructureFixed.get(), pdb)
        self._defineSourceRelation(self.inputStructureMoving.get(), pdb)

        logFile = os.path.abspath(self._getLogsPath()) + "/run.stdout"
        self._parseLogFile(logFile)
        self._store()

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram(SUPERPOSE)
        if not os.path.exists(program):
            errors.append("Cannot find " + program)

            # If there is any error at this point it is related to config variables
            errors.append("Check configuration file: " +
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set PHENIX_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("PHENIX_HOME = %s" % Plugin.getVar(PHENIX_HOME))
                errors.append("SUPERPOSE = %s" % SUPERPOSE)

        return errors

    def _summary(self):
        summary = []
        try:
            summary.append("RMSD between fixed and moving atoms (start): " +
                           str(self.startRMSD))
            summary.append("RMSD between fixed and moving atoms (final): " +
                           str(self.finalRMSD))
        except:
            summary.append("RMSD not yet computed")
        summary.append(
            "http://www.phenix-online.org/documentation/superpose_pdbs.htm")
        summary.append("Peter Zwart, Pavel Afonine, Ralf W. Grosse-Kunstleve")

        return summary

    # --------------------------- UTILS functions --------------------------

    def _parseLogFile(self, logFile):
        with open(logFile) as f:
            line = f.readline()
            while line:
                words = line.strip().split()
                if len(words) > 1:
                    if (words[0] == 'RMSD' and words[1] == 'between' and
                            words[6] == '(start):'):
                        self.startRMSD = Float(words[7])
                    elif (words[0] == 'RMSD' and words[1] == 'between' and
                          words[6] == '(final):'):
                        self.finalRMSD = Float(words[7])
                line = f.readline()

    def _runChangingCifFormatSuperpose(self, list_args):
        cwd = os.getcwd() + "/" + self._getExtraPath()
        try:
            if list_args[0].endswith(".cif") and list_args[1].endswith(".cif"):
                try:
                    # upgrade cifs
                    list_args1 = []
                    for i in range(0, 2):
                        list_args1.append(fromCIFTommCIF(list_args[i], list_args[i]))
                    args1 = list_args1[0] + " " + list_args1[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args1,
                                            extraEnvDict=None, cwd=cwd)
                except:
                    # convert cifs to pdbs
                    list_args2 = []
                    for i in range(0, 2):
                        list_args2.append(fromCIFToPDB(
                            list_args[i], list_args[i].replace('.cif', '.pdb')))
                    args2 = list_args2[0] + " " + list_args2[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args2,
                                            extraEnvDict=None, cwd=cwd)
            elif list_args[0].endswith(".cif") and list_args[1].endswith(".pdb"):
                try:
                    # pdbs: convert cif to pdb
                    list_args1 = []
                    list_args1.append(fromCIFToPDB(
                        list_args[0], list_args[0].replace('.cif', '.pdb')))
                    args1 = list_args1[0] + " " + list_args[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args1,
                                            extraEnvDict=None, cwd=cwd)
                except:
                    try:
                        # cifs: convert pdb to cif
                        list_args2 = []
                        list_args2.append(fromPDBToCIF(
                            list_args[1], list_args[1].replace('.pdb', '.cif')))
                        args2 = list_args[0] + " " + list_args2[0]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args2,
                                                extraEnvDict=None, cwd=cwd)
                    except:
                        # upgrade cif
                        list_args3 = []
                        list_args0 = args2.split()
                        for i in range(0, 2):
                            list_args3[i].append(fromCIFTommCIF(
                                list_args0[i], list_args0[i]))
                        args3 = list_args3[0] + " " + list_args3[1]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE),
                                                args3, extraEnvDict=None, cwd=cwd)
            elif list_args[0].endswith(".pdb") and list_args[1].endswith(".cif"):
                try:
                    # pdbs: convert cif to pdb
                    list_args1 = []
                    list_args1.append(fromCIFToPDB(
                        list_args[1], list_args[1].replace('.cif', '.pdb')))
                    args1 = list_args[0] + " " + list_args1[0]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args1,
                                            extraEnvDict=None, cwd=cwd)
                except:
                    try:
                        # cifs: convert pdb to cif
                        list_args2 = []
                        list_args2.append(fromPDBToCIF(
                            list_args[0], list_args[0].replace('.pdb', '.cif')))
                        args2 = list_args2[0] + " " + list_args[1]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE), args2,
                                                extraEnvDict=None, cwd=cwd)
                    except:
                        # upgrade cifs
                        list_args3 = []
                        list_args0 = args2.split()
                        for i in range(0, 2):
                            list_args3.append(fromCIFTommCIF(
                                list_args0[i], list_args0[i]))
                        args3 = list_args3[0] + " " + list_args3[1]
                        Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE),
                                                args3, extraEnvDict=None, cwd=cwd)
        except:
            # biopython conversion
            aSH = AtomicStructHandler()
            try:
                for i in range(0, 2):
                    aSH.read(list_args[i])
                    aSH.write(list_args[i])
                    args = list_args[0] + " " + list_args[1]
                    Plugin.runPhenixProgram(Plugin.getProgram(SUPERPOSE),
                                            args, extraEnvDict=None, cwd=cwd)
            except:
                print("CIF file standarization failed.")


