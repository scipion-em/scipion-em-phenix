# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

import os
import re
import sqlite3

from pwem.objects import AtomStruct
from pyworkflow.protocol.params import (StringParam,  IntParam,
                                        PointerParam, BooleanParam)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .protocol_refinement_base import PhenixProtRunRefinementBase
from pwem.convert.atom_struct import AtomicStructHandler
from pyworkflow.protocol import STEPS_PARALLEL
# from pwem.emlib.image import ImageHandler
from phenix import Plugin
from ccp4 import Plugin as PluginCCP4
from ccp4.convert import (runCCP4Program)
from ccp4.constants import CCP4_BINARIES
from phenix.constants import (REALSPACEREFINE)
from pwem.convert.atom_struct import retry

COOT = CCP4_BINARIES['COOT']
COOTSCRIPTFILENAME = "cootScript.py"
COOTPDBTEMPLATEFILENAME = "coot_%06d_Imol_%04d_version_%04d.pdb" # protId, modelID, counter
COOTPDBTEMPLATEFILENAMEINV = "coot_%06d_Imol_%04d_version_%04d_inv.pdb" # protId, modelID, counter
DATAFILE = 'db.sqlite3'
TABLE = 'myTable'

class PhenixProtSearchFit(PhenixProtRunRefinementBase):
    """given a chain of n alanines, a 3D map and
    a sequence search for the subsequence of n aminoacids
    that better fits in the density. Only works if the
    atomic structure has a single chain
    """
    _label = 'search fit'
    _program = ""
    # _version = VERSION_1_2
    FITTEDFILE = 'fitted.mrc'

    def __init__(self, **kwargs):
        super(PhenixProtSearchFit, self).__init__(**kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    # --------------------------- DEFINE param functions -------------------
    def _defineParams(self, form):
        super(PhenixProtSearchFit, self)._defineParams(form)
        param = form.getParam('numberOfThreads')
        param.default.set(1)
        param = form.getParam('inputStructure')
        param.help.set('Alanine chain used as template')
        form.addParam('inputSequence', PointerParam, pointerClass="Sequence",
                         label='Test sequence', important=True,
                         help="Input the aminoacid sequence to fit with the "
                              "alanine chain")
        form.addParam("firstaa", IntParam, important=True, default = -1,
                      label='First aa',
                      help='First aminoacid in the sequence to be considered,'
                           ' -1 -> first aa (1-> first aminoacid also)')
        form.addParam("lastaa", IntParam, important=True, default = -1,
                      label='Last aa',
                      help='Last aminoacid in the sequence to be considered,'
                           ' -1 -> first aa (1-> first aminoacid also)')
        form.addParam('extraCommands', StringParam,
                       label="Extra Params ",
                       default="",
                       expertLevel=LEVEL_ADVANCED,
                       help="This string will be added to the coot\n"
                            " script")
        # real space refine
        form.addParam("doSecondary", BooleanParam, label="Secondary structure",
                      default=True, expertLevel=LEVEL_ADVANCED,
                      help="Set to TRUE to use secondary structure "
                           "restraints.\nOnly for PHENIX versions higher than 1.13.")
        form.addParam("macroCycles", IntParam, label="Macro cycles",
                      default=5, expertLevel=LEVEL_ADVANCED,
                      help="Number of iterations of refinement.\nAlthough 5 "
                           "macro-cycles is usually sufficient, in cases in "
                           "which model geometry or/and model-to-map fit is "
                           "poor the use of more macro-cycles could be "
                           "helpful.\n")
        group = form.addGroup('Optimization strategy options')
        group.addParam('minimizationGlobal', BooleanParam,
                       label="Global minimization: ", default=True,
                       expertLevel=LEVEL_ADVANCED,
                       help="Phenix default parameter to look for the global "
                            "minimum of the model.\nGenerally, refinement "
                            "with all defaults is "
                            "sufficient.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('rigidBody', BooleanParam,
                       label="Rigid body: ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Refinement strategy that considers groups of "
                            "atoms that move (rotate and translate) as a "
                            "single body.\n")
        group.addParam('localGridSearch', BooleanParam,
                       label="Local grid search: ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Refinement strategy that considers "
                            "local rotamer fitting.\n\n Generally, refinement "
                            "with all defaults is sufficient.\n Including "
                            "local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('morphing', BooleanParam,
                       label="Morphing ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Morphing procedure distorts a model to match an "
                            "electron density map.\n\nGenerally, refinement "
                            "with all defaults is "
                            "sufficient.\n Including local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('simulatedAnnealing', BooleanParam,
                       label="Simulated annealing ", default=False,
                       expertLevel=LEVEL_ADVANCED,
                       help="Optimization technique known as molecular "
                            "dynamics refinement; it minimizes the energy of "
                            "the model.\n"
                            "Generally, refinement with all defaults is "
                            "sufficient.\n Including local fitting, morphing, "
                            "or simulated annealing "
                            "( local_grid_search+morphing+simulated_annealing) "
                            "into refinement may significantly increase "
                            "runtime.\nOther options "
                            "of use: run=minimization_global+local_grid_search"
                            "+morphing+simulated_annealing\n")
        group.addParam('adp', BooleanParam,
                       label="Atomic Displacement Parameters (ADPs) ",
                       default=True,
                       expertLevel=LEVEL_ADVANCED,
                       help="Phenix default parameter.\nGenerally, refinement "
                            "with all defaults is sufficient.\n\nADP ("
                            "B-factors) refinement against the map is "
                            "performed at the last macro-cycle only. ")
        form.addParallelSection(threads=0, mpi=1)


    # --------------------------- INSERT steps functions ---------------
    def _insertAllSteps(self):
        # compute alanine atom struct len
        inputPdb = self.inputStructure.get().getFileName()
        atomStruct = AtomicStructHandler(inputPdb)
        # we assume that there is a single model and a single chain
        atomStructSize = sum(1 for _ in atomStruct.getStructure().get_residues())
        chainName = next(atomStruct.getStructure().get_chains()).get_id()
        firstAAinChain =next(atomStruct.getStructure().get_residues()).id[1]


        # starting and ending aa
        firstaa = self.firstaa.get()
        lastaa = self.lastaa.get()

        # compute sequence size
        sequenceSize = len(self.inputSequence.get().getSequence())
        if lastaa == -1:
            lastaa = sequenceSize
        if firstaa == -1:
            firstaa = 1
        numberOfSteps = lastaa - firstaa + 1

        # steps
        prepareId = self._insertFunctionStep('convertInputStep', self.FITTEDFILE)
        # mutateChain
        mutateId = self._insertFunctionStep('mutateStep',
                                            firstaa,  # in seq
                                            firstAAinChain, # in struct
                                            atomStructSize,
                                            chainName,
                                            numberOfSteps,
                                            prerequisites=[prepareId])
        refineIdList=[]
        numberOfThreads = self.numberOfMpi.get()
        for start in range(numberOfThreads):
            refineId = self._insertFunctionStep('refineStep2',
                                                prerequisites=[mutateId])
            refineIdList.append(refineId)

        self._insertFunctionStep('createOutputStep', prerequisites=refineIdList)

    # --------------------------- STEPS functions --------------------------
    def createTable(self):
        """ Create table and clean it if needed"""
        conn = sqlite3.connect(os.path.abspath(self._getExtraPath(DATAFILE)))
        c = conn.cursor()
        create_table_sql = """CREATE TABLE IF NOT EXISTS %s (   id INTEGER PRIMARY KEY AUTOINCREMENT,
                                                           filename TEXT, 
                                                           done int DEFAULT 0,
                                                           model_to_map_fit float DEFAULT -1
                                                           )""" % TABLE
        print("create_table_sql", create_table_sql)
        c.execute(create_table_sql)
        c.close()
        conn.close()

    def mutateStep(self, firstaa, firstAAinChain,
                   atomStructSize,
                   chainName,
                   numberOfSteps):
        """ mutate atom struct inputStructure using
           aa in sequence  inputSequence starting at firstaa"""
        scriptFile = self._getExtraPath(COOTSCRIPTFILENAME)
        f = open(scriptFile, "w")
        fnAtomStruct = self.inputStructure.get().getFileName()
        f.write("# read atom structure (pdb) file\n")
        f.write("read_pdb('%s')\n" % os.path.abspath(fnAtomStruct))
        f.write("# mutation loop\n")
        database = os.path.abspath(self._getExtraPath(DATAFILE))
        f.write("import sqlite3\n")
        f.write("conn = sqlite3.connect('%s')\n" % database)
        f.write("cur = conn.cursor()\n")
        iMol = 0 # pdb id is  0 or 1 for inverse
        startMut = firstAAinChain  # 1
        endMut = firstAAinChain + atomStructSize -1
        self.createTable()

        for start in range(numberOfSteps):
            seq = self.inputSequence.get().getSequence()[firstaa + start : firstaa + start + atomStructSize]
            f.write("mutate_residue_range(%d, '%s', %d, %d, '%s')\n" % (iMol,
                                                                      chainName,
                                                                      startMut,
                                                                      endMut,
                                                                      seq))
            outFileName = self._getExtraPath(COOTPDBTEMPLATEFILENAME% (0,0,start))
            f.write("save_coordinates(0, '%s')\n" % outFileName)
            command = "INSERT INTO %s(filename) VALUES('%s')" % (TABLE,
                                                                 os.path.abspath(outFileName)
                                                                 )
            f.write('cur.execute("%s")\n' % command)
        # invert sequence
        ##iMol = 1 # pdb id is  0 or 1 for inverse
        ##f.write('reverse_direction_of_fragment(0, "%s", 1)\n'% chainName)
        ##f.write('db_mainchain(0, "%s", %d, %d,"forwards")\n' % (chainName, startMut, endMut))
        ##for start in range(numberOfSteps):
        ##    seq = self.inputSequence.get().getSequence()[firstaa + start : firstaa + start + atomStructSize]
        ##    f.write("mutate_residue_range(%d, '%s', %d, %d, '%s')\n" % (iMol,
        ##                                                              chainName,
        ##                                                              startMut,
        ##                                                              endMut,
        ##                                                              seq))
        ##    outFileName = self._getExtraPath(COOTPDBTEMPLATEFILENAMEINV% (0,0,start))
        ##    f.write("save_coordinates(1, '%s')\n" % outFileName)
        ##    command = "INSERT INTO %s(filename) VALUES('%s')" % (TABLE,
        ##                                                         os.path.abspath(outFileName)
        ##                                                         )
        ##    f.write('cur.execute("%s")\n' % command)
        if len(self.extraCommands.get()) > 0:
            f.write("\n#Extra Commands\n")
            f.write("%\n" % self.extraCommands.get())
        f.write("conn.commit()\n")
        f.write("cur.close()\n")
        f.write("conn.close()\n")
        f.write("exit(0)\n")
        f.close()

        args = ""
        args += " --no-graphics "
        args += " -s %s" % scriptFile
        # run
        self._log.info('Launching: ' + PluginCCP4.getProgram(COOT) + ' ' + args)
        runCCP4Program(PluginCCP4.getProgram(COOT), args)
    
    def extractNumber(self,filename):
        """ returns None if log file does not exist
            or the log file has not store the overall CC_mask value,
            otherwose returns the CC value"""
        if os.path.exists(filename):
            return None
        myfile = open(filename, "rt")
        contents = myfile.read()
        myfile.close()
        m = re.findall(
                r"verall(.+)\n(\*+)\nmodel-to-map fit, CC_mask: (\d+.\d+)", contents)[-1]
        # [-1] --> find last ocurrence
        if not m:
            return None
        return m[2]
      
    def refineStep2(self):
        # atomStruct = os.path.abspath(self.inputStructure.get().getFileName())
        ## vol = os.path.abspath(self._getInputVolume().getFileName())
        vol = os.path.abspath(self._getExtraPath(self.FITTEDFILE))
        cwd = os.getcwd() + "/" + self._getExtraPath()
        conn = sqlite3.connect(os.path.abspath(self._getExtraPath(DATAFILE)))
        c = conn.cursor()
        while 1:
            command = """SELECT filename from %s 
                         WHERE done=0 
                         LIMIT 1""" % TABLE
            c.execute(command)
            myrow = c.fetchone()
            if not myrow:
                #c.close()
                #conn.close()
                print("refineStep: no more available works")
                break  # break while if no job is available

            atomStructFn, = myrow
            accepted = c.execute("""UPDATE %s 
                                    SET done=1 
                                    WHERE filename='%s' AND done=0""" % (TABLE, atomStructFn))
            # ^^^^^^ This will return the number of rows updated.
            # Note that we only update if done is still '0', so if we get 1 updated
            # row, we're sure no one else took our job. This works because UPDATE is atomic.

            if not accepted:
                # Whoops this job was taken! Try again and get another one
                continue
            conn.commit()

            args = self._writeArgsRSR(atomStructFn, vol)
            # do not retry just execute the phenix program
            retry(Plugin.runPhenixProgram,
                 Plugin.getProgram(REALSPACEREFINE), args,
                 # cwd=os.path.abspath(self._getExtraPath()),
                 cwd=cwd,
                 listAtomStruct=[atomStructFn], log=self._log)
#            refinedFile = False
#            for item in os.listdir(self._getExtraPath()):
#                if item.endswith("_real_space_refined.cif"):
#                    refinedFile = True
#                    break
#            if refinedFile == False:
#                print("WARNING!!!\nPHENIX error:\n pdb_interpretation.clash_guard" \
#                      " failure: High number of nonbonded interaction distances " \
#                      "< 0.5. This error has been disable by running the same " \
#                      "command with the same following additional " \
#                      "argument:\npdb_interpretation.clash_guard." \
#                      "nonbonded_distance_threshold=None ")
#                args += " "
#                args += "pdb_interpretation.clash_guard." \
#                        "nonbonded_distance_threshold=None"
#                retry(Plugin.runPhenixProgram,
#                      Plugin.getProgram(REALSPACEREFINE), args,
#                      # cwd=os.path.abspath(self._getExtraPath()),
#                      cwd=cwd,
#                      listAtomStruct=[atomStructFn], log=self._log)
            logFileFn = atomStructFn[:-4] + "_real_space_refined.log"
            model_to_map_fit = self.extractNumber(logFileFn)
            # check model_to_map_fit, if it is None we need to retry
            # otherwise just continue.
            print("model_to_map_fit", model_to_map_fit)
            accepted = c.execute("""UPDATE %s 
                                    SET model_to_map_fit=%f 
                                    WHERE filename='%s'""" % (TABLE,
                                                              float(model_to_map_fit),
                                                              atomStructFn)
                                 )
            conn.commit()
        c.close()
        conn.close()

    def createOutputStep(self):
        # viewer: extract cc from database and plot it
        # make 5 pdbs with higher score available to scipion
        conn = sqlite3.connect(os.path.abspath(self._getExtraPath(DATAFILE)))
        c = conn.cursor()
        sqlCommand = """SELECT filename, model_to_map_fit
                        FROM   %s
                        ORDER BY model_to_map_fit DESC
                        LIMIT 5""" % TABLE
        c.execute(sqlCommand)
        rows = c.fetchall()

        argsOutput = {}
        for counter, row in enumerate(rows):
            atomStructFn = row[0][:-4] + "_real_space_refined.cif"
            atomStruct = AtomStruct()
            atomStruct.setFileName(atomStructFn)
            argsOutput["outpubAtomStruct_%d" % counter ] = atomStruct
            #self._defineSourceRelation(self.inputStructure.get(), atomStruct)
        c.close()
        conn.close()
        self._defineOutputs(**argsOutput)

    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        # Check that the input volume exist
        if self._getInputVolume() is None:
            errors.append("Error: You should provide a volume.\n")
        return errors

#    def _citations(self):
#        return ['Barad_2015']

    def _summary(self):
        summary=""
#        summary = PhenixProtRunRefinementBase._summary(self)
#        summary.append(
#            "https://www.phenix-online.org/documentation/reference/"
#            "real_space_refine.html")
        return summary

    def _writeArgsRSR(self, atomStruct, vol):
        args = "model_file=%s" % atomStruct
        args += " map_file=%s" % vol
        args += " resolution=%f" % self.resolution
        args += " secondary_structure.enabled=%s" % self.doSecondary
        args += " run="
        if self.minimizationGlobal == True:
            args += "minimization_global+"
        if self.rigidBody == True:
            args += "rigid_body+"
        if self.localGridSearch == True:
            args += "local_grid_search+"
        if self.morphing == True:
            args += "morphing+"
        if self.simulatedAnnealing == True:
            args += "simulated_annealing+"
        if self.adp == True:
            args += "adp+"
        args = args[:-1]
        # args += " run=minimization_global+local_grid_search+morphing+simulated_annealing"
        args += " macro_cycles=%d" % self.macroCycles
        args += " model_format=pdb+mmcif"
        # args += " write_pkl_stats=True"
        args += " %s " % self.extraParams.get()
        numberOfThreads = self.numberOfThreads.get()
        if numberOfThreads > 1:
            args += " nproc=%d" % numberOfThreads
        return args
