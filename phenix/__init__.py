# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
# *              Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
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
import pwem
import subprocess

import pyworkflow.utils as pwutils
from phenix.constants import *
from pwem.constants import MAXIT

_logo = "phenix.png"
_references = ['Adams_2010']
__version__ = "3.2.1"


class Plugin(pwem.Plugin):
    _homeVar = PHENIX_HOME
    _pathVars = [PHENIX_HOME]
    _supportedVersions = PHENIXVERSION

    @classmethod
    def _defineVariables(cls):
        if Plugin.getPhenixVersion() == PHENIXVERSION:
            cls._defineEmVar(PHENIX_HOME, 'phenix-1.18.2')
        else:
            cls._defineEmVar(PHENIX_HOME, ('phenix-' + (Plugin.getPhenixVersion())))

    @classmethod
    def getEnviron(cls, first=True):
        # environ = pwutils.Environ(os.environ)
        environ = pwutils.Environ()
        pos = pwutils.Environ.BEGIN if first else pwutils.Environ.END
        # add to variable
        if "DISPLAY" in os.environ:
            display = os.environ['DISPLAY']
        else:
            display = ''
        if "HOME" in os.environ:
            home = os.environ['HOME']
        else:
            home = ''

        environ.update({
            'DISPLAY': display,
            'HOME': home,
            'LIBTBX_BUILD': os.path.join(cls.getHome(), 'build'),
            'LIBTBX_OPATH': os.environ['PATH'],
            'PATH': os.path.join(Plugin.getHome(), 'build', 'bin') +
                    ':/usr/bin:'
                    '/bin'
        }, position=pos)
        return environ

    @classmethod
    def runPhenixProgram(cls, program, args=None, extraEnvDict=None, cwd=None):
        """ Internal shortcut function to launch a Phenix program. """
        env = cls.getEnviron()
        if extraEnvDict is not None:
            env.update(extraEnvDict)
        program = PHENIX_PYTHON + program
        pwutils.runJob(None, program, args, env=env, cwd=cwd)

    @classmethod
    def getProgram(cls, progName):
        """ Return the program binary that will be used. """
        return os.path.join(Plugin.getHome(),
                            mapBinarytoDirectory[progName],
                            os.path.basename(progName))

    @classmethod
    def getProgram(cls, progName):
        """ Return the program binary that will be used. """
        return os.path.join(Plugin.getHome(),
                            mapBinarytoDirectory[progName],
                            os.path.basename(progName))

    @classmethod
    def getPhenixVersion(cls):
        """ Return the program binary that will be used. """
        pid = subprocess.Popen(cls.getProgram("phenix.version"),
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        stdout, stderr = pid.communicate()
        searchKey = "Version: "
        stdout = stdout.decode("utf-8")
        place = stdout.find(searchKey)
        return stdout[place + len(searchKey):place + len(searchKey)+4]

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(PHENIXVERSION)

    @classmethod
    def defineBinaries(cls, env):
        pwem.Plugin.defineBinaries(env)
        env.getTarget(MAXIT).setDefault(True)
