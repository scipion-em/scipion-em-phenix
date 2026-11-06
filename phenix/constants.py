# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

PHENIX_HOME = 'PHENIX_HOME'
PHENIXVERSIONFILENAME = 'phenix_env.sh'
PHENIXVERSION = '2.1'
PHENIXVERSION210 = 2.1
PHENIXVERSIONdev_6085 = 'dev-6085'
# TO BE DELETED
# PHENIXVERSION = PHENIXVERSION210
# PHENIXVERSION18 = PHENIXVERSION
# PHENIXVERSION19 = PHENIXVERSION
# PHENIXVERSION20 = PHENIXVERSION
# PHENIXVERSION21 = PHENIXVERSION

# python used to run phenix scripts
# it is in the bin directory
PHENIX_PYTHON = 'phenix.python '  # keep the ending space

# phenix binaries are in several directories.
# Their location depends on the phenix version / installation layout, so these
# paths are configurable: set PHENIX_SCRIPT_PATH1 / PHENIX_SCRIPT_PATH2 in the
# [PLUGINS] section of scipion.conf to override the defaults below.
#
# NOTE: the paths differ between phenix major versions:
#   - phenix-1 (e.g. 1.21.2-5419, prebuilt "modules" layout):
#         PHENIX_SCRIPT_PATH1 = 'modules/cctbx_project/mmtbx/command_line'
#         PHENIX_SCRIPT_PATH2 = 'modules/phenix/phenix/command_line'
#   - phenix-2 (e.g. 2.2.1-6174, conda/python 3.11 site-packages layout):
#         PHENIX_SCRIPT_PATH1 = 'lib/python3.11/site-packages/mmtbx/command_line'
#         PHENIX_SCRIPT_PATH2 = 'lib/python3.11/site-packages/phenix/command_line'
#
# The defaults below correspond to the phenix-2 layout (current supported
# version). Users installing phenix-1 have to override them in scipion.conf.
#
# If the variables are not set, the values from scipion.conf (read by the
# scipion launcher into the environment) are used; otherwise these defaults
# apply.
#
#   PATH1 -> mmtbx scripts (molprobity, emringer, process_predicted_model, ...)
#   PATH2 -> phenix scripts (real_space_refine, superpose_pdbs, validation_cryoem, ...)
PHENIX_SCRIPT_PATH1 = os.environ.get('PHENIX_SCRIPT_PATH1', 'lib/python3.11/site-packages/mmtbx/command_line')
PHENIX_SCRIPT_PATH2 = os.environ.get('PHENIX_SCRIPT_PATH2', 'lib/python3.11/site-packages/phenix/command_line')
PHENIX_GETVERSION_PATH = 'bin'

# list of phenix scripts and corresponding binary directory
SUPERPOSE = 'superpose_pdbs.py'
REALSPACEREFINE = 'real_space_refine.py'
MOLPROBITY = 'molprobity.py'
# MOLPROBITY2 = 'molprobity.py'  # to be deleted
VALIDATION_CRYOEM = 'validation_cryoem.py'
EMRINGER = 'emringer.py'
GETVERSION = 'phenix.version'
DOCKINMAP = 'dock_in_map.py'
SYMMETRY = 'map_symmetry.py'
PROCESS = 'process_predicted_model.py'
DOCKPREDICTEDMODEL = 'dock_predicted_model.py'
REBUILDDOCKPREDICTEDMODEL = 'rebuild_predicted_model.py'
DOCKANDREBUILD = 'dock_and_rebuild.py'
PREDICTANDBUILD = 'predict_and_build.py'
mapBinarytoDirectory ={
    REALSPACEREFINE : PHENIX_SCRIPT_PATH2,
    SUPERPOSE : PHENIX_SCRIPT_PATH2,
    MOLPROBITY : PHENIX_SCRIPT_PATH1,
    # MOLPROBITY2 : PHENIX_SCRIPT_PATH1,
    VALIDATION_CRYOEM: PHENIX_SCRIPT_PATH2,
    EMRINGER : PHENIX_SCRIPT_PATH1,
    GETVERSION: PHENIX_GETVERSION_PATH,
    DOCKINMAP: PHENIX_SCRIPT_PATH2,
    SYMMETRY: PHENIX_SCRIPT_PATH2,
    PROCESS: PHENIX_SCRIPT_PATH1,
    DOCKPREDICTEDMODEL: PHENIX_SCRIPT_PATH2,
    REBUILDDOCKPREDICTEDMODEL: PHENIX_SCRIPT_PATH2,
    DOCKANDREBUILD: PHENIX_SCRIPT_PATH2,
    PREDICTANDBUILD: PHENIX_SCRIPT_PATH2
}
DISPLAY='display'

VERSION = "2.1"
