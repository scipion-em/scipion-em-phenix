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

PHENIX_HOME = 'PHENIX_HOME'
PHENIXVERSIONFILENAME = './phenix_env.sh'
PHENIXVERSION = '1.13' # plugin version
PHENIXVERSION18 = '1.18' # september 2020
PHENIXVERSION19 = '1.19' # june 2021
PHENIXVERSION20 = '1.20' # march 2022

#python used to run phenix scripts
PHENIX_PYTHON = 'phenix.python '  # keep the ending space

#phenix binaries are in several directories
PHENIX_SCRIPT_PATH1 = 'modules/cctbx_project/mmtbx/command_line'
PHENIX_SCRIPT_PATH2 = 'modules/phenix/phenix/command_line'
PHENIX_GETVERSION_PATH = 'build/bin'

# list of phenix scripts and corresponding binary directory
SUPERPOSE = 'superpose_pdbs.py'
REALSPACEREFINE = 'real_space_refine.py'
MOLPROBITY = 'molprobity.py'
MOLPROBITY2 = 'molprobity.py'
VALIDATION_CRYOEM = 'validation_cryoem.py'
EMRINGER = 'emringer.py'
GETVERSION = 'phenix.version'
DOCKINMAP = 'dock_in_map.py'
SYMMETRY = 'map_symmetry.py'
PROCESS = 'process_predicted_model.py'
DOCKPREDICTEDMODEL = 'dock_predicted_model.py'
REBUILDDOCKPREDICTEDMODEL = 'rebuild_predicted_model.py'
DOCKANDREBUILD = 'dock_and_rebuild.py'
mapBinarytoDirectory ={
    REALSPACEREFINE : PHENIX_SCRIPT_PATH2,
    SUPERPOSE : PHENIX_SCRIPT_PATH2,
    MOLPROBITY : PHENIX_SCRIPT_PATH1,
    MOLPROBITY2 : PHENIX_SCRIPT_PATH1,
    VALIDATION_CRYOEM: PHENIX_SCRIPT_PATH2,
    EMRINGER : PHENIX_SCRIPT_PATH1,
    GETVERSION: PHENIX_GETVERSION_PATH,
    DOCKINMAP: PHENIX_SCRIPT_PATH2,
    SYMMETRY: PHENIX_SCRIPT_PATH2,
    PROCESS: PHENIX_SCRIPT_PATH2,
    DOCKPREDICTEDMODEL: PHENIX_SCRIPT_PATH2,
    REBUILDDOCKPREDICTEDMODEL: PHENIX_SCRIPT_PATH2,
    DOCKANDREBUILD: PHENIX_SCRIPT_PATH2
}
DISPLAY='display'