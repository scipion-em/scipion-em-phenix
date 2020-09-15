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
mapBinarytoDirectory ={
    REALSPACEREFINE : PHENIX_SCRIPT_PATH2,
    SUPERPOSE : PHENIX_SCRIPT_PATH2,
    MOLPROBITY : PHENIX_SCRIPT_PATH1,
    MOLPROBITY2 : PHENIX_SCRIPT_PATH1,
    VALIDATION_CRYOEM: PHENIX_SCRIPT_PATH2,
    EMRINGER : PHENIX_SCRIPT_PATH1,
    GETVERSION: PHENIX_GETVERSION_PATH,
    DOCKINMAP: PHENIX_SCRIPT_PATH2,
    SYMMETRY: PHENIX_SCRIPT_PATH2
}
DISPLAY='display'

PHENIX_TO_SCIPION = {}
PHENIX_CYCLIC = 0  # SYM_CYCLIC = 0
PHENIX_DIHEDRAL_X = 1  # SYM_DIHEDRAL_X = SYM_DIHEDRAL = 1
PHENIX_TETRAHEDRAL = 2  # SYM_TETRAHEDRAL* = 2
PHENIX_OCTAHEDRAL = 3  # SYM_OCTAHEDRAL = 3
PHENIX_I = 4  # SYM_I* = 4

# symmetry dictionary
import pwem.constants as sciSym

PHENIX_TO_SCIPION[PHENIX_CYCLIC] = sciSym.SYM_CYCLIC
PHENIX_TO_SCIPION[PHENIX_DIHEDRAL_X] = sciSym.SYM_DIHEDRAL_X
PHENIX_TO_SCIPION[PHENIX_TETRAHEDRAL] = [sciSym.SYM_TETRAHEDRAL, sciSym.SYM_TETRAHEDRAL_Z3]
PHENIX_TO_SCIPION[PHENIX_OCTAHEDRAL] = sciSym.SYM_OCTAHEDRAL
PHENIX_TO_SCIPION[PHENIX_I] = [sciSym.SYM_I222, sciSym.SYM_I222r, sciSym.SYM_In25,
                               sciSym.SYM_In25r, sciSym.SYM_I2n3, sciSym.SYM_I2n3r,
                               sciSym.SYM_I2n5, sciSym.SYM_I2n5r]

PHENIX_SYM_NAME = dict()
PHENIX_SYM_NAME[PHENIX_CYCLIC] = 'Cn'
PHENIX_SYM_NAME[PHENIX_DIHEDRAL_X] = 'Dn'
PHENIX_SYM_NAME[PHENIX_TETRAHEDRAL] = 'T'
PHENIX_SYM_NAME[PHENIX_OCTAHEDRAL] = 'O'
PHENIX_SYM_NAME[PHENIX_I] = 'I'