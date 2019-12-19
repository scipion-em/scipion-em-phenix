# **************************************************************************
# *
# * Authors:     Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
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
# Horrible hack to release this plugin before scipion next version.
# TODO: remove when possible
from pyworkflow import LAST_VERSION, VERSION_2_0
if LAST_VERSION == VERSION_2_0 :
    from pyworkflow.utils import importFromPlugin
    retry = importFromPlugin('chimera.atom_struct', 'retry')
    fromCIFTommCIF = importFromPlugin('chimera.atom_struct', 'fromCIFTommCIF')
    fromCIFToPDB = importFromPlugin('chimera.atom_struct', 'fromCIFToPDB')
    fromPDBToCIF = importFromPlugin('chimera.atom_struct', 'fromPDBToCIF')
    AtomicStructHandler = importFromPlugin('chimera.atom_struct', 'AtomicStructHandler')
else:
    from pyworkflow.em.convert.atom_struct import retry
    from pyworkflow.em.convert.atom_struct import fromCIFTommCIF, fromCIFToPDB, fromPDBToCIF
    from pyworkflow.em.convert.atom_struct import AtomicStructHandler

from protocol_emringer import PhenixProtRunEMRinger
from protocol_molprobity import PhenixProtRunMolprobity
from protocol_real_space_refine import PhenixProtRunRSRefine
from protocol_refinement_base import PhenixProtRunRefinementBase
from protocol_superpose_pdbs import PhenixProtRunSuperposePDBs
from protocol_validation_cryoem import PhenixProtRunValidationCryoEM