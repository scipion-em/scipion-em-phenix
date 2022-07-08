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

from .test_protocol_emringer import (TestImportBase, TestImportData,
                                    TestEMRingerValidation2)
from .test_protocol_molprobity import (TestImportBase, TestImportData,
                                      TestMolprobityValidation2)
from .test_protocol_real_space_refine import (TestImportBase, TestImportData,
                                             TestPhenixRSRefine)
from .test_protocol_validation_cryoem import (TestImportBase, TestImportData,
                                             TestValCryoEM)
from .test_protocol_superpose_pdbs import (TestImportBase, TestProtSuperposePdbs,
                                          TestProtSuperposePdbs)
from .test_phenix import (TestVersion)
from .test_phenix_pdb_cif import (TestImportBase, TestImportData,
                                 TestPhenixPdbCif)
from .test_protocol_dock_in_map import TestProtDockInMap
from .test_protocol_search_fit import TestPhenixProtSearchFit
from .test_phenix_alphafold import TestAProtProcessDockBuildPredictedAlphaFold, \
    TestBProtProcessDockBuildPredictedAlphaFold, TestCProtProcessDockBuildPredictedAlphaFold
