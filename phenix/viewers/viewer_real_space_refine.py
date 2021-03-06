# *
# * Authors:     Roberto Marabini
# *              Marta Martinez
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

from phenix.protocols.protocol_real_space_refine import PhenixProtRunRSRefine
from .viewer_validation_cryoem import PhenixProtRunValidationCryoEMViewer

class PhenixProtRunRSRefineViewer(PhenixProtRunValidationCryoEMViewer):
    """ Viewer for Phenix program real space refine
    """
    _label = 'Real Space Refine viewer'
    _targets = [PhenixProtRunRSRefine]

    def __init__(self,  **kwargs):
         PhenixProtRunValidationCryoEMViewer.__init__(self, **kwargs)
         REALSPACEFILE = self.protocol._getExtraPath(
             self.protocol.REALSPACEFILE)
         self.protocol._parseFile(REALSPACEFILE)

    def _defineParams(self, form):
        PhenixProtRunValidationCryoEMViewer._defineParams(self,form)






