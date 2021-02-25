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

from phenix.protocols.protocol_validation_cryoem import PhenixProtRunValidationCryoEM
from .viewer_refinement_base import PhenixProtRefinementBaseViewer
from pyworkflow.protocol.params import LabelParam, EnumParam
from phenix import Plugin, PHENIXVERSION, PHENIXVERSION18
import matplotlib.pyplot as plt
import matplotlib.font_manager
from phenix import Plugin
import collections
import json


class PhenixProtRunValidationCryoEMViewer(PhenixProtRefinementBaseViewer):
    """ Viewer for Phenix program validation cryoEM
    """
    _label = 'Validation cryoEM viewer'
    _targets = [PhenixProtRunValidationCryoEM]

    VALIDATIONTMPFILE = 'tmpValidationFile.txt'
    CLASHESFILE = 'tmpClashesFile.py'
    CCCHAINFILE = "ccPerChain.py"
    CCCHAINFILE2 = "ccPerChain.txt"
    CCRESIDUESFILE = "ccPerResidue.py"
    CCRESIDUESFILE2 = "ccPerResidue.txt"
    FSCMODELMAPFILE = "FscModelMap.txt"
    FSCHALFMAPFILE = "FscHalfMap.txt"

    def __init__(self,  **kwargs):
         PhenixProtRefinementBaseViewer.__init__(self, **kwargs)
         if Plugin.getPhenixVersion() != PHENIXVERSION:
            # VALIDATIONCRYOEMFILE = self.protocol._getExtraPath(
            #     self.protocol.VALIDATIONCRYOEMFILE)

            self.VALIDATIONCRYOEMPKLFILENAME = self.protocol._getExtraPath(
                self.protocol.VALIDATIONCRYOEMPKLFILE)
            self._writePickleData2()
            self.dictOverall2 = json.loads(self.dictOverall2,
                                  object_pairs_hook=collections.OrderedDict)

    def _defineParams(self, form):
        if Plugin.getPhenixVersion() != PHENIXVERSION:
            if (self.protocol.inputVolume.get() \
                or self.protocol.inputStructure.get().getVolume()) \
                    is not None:
                if self.protocol.inputVolume.get() is not None:
                    self.vol = self.protocol.inputVolume.get()
                else:
                    self.vol = self.protocol.inputStructure.get().getVolume()
            PhenixProtRefinementBaseViewer._defineParams(self, form)
            form.addSection(label="Summary")
            group = form.addGroup('Model')
            group.addParam('showModelSummary', LabelParam,
                           label="Model Summary Table",
                           help="Model: Atomic structure.")
            group = form.addGroup('Data')
            group.addParam('showDataSummary', LabelParam,
                           label="Data Summary Table",
                           help="Data: 3D map obtained by cryo-EM.")
            group = form.addGroup('Model vs. Data')
            group.addParam('showModelVsDataSummary', LabelParam,
                           label="Model vs Data Summary Table",
                           help="Real space correlation values between 3D map and model.")

            form.addSection(label="MolProbity")
            if (self.dictOverall2['Clashes_n_outliers'] > 0):
                group = form.addGroup('Clashes: All-atom contact analysis')
                group.addParam('showClashes2', LabelParam, important = True,
                               label = "Bad contacts from PROBE (list)",
                               help = "This list summarizes all severe clashes "
                               "(more than 0.4 Angstroms non-H-bond "
                               "overlap) found by PROBE; you can view "
                               "these graphically in Coot. If no "
                               "hydrogens were present, REDUCE was "
                               "used to add them prior to running "
                               "PROBE.")
                group.addParam('exportFiles1', LabelParam,
                               label = 'Save list as text')
            else:
                group.addParam('showMesgNoClashes2', LabelParam,
                               label = "No bad contacts (> 0.4 Angstroms "
                                       "overlap) found.")
            group = form.addGroup('CaBLAM')
            if (self.dictOverall2['CaBLAM_outliers_n'] > 0):
                group.addParam('showCaBLAM', LabelParam, important = True,
                               label = 'CaBLAM evaluation',
                               help = 'C-based low-resolution annotation method')
            else:
                group.addParam('showMesgNoCaBLAM', LabelParam,
                               label="No CaBLAM outliers found.")
            group = form.addGroup('C-beta deviation analysis')
            if (self.dictOverall2['Cbeta_Outliers_n'] > 0):
                group.addParam('showCbetaOutliersTable2', LabelParam, important=True,
                               label='C-beta Outliers',
                               help="C-beta position outliers (position "
                                         "deviates from ideal by more than "
                                         "0.25A).\n\nIdeal CB position is "
                                         "determined from the average of the "
                                         "ideal C-N-CA-CB and N-C-CA-CB dihedrals."
                                         " This measure is more sensitive than"
                                         " individual measures to both sidechain "
                                         "and mainchain misfittings.\n")
            else:
                group.addParam('showMesgNoCbetaOutliers2', LabelParam,
                               label="No C-beta position outliers detected",
                               help="C-beta position outliers (position "
                                         "deviates from ideal by more than "
                                         "0.25A).\n\nIdeal CB position is "
                                         "determined from the average of the "
                                         "ideal C-N-CA-CB and N-C-CA-CB dihedrals."
                                         " This measure is more sensitive than"
                                         " individual measures to both sidechain "
                                         "and mainchain misfittings.\n")
            group = form.addGroup('Cis and twisted peptides')
            if (self.dictOverall2['Omega_outliers'] > 0):
                group.addParam('showCisAndTwistedPeptides2', LabelParam,
                                    important=True,
                                    label="Cis and Twisted peptides:",
                                    help="Cis conformations are observed in "
                                         "about 5% of Prolines.\n\nCis "
                                         "conformations are observed in about "
                                         "0.03% of general residues.\n\nTwisted "
                                         "peptides are almost certainly "
                                         "modeling errors.\n")
            else:
                group.addParam('showMesgNoNonTransPeptides2', LabelParam,
                                label="No non-trans peptides detected",
                                help="Cis conformations are observed in "
                                        "about 5% of Prolines.\n\nCis "
                                        "conformations are observed in about "
                                        "0.03% of general residues.\n\nTwisted "
                                        "peptides are almost certainly "
                                        "modeling errors.\n")
            group = form.addGroup('Rotamers')
            if (self.dictOverall2['Rota_Outliers_n'] > 0):
                group.addParam('showRotaOutliersTable2', LabelParam,
                                    important=True,
                                    label="Rotamer outlier list",
                                    help="Although a residue may lie in the "
                                     "favored regions of the Chi1-Chi2 plot, "
                                     "outliers are flagged based on the "
                                     "distribution of all non-branched Chi "
                                     "angles in a residue.\nZero outliers is "
                                     "not the goal. Rotamer outliers can be  "
                                     "justified by sufficiently strong " \
                                     "electron density, van der Waals "
                                     "packing, and/or hydrogen bonds.\n")
            else:
                group.addParam('showMesgNoRotaOutliers2', LabelParam,
                                    label="No Rotamer outliers detected")
            group.addParam('displayPlotRota', LabelParam, important=True,
                           label="Chi1-Chi2 graphs",
                           help="Visualization of Rotamer outliers in Chi1-Chi2 plots.")
            group = form.addGroup('Rhamachandran')
            if (self.dictOverall2['Rhama_Outliers_n'] > 0):
                group.addParam('showRhamaOutliersTable2', LabelParam,
                                    important=True,
                                    label="Rhamachandran outlier list",
                                    help="Ramachandran outliers are those aminoacids"
                                         " with non-favourable dihedral angles. "
                                         "Most of the time, Ramachandran outliers are "
                                         "a consequence of mistakes during the data "
                                         "processing.")
            else:
                group.addParam('showMesgNoRamaOutliers2', LabelParam,
                                    label="No Rhamachandran outliers detected")
            group.addParam('displayPlotRhama', LabelParam, important=True,
                           label="Rhamachandran graphs",
                           help="Visualization of outliers in Rhamachandran graphs.")
            group = form.addGroup('Geometry Restraints')
            group.addParam('showHelp', LabelParam, label='',
                           help="This section reports statistics for geometry restraints "
                                "used in refinement. As a general rule, a fully refined "
                                "structure should not have any outliers unless these are "
                                "exceptionally clear in the electron density (usually at "
                                "very high resolution). Be sure to also check the MolProbity "
                                "validation results for this structure, as they are more "
                                "sensitive to the geometric properties of proteins and nucleic "
                                "acids (especially in the case of dihedral angles).")
            group.addParam('showBLrestraints2', LabelParam,
                           important=True,
                           label="Bond length",
                           help="Check here the number of outlier pairs of atoms "
                                "according to the bond length restraints "
                                "between pairs of linked atoms.\nWarning!!!: "
                                "Refined structures should not have any outliers"
                                " except those are obvious in high "
                                "resolution electron density maps.\n")
            if self.dictOverall2['Length_outliers'] > 0:
                group.addParam('showBLoutliers2', LabelParam,
                               label="List of bond length outliers",
                               help="List of outlier pairs of atoms (sorted by deviation) "
                                    "according to the bond length restraints.\n")
            group.addParam('showBArestraints2', LabelParam, important=True,
                           label="Bond angle",
                           help="Check here the number of outlier triplets of atoms "
                                "according to the bond angle restraints.\n"
                                "Warning!!!: Refined structures should not "
                                "have any outliers except those are obvious in "
                                "high resolution electron density maps.")
            if self.dictOverall2['Angles_outliers'] > 0:
                group.addParam('showBAoutliers2', LabelParam,
                               label="List of bond angle outliers",
                               help="List of outlier triplets of atoms (sorted by "
                                    "deviation) according to the bond angle "
                                    "restraints")
            group.addParam('showDArestraints2', LabelParam, important=True,
                           label="Dihedral angle",
                           help="Check here the number of outlier tetrads of atoms "
                                "according "
                                "to the side chain dihedral torsion (chi) angle "
                                "restraints.\n"
                                "Warning!!!: Refined structures should not "
                                "have any outliers except those are obvious in "
                                "high resolution electron density maps.")
            if self.dictOverall2['Dihedral_outliers'] > 0:
                group.addParam('showDAoutliers2', LabelParam,
                               label="List of dihedral angle outliers",
                               help="List of outlier tetrads of atoms ("
                                    "sorted by deviation) "
                                    "according to the dihedral "
                                    "angle restraints")
            group.addParam('showCHILrestraints2', LabelParam, important=True,
                           label="Chirality",
                           help="Check here the number of outlier tetrads of atoms "
                                "according to the volume chirality "
                                "restraints.\n"
                                "Warning!!!: Refined structures should not "
                                "have any outliers except those are obvious "
                                "in high resolution electron density maps.")
            group.addParam('showPLANARrestraints2', LabelParam, important=True,
                           label="Planarity",
                           help="Check here the number of outliers of planar "
                                "groups, such as aromatic rings,"
                                "according to the planar "
                                "restraints.\n"
                                "Warning!!!: Refined structures should not "
                                "have any outliers except those are obvious "
                                "in high resolution electron density maps.")
            group.addParam('showPARALrestraints', LabelParam, important=True,
                           label="Parallelity",
                           help="\n"
                                "Warning!!!: Refined structures should not "
                                "have any outliers except those are obvious "
                                "in high resolution electron density maps.")
            group.addParam('showNONBondDistancerestraints', LabelParam,
                           label="Non-bonded distance", important=True,
                           help="\n"
                                "Warning!!!: Refined structures should not "
                                "have any outliers except those are obvious "
                                "in high resolution electron density maps.")
            group = form.addGroup('Display of rotamer and Rhamachandran outliers and clashes')
            group.addParam('showCootOutliers', LabelParam,
                           important=True,
                           label="Open in Coot",
                           help="Interactive visualization of outliers and clashes"
                                " with Coot:\n\nRamachandran outliers\n"
                                "Rotamer outliers\nC-beta outliers\n"
                                "Severe clashes ")
            form.addSection(label="Model vs. Data")
            form.addParam('showCorCoefTable', LabelParam, important=True,
                           label="Overall correlation coefficients",
                           help="Real-space correlation coefficients\n\nFor a detailed "
                                "definition of global correlation metrics CC (mask), "
                                "CC (box), CC (volume), and CC (peaks) see Afonine, P. V., "
                                "Klaholz, B. P., Moriarty, N. W., Poon, B. K., Sobolev, "
                                "O. V., Terwilliger, T. C., Adams, P. D. "
                                "& Urzhumtsev, A. (2018).\n"
                                "https://doi.org/10.1107/S20597983180093.\n\n"
                                "CC (mask): Model map vs. experimental map correlation "
                                "coefficient calculated considering map values inside "
                                "a mask calculated around the macromolecule"
                                ".\n\nCC (box): Model map vs. experimental map correlation "
                                "coefficient calculated considering all grid points of the "
                                "box.\n\nCC (volume) and CC (peaks) compare only "
                                "map regions with the highest density values "
                                "and regions below a certain contouring "
                                "threshold level are ignored.\nCC (volume): "
                                "The map region considered is defined by "
                                "the N highest points inside the molecular "
                                "mask.\n\nCC (peaks): In this case, calculations "
                                "consider the union of regions defined by "
                                "the N highest peaks in the model-calculated "
                                "map and the N highest peaks in the "
                                "experimental map.\n\nLocal real-space correlation coefficients "
                                "CC (main chain) and CC (side chain) involve main skeleton chain "
                                "and lateral chains, respectively.\n\n")
            group = form.addGroup("Correlation graphs")
            if len(self.dictOverall2['Chain_list']) > 0:
                self.Chains_list = []
                for item in self.dictOverall2['Chain_list']:
                    self.Chains_list.append(item[0])
                group.addParam('displayPlotChains', LabelParam, important=True,
                               label="Plot CC vs. Chain ID")
                group.addParam('exportFiles2', LabelParam,
                               label='Save chain CC as text')
                group.addParam('selectChain', EnumParam, important=True,
                               choices=self.Chains_list, default=0,
                               label="Chain",
                               help="Select a chain of the macromolecule")
                group.addParam('displayPlotResidues', LabelParam, important=True,
                               label="Plot CC vs. Residue number of the selected Chain",
                               help="Dashed line indicates the average value of CC for "
                                    "the whole chain.")
                group.addParam('exportFiles3', LabelParam,
                               label='Save residue CC of the selected Chain as text')
            form.addSection(label="Data")
            group = form.addGroup('Summary')
            group.addParam("showBox", LabelParam, important=True,
                           label="Box info (unit cell)\n",
                           help="Cell dimensions of the map (pixels).")
            group.addParam("showBoxDetails", LabelParam,
                           label="Unit cell: " + self.dictOverall2['Unit cell'] +
                                 "\n" + "Space group: " +
                                 self.dictOverall2['Space group'] + "\n")
            group.addParam('showMapResolution', LabelParam, important=True,
                           label="Map Resolution Estimates (Angstroms)",
                           help="Resolution estimates of the map considering "
                                "both map experimental data and model-derived "
                                "information.\n\nd99: Resolution cutoff beyond "
                                "which Fourier map coefficients are negligibly "
                                "small. Calculated from the full map or from each"
                                " one of half maps [d99 (half map 1), d99 (half "
                                "map 2)].\nOverall Biso: Overall isotropic B-value."
                                "\nd_model: Resolution cutoff at which the model "
                                "map is the most similar to the target (experimental)"
                                " map. Requires map and model. For d_model to be "
                                "meaningful, model is expected to fit the map as good"
                                " as possible. d_model (B factors = 0) tries to avoid"
                                "the blurring of the map.\nFSC (model): d_FSC_model; "
                                "Resolution cutoff up to which the model and map "
                                "Fourier coefficients are similar at FSC values of 0, "
                                "0.143, 0.5.\nFSC (half map 1, 2) = 0.143 (d_fsc): d_FSC; "
                                "Highest resolution at which the experimental data are "
                                "confident. Obtained from FSC curve calculated using"
                                " two half-maps and taken at FSC=0.143. The two half"
                                " maps are required to compute this value.\nMask "
                                "smoothing radius (Angstroms): Radius of the default "
                                "soft mask used since sharp edges resulting from "
                                "applying a binary map may introduce Fourier artifacts.")
            group.addParam('showMapStatistics', LabelParam, important=True,
                           label="Map Statistics",
                           help="Origin: Coordinates of map origin.\nAll: Map box size "
                                "in pixels.\nMin, Max, Mean: Statistics of the map "
                                "electron density values.")
            if self.vol.getHalfMaps():
                group.addParam('showHalfMapCC', LabelParam, important=True,
                               label="Half-map CC: " +
                                     ("%.4f" % self.dictOverall2["HalfMapCC"]),
                               help="Correlation coefficient between both half maps.")
            group.addParam('showMapHistogram', LabelParam,
                           label="Histogram of Map Values",
                           help="Plot that shows the number of map grid points regarding"
                                " the electron density values. Histograms are "
                                "usually skewed to the right.")
            group.addParam('showMapHistValues', LabelParam,
                           label="Map Values (slots) vs. Map (count)",
                           help="Table of intervals of map electron density (slots of map"
                                " values) and number of grid points involved.")
            if self.vol.getHalfMaps():
                group = form.addGroup('FSC (Half-maps)')
                group.addParam("showPlotFSC1", LabelParam, important=True,
                               label="Plot FSC vs. resolution (Angstroms)\n",
                               help="FSC curve calculated using two half maps regarding "
                                    "the spatial frequency (1/Angstroms) and resolution "
                                    "(Angstroms).")
                group.addParam('exportFiles4', LabelParam,
                            label='Save FSC plot data as text')
            group = form.addGroup('FSC (Model-map)')
            group.addParam("showPlotFSC2", LabelParam, important=True,
                           label="Plot FSC vs. resolution (Angstroms)\n",
                           help="FSC curve calculated using the full map and the "
                                "model-derived map regarding the spatial frequency "
                                "(1/Angstroms) and resolution (Angstroms).")
            group.addParam('exportFiles5', LabelParam,
                           label='Save FSC plot data as text')
        elif Plugin.getPhenixVersion() == PHENIXVERSION:
            PhenixProtRefinementBaseViewer._defineParams(self, form)

    if Plugin.getPhenixVersion() != PHENIXVERSION:
        def _getVisualizeDict(self):
            return{
                   'displayMapModel': self._displayMapModel,
                   'showModelSummary': self._showModelSummary,
                   'showDataSummary': self._showDataSummary,
                   'showModelVsDataSummary': self._showModelVsDataSummary,
                   'showClashes2': self._showClashes2,
                   'exportFiles1': self._exportFiles1,
                   'showCaBLAM': self._showCaBLAM,
                   'showCbetaOutliersTable2': self._showCbetaOutliersTable2,
                   'showCisAndTwistedPeptides2': self._showCisAndTwistedPeptides2,
                   'showRotaOutliersTable2': self._showRotaOutliersTable2,
                   'displayPlotRota': self._displayPlotRota,
                   'showCootOutliers': self._showCootOutliers,
                   'showRhamaOutliersTable2': self._showRhamaOutliersTable2,
                   'displayPlotRhama': self._displayPlotRhama,
                   'showBLrestraints2':self._showBLrestraints2,
                   'showBLoutliers2': self._showBLoutliers,
                   'showBArestraints2': self._showBArestraints2,
                   'showBAoutliers2': self._showBAoutliers,
                   'showDArestraints2': self._showDArestraints2,
                   'showDAoutliers2': self._showDAoutliers,
                   'showCHILrestraints2': self._showCHILrestraints2,
                   'showPLANARrestraints2': self._showPLANARrestraints2,
                   'showPARALrestraints': self._showPARALrestraints,
                   'showNONBondDistancerestraints': self._showNONBondDistancerestraints,
                   'showCorCoefTable': self._showCorCoefTable,
                   'displayPlotChains': self._displayPlotChains,
                   'exportFiles2': self._exportFiles2,
                   'displayPlotResidues': self._displayPlotResidues,
                   'exportFiles3': self._exportFiles3,
                   'showMapResolution': self._showMapResolution,
                   'showMapStatistics': self._showMapStatistics,
                   'showMapHistogram': self._showMapHistogram,
                   'showMapHistValues': self._showMapHistValues,
                   'showPlotFSC1': self._showPlotFSC1,
                   'exportFiles4': self._exportFiles4,
                   'showPlotFSC2': self._showPlotFSC2,
                   'exportFiles5': self._exportFiles5
                  }
    if Plugin.getPhenixVersion() == PHENIXVERSION:
        def _getVisualizeDict(self):
            return{
                   'displayMapModel': self._displayMapModel,
                   'showMolProbityResults': self._visualizeMolProbityResults,
                   'showCootOutliers': self._showCootOutliers,
                   'showMissingAtoms': self._showMissingAtoms,
                   'showBLrestraints': self._showBLrestraints,
                   'showBLoutliers': self._showBLoutliers,
                   'showBArestraints': self._showBArestraints,
                   'showBAoutliers': self._showBAoutliers,
                   'showDArestraints': self._showDArestraints,
                   'showDAoutliers': self._showDAoutliers,
                   'showCHILrestraints': self._showCHILrestraints,
                   'showCHILoutliers': self._showCHILoutliers,
                   'showPLANARrestraints': self._showPLANARrestraints,
                   'showPLANARoutliers': self._showPLANARoutliers,
                   'showPlotType': self._showPlotType,
                   'showRamaOutliersTable': self._showRamaOutliersTable,
                   'showRotaOutliersTable': self._showRotaOutliersTable,
                   'showCbetaOutliersTable': self._showCbetaOutliersTable,
                   'showBackAsnGlnHisSidechains': self._showBackAsnGlnHisSidechains,
                   'showCisAndTwistedPeptides': self._showCisAndTwistedPeptides,
                   'showMultiCriterionPlot': self._showMultiCriterionPlot,
                   'showOverallRSCResults': self._showOverallRSCResults,
                   'showClashes': self._showClashes,
                   'showCCTable': self._showCCTable,
                   'displayFSCplot': self._displayFSCplot,
                   'showOccupancies' : self._showOccupancies,
                   'showIsotropicB': self._showIsotropicB,
                   'showSuspiciousBfactors': self. _showSuspiciousBfactors
                   }

    def _showModelSummary(self, e = None):
        headerList = ['Item', 'Value']
        dataList1_1 = ['Composition (#)', '     Chains', '     Atoms', '     Residues',
                       '     Water', '     Ligands', 'Bonds (RMSD)',
                       '     Length (Angstroms) (# > 4sigma)',
                       '     Angles (degrees) (# > 4sigma)', 'MolProbity score',
                       'Clash score', 'Ramachandran plot (%)', '     Outliers',
                       '     Allowed', '     Favored']
        if Plugin.getPhenixVersion() >= PHENIXVERSION18:
            dataList1_2 = ['Rama-Z (Ramachandran plot Z-score, RMSD)',
                           '     whole (N = ' + str(self.dictOverall2['Rama_Z_whole_n']) + ')',
                           '     helix (N = ' + str(self.dictOverall2['Rama_Z_helix_n']) + ')',
                           '     sheet (N = ' + str(self.dictOverall2['Rama_Z_sheet_n']) + ')',
                           '     loop (N = ' + str(self.dictOverall2['Rama_Z_loop_n']) + ')']
        dataList1_3 = ['Rotamer outliers (%)',
                       'Cbeta outliers (%)', 'Peptide plane (%)', '     Cis proline/general',
                       '     Twisted proline/general', 'CaBLAM outliers (%)', 'ADP (B-factors)',
                       '     Iso/Aniso (#)', '     min/max/mean', '           Protein',
                       '           Nucleotide', '           Ligand', '           Water',
                       'Occupancy', '     Mean', '     occ = 1 (%)', '     0 < occ < 1 (%)',
                       '     occ > 1 (%)']
        if Plugin.getPhenixVersion() >= PHENIXVERSION18:
            dataList1 = dataList1_1 + dataList1_2 + dataList1_3
        else:
            dataList1 = dataList1_1 + dataList1_3

        dataList2_1 = ["", self.dictOverall2['Chains'],
                     str(self.dictOverall2["Atoms"])+ " " + "(Hydrogens: " +
                     str(self.dictOverall2["Hydrogens"]) + ")",
                     "Protein: " + str(self.dictOverall2['Protein_residues']) + " " +
                     "Nucleotide: " + str(self.dictOverall2['Nucleotide_residues']),
                     self.dictOverall2['Water'], self.dictOverall2['Ligands'], '',
                     self.dictOverall2['Length'] + " (" +
                     str(self.dictOverall2['Length_outliers']) + ")",
                     self.dictOverall2['Angles'] + " (" +
                     str(self.dictOverall2['Angles_outliers']) + ")",
                     self.dictOverall2['MolProbity_score'],
                     self.dictOverall2['Clash_score'],
                     "", self.dictOverall2['Rhama_Outliers'],
                     self.dictOverall2['Rhama_Allowed'],
                     self.dictOverall2['Rhama_Favored']]
        if Plugin.getPhenixVersion() >= PHENIXVERSION18:
            dataList2_2 = ["", self.dictOverall2['Rama_Z_whole_value'] + " (" +
                         self.dictOverall2['Rama_Z_whole_std'] + ")",
                         self.dictOverall2['Rama_Z_helix_value'] + " (" +
                         self.dictOverall2['Rama_Z_helix_std'] + ")",
                         self.dictOverall2['Rama_Z_sheet_value'] + " (" +
                         self.dictOverall2['Rama_Z_sheet_std'] + ")",
                         self.dictOverall2['Rama_Z_loop_value'] + " (" +
                         self.dictOverall2['Rama_Z_loop_std'] + ")"]
        dataList2_3 = [ self.dictOverall2['Rota_Outliers'],
                     self.dictOverall2['Cbeta_Outliers'],
                     "", self.dictOverall2['Cis_proline'] + "/" +
                     self.dictOverall2['Cis_general'],
                     self.dictOverall2['Twisted_proline'] + "/" +
                     self.dictOverall2['Twisted_general'],
                     self.dictOverall2['CaBLAM_outliers'], "",
                     str(self.dictOverall2['n_iso']) + "/" +
                     str(self.dictOverall2['n_aniso']),
                     "", self.dictOverall2['protein_min'] + "/" +
                     self.dictOverall2['protein_max'] + "/" +
                     self.dictOverall2['protein_mean'],
                     self.dictOverall2['nucleotide_min'] + "/" +
                     self.dictOverall2['nucleotide_max'] + "/" +
                     self.dictOverall2['nucleotide_mean'],
                     self.dictOverall2['other_min'] + "/" +
                     self.dictOverall2['other_max']
                     + "/" + self.dictOverall2['other_mean'],
                     self.dictOverall2['water_min'] + "/" +
                     self.dictOverall2['water_max'] + "/" +
                     self.dictOverall2['water_mean'],
                     "",
                     self.dictOverall2['occupancy_mean'],
                     self.dictOverall2['occupancy_occ_1'],
                     self.dictOverall2['occupancy_0_occ_1'],
                     self.dictOverall2['occupancy_occ_higher_1']]

        if Plugin.getPhenixVersion() >= PHENIXVERSION18:
            dataList2 = dataList2_1 + dataList2_2 + dataList2_3
        else:
            dataList2 = dataList2_1 + dataList2_3

        dataList = []
        for a1, a2 in zip(dataList1, dataList2):
            dataList.append((a1, a2))

        mesg = "Model"
        title = "Summary Table"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showDataSummary(self, e = None):
        headerList = ['Item', '', '']
        dataList1 = ['Box', '     Length (Angstroms)', '     Angles (degrees)',
                     'Supplied Resolution (Angstroms)', 'Resolution Estimates (Angstroms)',
                     '     d FSC (half maps; 0.143)', '     d 99 (full/half1/half2)',
                     '     d model', '     d FSC model (0/0.143/0.5)', 'Map min/max/mean']
        dataList2 = ['', self.dictOverall2['Box_lengths'],
                     self.dictOverall2['Box_angles'],
                     ("%.1f" % self.dictOverall2['Supplied_Resolution']),
                     "Masked",
                     self.dictOverall2['dFSC_half_maps_0.143_masked'],
                     self.dictOverall2['d99_full_masked'] + "/" +
                     self.dictOverall2['d99_half1_masked'] + "/" +
                     self.dictOverall2['d99_half2_masked'],
                     self.dictOverall2['dmodel_masked'],
                     self.dictOverall2['dFSCmodel_0_masked'] + "/" +
                     self.dictOverall2['dFSCmodel_0.143_masked'] + "/" +
                     self.dictOverall2['dFSCmodel_0.5_masked'],
                     self.dictOverall2['Map_min'] + "/" +
                     self.dictOverall2['Map_max'] + "/" +
                     self.dictOverall2['Map_mean']]
        dataList3 = ['', '', '', '',
                     "Unmasked",
                     self.dictOverall2['dFSC_half_maps_0.143_unmasked'],
                     self.dictOverall2['d99_full_unmasked'] + "/" +
                     self.dictOverall2['d99_half1_unmasked'] + "/" +
                     self.dictOverall2['d99_half2_unmasked'],
                     self.dictOverall2['dmodel_unmasked'],
                     self.dictOverall2['dFSCmodel_0_unmasked'] + "/" +
                     self.dictOverall2['dFSCmodel_0.143_unmasked'] + "/" +
                     self.dictOverall2['dFSCmodel_0.5_unmasked'], '']

        dataList = []
        for a1, a2, a3 in zip(dataList1, dataList2, dataList3):
            dataList.append((a1, a2, a3))

        mesg = "Data"
        title = "Summary Table"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showModelVsDataSummary(self, e = None):
        headerList = ['Item', 'Value']
        dataList1 = ['CC (mask)', 'CC (box)', 'CC (peaks)', 'CC (volume)',
                     'Mean CC for ligands']
        dataList2 = [("%.2f" % self.dictOverall2['CC_mask']),
                     ("%.2f" % self.dictOverall2['CC_box']),
                     ("%.2f" % self.dictOverall2['CC_peaks']),
                     ("%.2f" % self.dictOverall2['CC_volume']),
                     self.dictOverall2['Ligand_CC']]

        dataList = []
        for a1, a2 in zip(dataList1, dataList2):
            dataList.append((a1, a2))

        mesg = "Model vs. Data"
        title = "Summary Table"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showClashes2(self, e=None):
        headerList = self.dictOverall2['Clashes_header']
        dataList = self.dictOverall2['Clashes_table']
        mesg = "Bad contacts from PROBE: %d overlapping atom pairs" \
               % len(dataList)
        title = "All atom-contact analyis"
        self._showOutliers(headerList, dataList, mesg, title)

    def _exportFiles1(self, e=None):
        CLASHESFILENAME = self.protocol._getExtraPath(self.CLASHESFILE)

        def onSelect(obj):
            dirName = obj.getPath()
            command = """import pickle
def pickleData(file):
    with open(file,"r") as f:
        return pickle.load(f)

# process file {VALIDATIONCRYOEMPKLFILENAME}"
data = pickleData('{VALIDATIONCRYOEMPKLFILENAME}')
data.model.geometry.clash.clashes.save_table_data("{dirName}/clashes.txt")
""".format(VALIDATIONCRYOEMPKLFILENAME=self.VALIDATIONCRYOEMPKLFILENAME,
                   dirName=dirName)
            with open(CLASHESFILENAME, "w") as f:
                f.write(command)
            # execute file with phenix.python
            Plugin.runPhenixProgram("", CLASHESFILENAME)

        self._openBrowser(onSelect)

    def _showCaBLAM(self, e=None):
        headerList = self.dictOverall2['CaBLAM_header']
        dataList = self.dictOverall2['CaBLAM_table']
        dataListNew = []
        for i in dataList:
            i1 = []
            for j in i[:9]:
                if type(j) != float:
                    i1.append(j)
                elif type(j) == float:
                    i1.append("%6.5f" % j)
            dataListNew.append(i1)

        mesg = "Outliers (%%): %s  Disfavored (%%): %.2f  Calpha outliers (%%): %.2f" \
               % (self.dictOverall2['CaBLAM_outliers'],
                  self.dictOverall2['CaBLAM_disfavored'],
                  self.dictOverall2['CaBLAM_Calpha_outliers'])
        title = "CaBLAM"
        self._showOutliers(headerList, dataListNew, mesg, title)

    def _showCbetaOutliersTable2(self, e=None):
        headerList = self.dictOverall2['Cbeta_header']
        dataList = self.dictOverall2['Cbeta_table']
        dataListNew = []
        for i in dataList:
            i1 = []
            for j in i:
                if type(j) != float:
                    i1.append(j)
                elif type(j) == float:
                    i1.append("%6.3f" % j)
            dataListNew.append(i1)
        mesg = "C-beta outliers"
        title = "C-beta deviation analysis"
        self._showOutliers(headerList, dataListNew, mesg, title)

    def _showCisAndTwistedPeptides2(self, e=None):
        headerList = self.dictOverall2['Omega_headers']
        dataList = self.dictOverall2['Omega_table']
        dataListNew = []
        for i in dataList:
            i1 = []
            for j in i:
                if type(j) != float:
                    i1.append(j)
                elif type(j) == float:
                    i1.append("%6.2f" % j)
            dataListNew.append(i1)
        mesg = "Cis conformations are observed in about 5% of Prolines.\n" \
               "Cis conformations are observed in about 0.03% of general residues.\n" \
               "Twisted peptides are almost certainly modeling errors."
        title = "Cis and twisted peptides"
        self._showOutliers(headerList, dataListNew, mesg, title)

    def _showRotaOutliersTable2(self, e=None):
        headerList = self.dictOverall2['Rota_header']
        dataList = self.dictOverall2['Rota_table']
        dataListNew = []
        for i in dataList:
            i1 = []
            for j in i:
                if type(j) != float:
                    i1.append(j)
                elif j == i[2] and type(j) == float:
                    i1.append("%6.2f" % j)
                elif j != i[2] and type(j) == float:
                    i1.append("%6.1f" % j)
            dataListNew.append(i1)
        mesg = "Outlier list"
        title = "Rotamer analysis"
        self._showOutliers(headerList, dataListNew, mesg, title)

    def _displayPlotRota(self, e=None):
        self._displayPlotBase(self.ROTATMPFILE, 1)

    def _displayPlotBase(self, file, listNumber):
        FILENAME = self.protocol._getExtraPath(file)
        self._writeCommand2(listNumber)
        with open(FILENAME, "w") as f:
            f.write(self.command)
        # execute file with phenix.python
        Plugin.runPhenixProgram("", FILENAME)

    def _writeCommand2(self, listNumber):
        self.command ="""import pickle     
def pickleData(file):
    with open(file,"r") as f:
        return pickle.load(f)
        
# process file {VALIDATIONCRYOEMPKLFILENAME}"
data = pickleData('{VALIDATIONCRYOEMPKLFILENAME}')
""".format(VALIDATIONCRYOEMPKLFILENAME=self.VALIDATIONCRYOEMPKLFILENAME)
        if (listNumber == 0 or listNumber == 1):
            if (listNumber == 0):
                self.command += """# Rhamachandran plot
if data.model.geometry.ramachandran.ramalyze is not None:
"""
            elif (listNumber == 1):
                self.command += """# Rotamer plot
if data.model.geometry.rotamer.rotalyze is not None:
"""
            self.command +="""
    try :
        import wxtbx.app
    except ImportError, e :
        raise Sorry("wxPython not available.")
    app = wxtbx.app.CCTBXApp(0)
"""
            if (listNumber == 0):
                self.command += """
    data.model.geometry.ramachandran.ramalyze.display_wx_plots() 
"""
            elif (listNumber == 1):
                self.command += """
    data.model.geometry.rotamer.rotalyze.display_wx_plots() 
"""
            self.command +="""
    app.MainLoop()
"""
        elif (listNumber == 2):
            self.command += """# Chain plot
if data.model_vs_data.cc is not None:
    import matplotlib.pyplot as plt
    chain_names = []
    chain_cc = []
    for item in data.model_vs_data.cc.cc_per_chain:
        if item.chain_id not in chain_names:
            chain_names.append(item.chain_id)
            chain_cc.append(item.cc)
    assert(len(chain_names) == len(chain_cc))
    fig = plt.figure(figsize=(10,4))
    p = fig.add_subplot(111)
    fig.subplots_adjust(left = 0.10, right = 0.975, bottom = 0.25, top = 0.80)
    p.bar(chain_names, chain_cc, width = 0.9)
    p.set_xlabel('Chain ID')
    p.set_ylabel('CC')
    p.set_ylim(0, 1)
    p.set_title('Plot of correlation coefficients regarding chain IDs', fontsize = 14)
    x = range(len(chain_cc))
    p.set_xticks(x)
    p.set_xticklabels(chain_names)
    fig.canvas.draw()
    plt.show()
"""
        elif (listNumber == 3):
            chain = self.Chains_list[self.selectChain.get()]
            self.command += """# Chain plot
if data.model_vs_data.cc is not None:
    import matplotlib.pyplot as plt
    resseq_list = []
    residue_cc = []
    for item in data.model_vs_data.cc.cc_per_residue:
        if item.chain_id == "%s":
            resseq_list.append(item.resseq)
            residue_cc.append(item.cc)
    chain_mean = sum(residue_cc) / len(residue_cc) 
    assert(len(resseq_list) == len(residue_cc))
    fig = plt.figure(figsize=(15,4))
    p = fig.add_subplot(111)
    fig.subplots_adjust(left = 0.10, right = 0.975, bottom = 0.25, top = 0.80)
    p.plot(resseq_list, residue_cc, 'b-', linewidth=1)
    p.set_xlabel('Residue number')
    p.set_ylabel('CC')
    p.set_ylim(-0.1, 1)
    p.set_title('Plot of correlation coefficients regarding residue numbers of chain "%s"', fontsize = 14)
    intv = int(round((len(resseq_list))/5))
    x_axis = [0, intv, 2 * intv, 3 * intv, 4 * intv, 5 * intv]
    p.set_xticks(x_axis)
    p.set_xticklabels(x_axis)
    # draw line for mean
    p.axhline(y = chain_mean, color = 'k', linestyle = 'dashed')
    fig.canvas.draw()
    plt.show()
"""% (chain, chain)

    def _showRhamaOutliersTable2(self, e=None):
        headerList = self.dictOverall2['Rhama_header']
        dataList = self.dictOverall2['Rhama_table']
        dataListNew = []
        for i in dataList:
            i1 = []
            for j in i:
                if type(j) != float:
                    i1.append(j)
                elif j == i[3] and type(j) == float:
                    i1.append("%6.2f" % j)
                elif j != i[3] and type(j) == float:
                    i1.append("%6.1f" % j)
            dataListNew.append(i1)
        mesg = "Outlier list"
        title = "Rhamachandran analysis"
        self._showOutliers(headerList, dataListNew, mesg, title)

    def _displayPlotRhama(self, e=None):
        self._displayPlotBase(self.RAMATMPFILE, 0)

    def _showBLrestraints2(self, e=None):
        headerList = ['measure', 'value']
        dictBL = collections.OrderedDict()
        dictBL['Number of restraints:'] = self.dictOverall2['BLRestraints']
        dictBL['RMS (deviation):'] = self.dictOverall2['BLMean']
        dictBL['Max. deviation:'] = self.dictOverall2['BLMax']
        dictBL['Min. deviation:'] = self.dictOverall2['BLMin']
        dictBL['Number of outliers > 4sigma:'] = self.dictOverall2['Length_outliers']
        val = 0.3
        mesg = "Bond Length Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Geometry Restraints"
        self._showResults(headerList, dictBL, val, mesg, title)

    def _showBArestraints2(self, e=None):
        headerList = ['measure', 'value']
        dictBA = collections.OrderedDict()
        dictBA['Number of restraints:'] = self.dictOverall2['BARestraints']
        dictBA['RMS (deviation):'] = self.dictOverall2['BAMean']
        dictBA['Max. deviation:'] = self.dictOverall2['BAMax']
        dictBA['Min. deviation:'] = self.dictOverall2['BAMin']
        dictBA['Number of outliers > 4sigma:'] = self.dictOverall2['Angles_outliers']
        val = 0.3
        mesg = "Bond Angle Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Geometry Restraints"
        self._showResults(headerList, dictBA, val, mesg, title)

    def _showDArestraints2(self, e=None):
        headerList = ['measure', 'value']
        dictDA = collections.OrderedDict()
        dictDA['Number of restraints:'] = self.dictOverall2['DARestraints']
        dictDA['RMS (deviation):'] = self.dictOverall2['DAMean']
        dictDA['Max. deviation:'] = self.dictOverall2['DAMax']
        dictDA['Min. deviation:'] = self.dictOverall2['DAMin']
        dictDA['Number of outliers > 4sigma:'] = self.dictOverall2['Dihedral_outliers']
        val = 0.3
        mesg = "Dihedral Angle Restraints\n(Deviations from ideal values)"
        title = "MolProbity: Geometry Restraints"
        self._showResults(headerList, dictDA, val, mesg, title)

    def _showCHILrestraints2(self, e=None):
        headerList = ['measure', 'value']
        dictChil = collections.OrderedDict()
        dictChil['Number of restraints:'] = self.dictOverall2['ChiralityRestraints']
        dictChil['RMS (deviation):'] = self.dictOverall2['ChiralityMean']
        dictChil['Max. deviation:'] = self.dictOverall2['ChiralityMax']
        dictChil['Min. deviation:'] = self.dictOverall2['ChiralityMin']
        val = 0.3
        mesg = "Chirality restraints\n(Deviations from ideal values)"
        title = "MolProbity: Geometry Restraints"
        self._showResults(headerList, dictChil, val, mesg, title)

    def _showPLANARrestraints2(self, e=None):
        headerList = ['measure', 'value']
        dictPLANAR = collections.OrderedDict()
        dictPLANAR['Number of restraints:'] = self.dictOverall2['PlanarityRestraints']
        dictPLANAR['RMS (deviation):'] = self.dictOverall2['PlanarityMean']
        dictPLANAR['Max. deviation:'] = self.dictOverall2['PlanarityMax']
        dictPLANAR['Min. deviation:'] = self.dictOverall2['PlanarityMin']
        val = 0.3
        mesg = "Planarity restraints\n(Deviations from ideal values)"
        title = "MolProbity: Geometry Restraints"
        self._showResults(headerList, dictPLANAR, val, mesg, title)

    def _showPARALrestraints(self, e=None):
        headerList = ['measure', 'value']
        dictPARAL = collections.OrderedDict()
        dictPARAL['Number of restraints:'] = self.dictOverall2['ParallelityRestraints']
        dictPARAL['RMS (deviation):'] = self.dictOverall2['ParallelityMean']
        dictPARAL['Max. deviation:'] = self.dictOverall2['ParallelityMax']
        dictPARAL['Min. deviation:'] = self.dictOverall2['ParallelityMin']
        val = 0.3
        mesg = "Parallelity restraints\n(Deviations from ideal values)"
        title = "MolProbity: Geometry Restraints"
        self._showResults(headerList, dictPARAL, val, mesg, title)

    def _showNONBondDistancerestraints(self, e=None):
        headerList = ['measure', 'value']
        dictNONBondDistance = collections.OrderedDict()
        dictNONBondDistance['Number of restraints:'] = self.dictOverall2['NonbondedRestraints']
        dictNONBondDistance['RMS (deviation):'] = self.dictOverall2['NonbondedMean']
        dictNONBondDistance['Max. deviation:'] = self.dictOverall2['NonbondedMax']
        dictNONBondDistance['Min. deviation:'] = self.dictOverall2['NonbondedMin']
        val = 0.3
        mesg = "Non-bonded distance restraints\n(Deviations from ideal values)"
        title = "MolProbity: Geometry Restraints"
        self._showResults(headerList, dictNONBondDistance, val, mesg, title)

    def _showCorCoefTable(self, e=None):
        headerList = ['statistic', 'value']
        dictCC = collections.OrderedDict()
        dictCC["CC (mask): "] = self.dictOverall2['CC_mask']
        dictCC["CC (box): "] = self.dictOverall2['CC_box']
        dictCC["CC (volume): "] = self.dictOverall2['CC_volume']
        dictCC["CC (peaks): "] = self.dictOverall2['CC_peaks']
        dictCC["CC (main chain): "] = self.dictOverall2['CC_main_chain']
        dictCC["CC (side chain): "] = self.dictOverall2['CC_side_chain']
        val = 0.2
        mesg = "Overall correlation coefficients"
        title = "Model vs. Data"
        self._showResults(headerList, dictCC, val, mesg, title)

    def _displayPlotChains(self, e=None):
        self._displayPlotBase(self.CCCHAINFILE, 2)

    def _displayPlotResidues(self, e=None):
        self._displayPlotBase(self.CCRESIDUESFILE, 3)

    def _exportFiles2(self, e=None):
        Chain_list = self.dictOverall2['Chain_list']
        CCCHAINFILE2NAME = self.CCCHAINFILE2
        def onSelect(obj):
            dirName = obj.getPath()
            Path = ("{dirName}/{fileName}").format(dirName=dirName, fileName=CCCHAINFILE2NAME)
            with open(Path, "w") as f:
                f.write('CC per chain\n')
                for item in Chain_list:
                    f.write('   ' + item[0] + '   ' + str(item[1]) + '\n')

        self._openBrowser(onSelect)

    def _exportFiles3(self, e=None):
        Residue_list = self.dictOverall2['Residue_list']
        CCRESIDUESFILE2NAME = self.CCRESIDUESFILE2
        chain = self.Chains_list[self.selectChain.get()]
        def onSelect(obj):
            dirName = obj.getPath()
            Path = ("{dirName}/{fileName}").format(dirName=dirName, fileName=CCRESIDUESFILE2NAME)
            with open(Path, "w") as f:
                f.write('CC per residue\n')
                for item in Residue_list:
                    if item[1] == chain:
                        f.write('   ' + item[0] + '   ' + item[1] + '   ' + item[2] +
                                '   ' + item[3] + '   ' + str(item[4]) + '\n')

        self._openBrowser(onSelect)

    def _openBrowser(self, onSelect):
        from pyworkflow.gui.browser import FileBrowserWindow
        browser = FileBrowserWindow("Select a folder to save file",
                                     master=self.getWindow(),
                                     path=".",
                                     onSelect=onSelect,
                                     onlyFolders=True)
        browser.show()

    def _showMapResolution(self, e=None):
        headerList = ['', 'Masked', 'Unmasked']
        dataList1_1 = ["Using map alone (d99): ", "Overall Biso: ", "d_model: ",
                     "d_model (B factors = 0): ", "FSC (model) = 0: ",
                     "FSC (model) = 0.143: ", "FSC (model) = 0.5: "]
        if self.vol.getHalfMaps():
            dataList1_2 = ["d99 (half map 1): ", "d99 (half map 2): ",
                           "FSC (half map 1,2) = 0.143 (d_fsc): "]
        dataList1_3 = ["Mask smoothing radius (Angstroms): "]
        dataList2_1 = [self.dictOverall2['*d99_full_masked'],
                       self.dictOverall2["overall_b_iso_masked"],
                       self.dictOverall2['*dmodel_masked'],
                       self.dictOverall2['d_model_b0_masked'],
                       self.dictOverall2['*dFSCmodel_0_masked'],
                       self.dictOverall2['*dFSCmodel_0.143_masked'],
                       self.dictOverall2['*dFSCmodel_0.5_masked']]
        if self.vol.getHalfMaps():
            dataList2_2 = [self.dictOverall2['d99_half1_masked2'],
                           self.dictOverall2['d99_half2_masked2'],
                           self.dictOverall2['*dFSC_half_maps_0.143_masked']]
        dataList2_3 = [self.dictOverall2['mask_smoothing_radius']]
        dataList3_1 = [self.dictOverall2['*d99_full_unmasked'],
                       self.dictOverall2["overall_b_iso_unmasked"],
                       self.dictOverall2['*dmodel_unmasked'],
                       self.dictOverall2['d_model_b0_unmasked'],
                       self.dictOverall2['*dFSCmodel_0_unmasked'],
                       self.dictOverall2['*dFSCmodel_0.143_unmasked'],
                       self.dictOverall2['*dFSCmodel_0.5_unmasked']]
        if self.vol.getHalfMaps():
            dataList3_2 = [self.dictOverall2['d99_half1_unmasked2'],
                           self.dictOverall2['d99_half2_unmasked2'],
                           self.dictOverall2['*dFSC_half_maps_0.143_unmasked']]
        dataList3_3 = ['']
        if self.vol.getHalfMaps():
            dataList1 = dataList1_1 + dataList1_2 + dataList1_3
            dataList2 = dataList2_1 + dataList2_2 + dataList2_3
            dataList3 = dataList3_1 + dataList3_2 + dataList3_3
        else:
            dataList1 = dataList1_1  + dataList1_3
            dataList2 = dataList2_1  + dataList2_3
            dataList3 = dataList3_1  + dataList3_3
        dataList = []
        for a1, a2, a3 in zip(dataList1, dataList2, dataList3):
            dataList.append((a1, a2, a3))
        mesg = "Data"
        title = "Map Resolution (Angstroms) Estimates"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showMapStatistics(self, e=None):
        headerList = ['', '']
        dataList1 = ["Origin: ", "Last: ", "Focus: ",
                     "All: ", "Min, Max, Mean: "]
        dataList2 = [str(tuple(self.dictOverall2['Map_origin'])),
                     str(tuple(self.dictOverall2['Map_last'])),
                     str(tuple(self.dictOverall2['Map_focus'])),
                     str(tuple(self.dictOverall2['Map_all'])),
                     str("(" + self.dictOverall2['Map_min_max_mean'] + ")")]
        dataList = []
        for a1, a2 in zip(dataList1, dataList2):
            dataList.append((a1, a2))
        mesg = "Map Statistics"
        title = "Data"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showMapHistogram(self, e=None):
        a = [pow(10, i) for i in range(10)]
        fig = plt.figure(figsize=(15, 4))
        fig.subplots_adjust(left=0.10, right=0.975, bottom=0.25, top=0.80)

        if self.vol.getHalfMaps():
            splots = []
            assert (len(self.dictOverall2['Map_value']) ==
                    len(self.dictOverall2['HalfMap1_count']))
            assert (len(self.dictOverall2['Map_value']) ==
                    len(self.dictOverall2['HalfMap2_count']))

            p1 = fig.add_subplot(121)
            p1.plot(self.dictOverall2['Map_value'], self.dictOverall2['Map_count'],
                   'b-', linewidth=1, label="Map")
            splots.append(p1)

            p2 = fig.add_subplot(122)
            p2.plot(self.dictOverall2['Map_value'], self.dictOverall2['HalfMap1_count'],
                   'b-', linewidth=1, label="Half-map 1")
            p2.plot(self.dictOverall2['Map_value'], self.dictOverall2['HalfMap2_count'],
                   'r-', linewidth=1, label="Half-map 2")
            splots.append(p2)
            for p in splots:
                # Place the legend to the upper right of this subplot
                p.legend()
                p.set_xlabel('Map value')
                p.set_ylabel('Count')
                p.set_yscale('log')
            fig.suptitle('Histogram of Map Values', fontsize=14)

        else:
            assert (len(self.dictOverall2['Map_value']) == len(self.dictOverall2['Map_count']))
            p = fig.add_subplot(111)
            fig.subplots_adjust(left=0.10, right=0.975, bottom=0.25, top=0.80)
            p.plot(self.dictOverall2['Map_value'], self.dictOverall2['Map_count'],
                'b-', linewidth=1, label="Map")
            # Place the legend to the upper right of this subplot
            p.legend()
            p.set_xlabel('Map value')
            p.set_ylabel('Count')
            p.set_yscale('log')
            p.set_title('Histogram of Map Values', fontsize=14)
        fig.canvas.draw()
        plt.show()

    def _showMapHistValues(self, e=None):
        dataList = []
        Left = []
        Right = []
        Map_count = self.dictOverall2['Map_count']
        slots = self.dictOverall2['Map_slots'].split('\n')
        for item in range(len(slots)):
            slots_parts = slots[item].split()
            Left.append("%.4f" % round(float(slots_parts[0]), 4))
            Right.append("%.4f" % round(float(slots_parts[2].split(':')[0]), 4))

        if self.vol.getHalfMaps():
            headerList = ['Map_values', 'Map', 'Half-map 1', 'Half-map 2']
            HalfMap1_count = self.dictOverall2['HalfMap1_count']
            HalfMap2_count = self.dictOverall2['HalfMap2_count']
            for a1, a2, a3, a4, a5 in zip(
                    Left, Right, Map_count, HalfMap1_count, HalfMap2_count):
                dataList.append((str(a1) + "  -  " + str(a2), a3, a4, a5))
        else:
            headerList = ['Map_values', 'Map']
            for a1, a2, a3 in zip(Left, Right, Map_count):
                dataList.append((str(a1) + "  -  " + str(a2), a3))

        mesg = "Map Values (slots) vs. Map (count)"
        title = "Data"
        self._showOutliers(headerList, dataList, mesg, title)

    def _showPlotFSC1(self, e=None):
        fsc = self.dictOverall2['FSC_Masked']
        x_d_inv = self.dictOverall2['d_inv_Masked']
        fsc_u = self.dictOverall2['FSC_Unmasked']
        x_d_inv_u = self.dictOverall2['d_inv_Unmasked']
        d = 0.143
        self._show_plot(x_d_inv, fsc, x_d_inv_u, fsc_u, d,
                        label1="Masked", label2="Unmasked")

    def _showPlotFSC2(self, e=None):
        fsc = self.dictOverall2['FSC_Model_Map_Masked']
        x_d_inv = self.dictOverall2['d_inv_Model_Map_Masked']
        fsc_u = self.dictOverall2['FSC_Model_Map_Unmasked']
        x_d_inv_u = self.dictOverall2['d_inv_Model_Map_Unmasked']
        d = 0.5
        self._show_plot(x_d_inv, fsc, x_d_inv_u, fsc_u, d,
                        label1="Masked", label2="Unmasked")

    def _show_plot(self, x1, y1, x2, y2, d, label1, label2):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twiny()
        l1, = ax1.plot(x1, y1, 'r-', linewidth=1, label=label1)
        l2, = ax2.plot(x2, y2, 'b-', linewidth=1, label=label2)

        ax1.set_ylim(-0.2, 1.05)
        ax1.set_ylabel("FSC",
                       fontproperties=matplotlib.font_manager.FontProperties(
                           family=['Helvetica', 'sans-serif'],
                           weight='bold',
                           size=12))

        xticks = []
        max_x = max(max(x1), max(x2))
        for i in range(int((round(max_x, 1) + 0.05) / 0.05) + 1):
            xticks.append(round(i * 0.05, 2))
        ax1.set_xticks(xticks)
        ax2.set_xticks(xticks)

        ax1.set_xlabel('1/resolution (1/$\AA$)',
                       fontproperties=matplotlib.font_manager.FontProperties(
                           family=['Helvetica', 'sans-serif'],
                           weight='bold',
                           size=12))
        labels1 = []
        for i in range(len(xticks)):
            labels1.append('%.2f' % xticks[i])
        ax1.set_xticklabels(labels1,
                            fontproperties=matplotlib.font_manager.FontProperties(
                                family=['Courier', 'Monaco', 'monospace'],
                                weight='normal',
                                size=12))
        ax1.axhline(y=d, linewidth=1, color='k', linestyle='dashed')
        idx = float()
        if d == 0.5 or \
                (d == 0.143 and self.dictOverall2["HalfMapCC"] != 1.):
            for i in range(len(y1) - 1):
                a = y1[i] - d
                b = y1[i + 1] - d
                if ((a > 0 and b < 0) or (a == 0 and b < 0) or (a > 0 and b == 0)):
                    idx = (x1[i] + x1[i + 1]) / 2
                else:
                    pass
        if idx is not None:
            ax1.axvline(x=idx, color='k', linestyle='dashed')

        ax2.set_xlabel('resolution ($\AA$)',
                       fontproperties=matplotlib.font_manager.FontProperties(
                           family=['Helvetica', 'sans-serif'],
                           weight='bold',
                           size=12))

        labels2 = []
        for i in range(len(xticks)):
            if (xticks[i] > 0.0):
                x = 1.0 / xticks[i]
                labels2.append('%.1f' % x)
            else:
                labels2.append('\u221E')
        ax2.set_xticklabels(labels2,
                            fontproperties=matplotlib.font_manager.FontProperties(
                                family=['Courier', 'Monaco', 'monospace'],
                                weight='normal',
                                size=12))
        if d == 0.5 or \
                (d == 0.143 and self.dictOverall2["HalfMapCC"] != 1.):
            for i in range(len(y2) - 1):
                a = y2[i] - d
                b = y2[i + 1] - d
                if ((a > 0 and b < 0) or (a == 0 and b < 0) or (a > 0 and b == 0)):
                    idx = (x2[i] + x2[i + 1]) / 2
                else:
                    pass
        if idx is not None:
            ax2.axvline(x=idx, color='k', linestyle='dashed')
        fig.legend((l2, l1), (label2, label1), 'upper right')
        plt.tight_layout()
        plt.show()

    def _exportFiles4(self, e=None):
        fsc = self.dictOverall2['FSC_Masked']
        x_d_inv = self.dictOverall2['d_inv_Masked']
        fsc_u = self.dictOverall2['FSC_Unmasked']
        x_d_inv_u = self.dictOverall2['d_inv_Unmasked']
        FILENAME = self.FSCHALFMAPFILE
        self._export_file(fsc, x_d_inv, fsc_u, x_d_inv_u, FILENAME)

    def _exportFiles5(self, e=None):
        fsc = self.dictOverall2['FSC_Model_Map_Masked']
        x_d_inv = self.dictOverall2['d_inv_Model_Map_Masked']
        fsc_u = self.dictOverall2['FSC_Model_Map_Unmasked']
        x_d_inv_u = self.dictOverall2['d_inv_Model_Map_Unmasked']
        FILENAME = self.FSCMODELMAPFILE
        self._export_file(fsc, x_d_inv, fsc_u, x_d_inv_u, FILENAME)

    def _export_file(self, fsc, x_d_inv, fsc_u, x_d_inv_u, fileName):

        def onSelect(obj):
            dirName = obj.getPath()
            Path = ("{dirName}/{fileName}").format(dirName=dirName, fileName=fileName)
            with open(Path, "w") as f:
                f.write('Unmasked\n')
                f.write(' ' + 'd_inv' + '            ' + 'fsc' + '\n')
                for a1, a2 in zip(x_d_inv_u, fsc_u):
                    f.write((' ' + str("%.12f") + '   ' + str("%.12f") + '\n') % (a1, a2))
                f.write('\n\nMasked\n')
                f.write(' ' + 'd_inv' + '            ' + 'fsc' + '\n')
                for a1, a2 in zip(x_d_inv, fsc):
                    f.write((' ' + str("%.12f") + '   ' + str("%.12f") + '\n') % (a1, a2))

        self._openBrowser(onSelect)

    def _writePickleData2(self):
        VALIDATIONTMPFILENAME = self.protocol._getExtraPath(
            self.VALIDATIONTMPFILE)
        command = """
import pickle
import collections
import json

def pickleData(file):
    with open(file,"r") as f:
        return pickle.load(f)
        
# process file {VALIDATIONCRYOEMPKLFILENAME}"
data = pickleData('{VALIDATIONCRYOEMPKLFILENAME}')
dictOverall2 = collections.OrderedDict()

# composition
chain_names = []
for item in data.model_vs_data.cc.cc_per_chain:
    if item.chain_id not in chain_names:
        chain_names.append(item.chain_id)   
dictOverall2['Chains'] = len(chain_names)
dictOverall2['Atoms'] = data.model.composition.n_atoms
dictOverall2['Hydrogens'] = data.model.composition.n_hd
dictOverall2['Protein_residues'] = data.model.composition.n_protein
dictOverall2['Nucleotide_residues'] = data.model.composition.n_nucleotide
dictOverall2['Water'] = data.model.composition.n_water

Ligands = []
if data.model.composition.n_other != 0:
    for k, v  in data.model.composition.other_cnts.iteritems():
        Ligands.append(str(k) + ': ' + str(v))
        dictOverall2['Ligands'] = Ligands
else:
    dictOverall2['Ligands'] = '---'
        
# Bonds (RMSD)
dictOverall2['Length'] = "%.3f" % (round(data.model.geometry.bond.mean, 3))
dictOverall2['Length_outliers'] = len(data.model.geometry.bond.outliers)
dictOverall2['BLRestraints'] = data.model.geometry.bond.n
dictOverall2['BLMean'] = data.model.geometry.bond.mean
dictOverall2['BLMax'] = data.model.geometry.bond.max
dictOverall2['BLMin'] = data.model.geometry.bond.min
dictOverall2['Angles'] = "%.3f" % (round(data.model.geometry.angle.mean, 3))
dictOverall2['Angles_outliers'] = len(data.model.geometry.angle.outliers)
dictOverall2['BARestraints'] = data.model.geometry.angle.n
dictOverall2['BAMean'] = data.model.geometry.angle.mean
dictOverall2['BAMax'] = data.model.geometry.angle.max
dictOverall2['BAMin'] = data.model.geometry.angle.min
dictOverall2['Dihedral_outliers'] = len(data.model.geometry.dihedral.outliers)
dictOverall2['DARestraints'] = data.model.geometry.dihedral.n
dictOverall2['DAMean'] = data.model.geometry.dihedral.mean
dictOverall2['DAMax'] = data.model.geometry.dihedral.max
dictOverall2['DAMin'] = data.model.geometry.dihedral.min
dictOverall2['ChiralityRestraints'] = data.model.geometry.chirality.n
dictOverall2['ChiralityMean'] = data.model.geometry.chirality.mean
dictOverall2['ChiralityMax'] = data.model.geometry.chirality.max
dictOverall2['ChiralityMin'] = data.model.geometry.chirality.min
dictOverall2['PlanarityRestraints'] = data.model.geometry.planarity.n
dictOverall2['PlanarityMean'] = data.model.geometry.planarity.mean
dictOverall2['PlanarityMax'] = data.model.geometry.planarity.max
dictOverall2['PlanarityMin'] = data.model.geometry.planarity.min
dictOverall2['ParallelityRestraints'] = data.model.geometry.parallelity.n
dictOverall2['ParallelityMean'] = data.model.geometry.parallelity.mean
dictOverall2['ParallelityMax'] = data.model.geometry.parallelity.max
dictOverall2['ParallelityMin'] = data.model.geometry.parallelity.min
dictOverall2['NonbondedRestraints'] = data.model.geometry.nonbonded.n
dictOverall2['NonbondedMean'] = data.model.geometry.nonbonded.mean
dictOverall2['NonbondedMax'] = data.model.geometry.nonbonded.max
dictOverall2['NonbondedMin'] = data.model.geometry.nonbonded.min
        
# MolProbity score
dictOverall2['MolProbity_score'] = "%.2f" % (round(data.model.geometry.molprobity_score, 2))
        
# Clash score
dictOverall2['Clash_score'] = "%.2f" % (round(data.model.geometry.clash.score, 2))
dictOverall2['Clashes_n_outliers'] = data.model.geometry.clash.clashes.n_outliers
dictOverall2['Clashes_header'] = data.model.geometry.clash.clashes.gui_list_headers
dictOverall2['Clashes_table'] = data.model.geometry.clash.clashes.as_gui_table_data()
        
# Rhamachandran plot (%)
dictOverall2['Rhama_Outliers_n'] = data.model.geometry.ramachandran.outliers
dictOverall2['Rhama_Outliers'] = "%.2f" % (round(data.model.geometry.ramachandran.outliers, 2))
dictOverall2['Rhama_Allowed'] = "%.2f" % (round(data.model.geometry.ramachandran.allowed, 2))
dictOverall2['Rhama_Favored'] = "%.2f" % (round(data.model.geometry.ramachandran.favored, 2))
dictOverall2['Rhama_header'] = data.model.geometry.ramachandran.ramalyze.gui_list_headers
dictOverall2['Rhama_table'] = data.model.geometry.ramachandran.ramalyze.as_gui_table_data()
        
# Rotamer outliers (%)
dictOverall2['Rota_Outliers_n'] = data.model.geometry.rotamer.outliers
dictOverall2['Rota_Outliers'] = "%.2f" % (round(data.model.geometry.rotamer.outliers, 2))
dictOverall2['Rota_header'] = data.model.geometry.rotamer.rotalyze.gui_list_headers
dictOverall2['Rota_table'] = data.model.geometry.rotamer.rotalyze.as_gui_table_data()
        
# Cbeta outliers (%)
dictOverall2['Cbeta_Outliers_n'] = data.model.geometry.c_beta.cbetadev.n_outliers
dictOverall2['Cbeta_Outliers'] = "%.2f" % (round(data.model.geometry.c_beta.outliers, 2))
dictOverall2['Cbeta_header'] = data.model.geometry.c_beta.cbetadev.gui_list_headers
dictOverall2['Cbeta_table'] = data.model.geometry.c_beta.cbetadev.as_gui_table_data()
        
# Peptide plane (%)
dictOverall2['Cis_proline'] = "%.1f" % (round(data.model.geometry.omega.cis_proline, 1))
dictOverall2['Cis_general'] = "%.1f" % (round(data.model.geometry.omega.cis_general, 1))
dictOverall2['Twisted_proline'] = "%.1f" % (round(data.model.geometry.omega.twisted_proline, 1))
dictOverall2['Twisted_general'] = "%.1f" % (round(data.model.geometry.omega.twisted_general, 1))
dictOverall2['Omega_outliers'] = data.model.geometry.omega.omegalyze.n_outliers
dictOverall2['Omega_headers'] = data.model.geometry.omega.omegalyze.gui_list_headers
dictOverall2['Omega_table'] = data.model.geometry.omega.omegalyze.as_gui_table_data()
      
# CaBLAM outliers (%)
dictOverall2['CaBLAM_outliers_n'] = data.model.geometry.cablam.outliers
dictOverall2['CaBLAM_outliers'] = "%.2f" % (round(data.model.geometry.cablam.outliers, 2))
dictOverall2['CaBLAM_disfavored'] = round(data.model.geometry.cablam.disfavored, 2)
dictOverall2['CaBLAM_Calpha_outliers'] = round(data.model.geometry.cablam.ca_outliers, 2)
dictOverall2['CaBLAM_header'] = data.model.geometry.cablam.gui_table()['column_labels']
dictOverall2['CaBLAM_table'] = data.model.geometry.cablam.gui_table()['data']
        
# ADP (B-factors)
dictOverall2['n_iso'] = data.model.adp.overall.n_iso
dictOverall2['n_aniso'] = data.model.adp.overall.n_aniso
#min/max/mean
if data.model.adp.protein is not None:
    dictOverall2['protein_min'] = "%.2f" % (round(data.model.adp.protein.min, 2))
    dictOverall2['protein_max'] = "%.2f" % (round(data.model.adp.protein.max, 2))
    dictOverall2['protein_mean'] = "%.2f" % (round(data.model.adp.protein.mean, 2))
else:
    dictOverall2['protein_min'] = '---'
    dictOverall2['protein_max'] = '---'
    dictOverall2['protein_mean'] = '---'
if data.model.adp.nucleotide is not None:
    dictOverall2['nucleotide_min'] = "%.2f" % (round(data.model.adp.nucleotide.min, 2))
    dictOverall2['nucleotide_max'] = "%.2f" % (round(data.model.adp.nucleotide.max, 2))
    dictOverall2['nucleotide_mean'] = "%.2f" % (round(data.model.adp.nucleotide.mean, 2))
else:
    dictOverall2['nucleotide_min'] = '---'
    dictOverall2['nucleotide_max'] = '---'
    dictOverall2['nucleotide_mean'] = '---'
if data.model.adp.other is not None:
    dictOverall2['other_min'] = "%.2f" % (round(data.model.adp.other.min, 2))
    dictOverall2['other_max'] = "%.2f" % (round(data.model.adp.other.max, 2))
    dictOverall2['other_mean'] = "%.2f" % (round(data.model.adp.other.mean, 2))
else:
    dictOverall2['other_min'] = '---'
    dictOverall2['other_max'] = '---'
    dictOverall2['other_mean'] = '---'
if data.model.adp.water is not None:
   dictOverall2['water_min'] = "%.2f" % (round(data.model.adp.water.min, 2))
   dictOverall2['water_max'] = "%.2f" % (round(data.model.adp.water.max, 2))
   dictOverall2['water_mean'] = "%.2f" % (round(data.model.adp.water.mean, 2))
else:
   dictOverall2['water_min'] = '---'
   dictOverall2['water_max'] = '---'
   dictOverall2['water_mean'] = '---'     

# Occupancy
dictOverall2['occupancy_mean'] = "%.2f" % (round(data.model.occupancy.mean, 2))
dictOverall2['occupancy_occ_1'] = "%.2f" % (round(data.model.occupancy.equal_to_1_fraction, 2))
dictOverall2['occupancy_0_occ_1'] = "%.2f" % (round(data.model.occupancy.between_0_and_1_fraction, 2))
dictOverall2['occupancy_occ_higher_1'] = "%.2f" % (round(data.model.occupancy.greater_than_1_fraction, 2))
        
# Box
lengths = []
for item in data.data.crystal_symmetry.unit_cell().parameters()[:3]:
    lengths.append("%.2f" % item)
dictOverall2['Box_lengths'] = str(lengths[0]) + ", " + str(lengths[1]) + ", " + str(lengths[2])
angles = []
for item in data.data.crystal_symmetry.unit_cell().parameters()[3:]:
    angles.append("%.2f" % item)
dictOverall2['Box_angles'] = str(angles[0]) + ", " + str(angles[1]) + ", " + str(angles[2])
dictOverall2['Unit cell'] = str(data.data.crystal_symmetry.unit_cell().parameters())
dictOverall2['Space group'] = str(data.data.crystal_symmetry.space_group().info().symbol_and_number())
        
# Supplied Resolution
dictOverall2['Supplied_Resolution'] = round(data.model_vs_data.cc.atom_radius, 1)
        
# Resolution Estimates      Masked      Unmasked
if data.data.masked.d_fsc is not None:
    dictOverall2['dFSC_half_maps_0.143_masked'] = "%.1f" % (round(data.data.masked.d_fsc, 1))
    dictOverall2['*dFSC_half_maps_0.143_masked'] = "%.2f" % (round(data.data.masked.d_fsc, 2))
else:
    dictOverall2['dFSC_half_maps_0.143_masked'] = '---'
    dictOverall2['*dFSC_half_maps_0.143_masked'] = '---'
if data.data.unmasked.d_fsc is not None:
    dictOverall2['dFSC_half_maps_0.143_unmasked'] = "%.1f" % (round(data.data.unmasked.d_fsc, 1))
    dictOverall2['*dFSC_half_maps_0.143_unmasked'] = "%.2f" % (round(data.data.unmasked.d_fsc, 2))
else:
    dictOverall2['dFSC_half_maps_0.143_unmasked'] = '---'
    dictOverall2['*dFSC_half_maps_0.143_unmasked'] = '---'
if data.data.masked.d99 is not None:
    dictOverall2['d99_full_masked'] = "%.1f" % (round(data.data.masked.d99, 1))
    dictOverall2['*d99_full_masked'] = "%.2f" % (round(data.data.masked.d99, 2))
else:
    dictOverall2['d99_full_masked'] = '---'
    dictOverall2['*d99_full_masked'] = '---'
if data.data.unmasked.d99 is not None:
    dictOverall2['d99_full_unmasked'] = "%.1f" % (round(data.data.unmasked.d99, 1))
    dictOverall2['*d99_full_unmasked'] = "%.2f" % (round(data.data.unmasked.d99, 2))
else:
    dictOverall2['d99_full_unmasked'] = '---'
    dictOverall2['*d99_full_unmasked'] = '---'
if data.data.masked.d99_1 is not None:
    dictOverall2['d99_half1_masked'] = "%.1f" % (round(data.data.masked.d99_1, 1))
    dictOverall2['d99_half1_masked2'] = "%.2f" % (round(data.data.masked.d99_1, 2))
else:
    dictOverall2['d99_half1_masked'] = '---'
if data.data.unmasked.d99_1 is not None:
    dictOverall2['d99_half1_unmasked'] = "%.1f" % (round(data.data.unmasked.d99_1, 1))
    dictOverall2['d99_half1_unmasked2'] = "%.2f" % (round(data.data.unmasked.d99_1, 2))
else:
    dictOverall2['d99_half1_unmasked'] = '---'
if data.data.masked.d99_2 is not None:
    dictOverall2['d99_half2_masked'] = "%.1f" % (round(data.data.masked.d99_2, 1))
    dictOverall2['d99_half2_masked2'] = "%.2f" % (round(data.data.masked.d99_2, 2))
else:
    dictOverall2['d99_half2_masked'] = '---'
if data.data.unmasked.d99_2 is not None:
    dictOverall2['d99_half2_unmasked'] = "%.1f" % (round(data.data.unmasked.d99_2, 1))
    dictOverall2['d99_half2_unmasked2'] = "%.2f" % (round(data.data.unmasked.d99_2, 2))
else:
    dictOverall2['d99_half2_unmasked'] = '---'
if data.data.masked.d_model is not None:
    dictOverall2['dmodel_masked'] = "%.1f" % (round(data.data.masked.d_model, 1))
    dictOverall2['*dmodel_masked'] = "%.2f" % (round(data.data.masked.d_model, 2))
else:
    dictOverall2['dmodel_masked'] = '---'
    dictOverall2['*dmodel_masked'] = '---'
if data.data.unmasked.d_model is not None:
    dictOverall2['dmodel_unmasked'] = "%.1f" % (round(data.data.unmasked.d_model, 1))
    dictOverall2['*dmodel_unmasked'] = "%.2f" % (round(data.data.unmasked.d_model, 2))
else:
    dictOverall2['dmodel_masked'] = '---'
    dictOverall2['*dmodel_unmasked'] = '---'
if data.data.masked.d_model_b0 is not None:
   dictOverall2['d_model_b0_masked'] = "%.2f" % (round(data.data.masked.d_model_b0, 2))
else:
   dictOverall2['d_model_b0_masked'] = '---'
if data.data.unmasked.d_model_b0 is not None:
   dictOverall2['d_model_b0_unmasked'] = "%.2f" % (round(data.data.unmasked.d_model_b0, 2))
else:
   dictOverall2['d_model_b0_unmasked'] = '---'
if data.data.masked.d_fsc_model_0 is not None:
    dictOverall2['dFSCmodel_0_masked'] = "%.1f" % (round(data.data.masked.d_fsc_model_0, 1))
    dictOverall2['*dFSCmodel_0_masked'] = "%.2f" % (round(data.data.masked.d_fsc_model_0, 2))
else:
    dictOverall2['dFSCmodel_0_masked'] = '---'
    dictOverall2['*dFSCmodel_0_masked'] = '---'
if data.data.unmasked.d_fsc_model_0 is not None:
    dictOverall2['dFSCmodel_0_unmasked'] = "%.1f" % (round(data.data.unmasked.d_fsc_model_0, 1))
    dictOverall2['*dFSCmodel_0_unmasked'] = "%.2f" % (round(data.data.unmasked.d_fsc_model_0, 2))
else:
    dictOverall2['dFSCmodel_0_unmasked'] = '---'
    dictOverall2['*dFSCmodel_0_unmasked'] = '---'
if data.data.masked.d_fsc_model_0143 is not None:
    dictOverall2['dFSCmodel_0.143_masked'] = "%.1f" % (round(data.data.masked.d_fsc_model_0143, 1))
    dictOverall2['*dFSCmodel_0.143_masked'] = "%.2f" % (round(data.data.masked.d_fsc_model_0143, 2))
else:
    dictOverall2['dFSCmodel_0.143_masked'] = '---'
    dictOverall2['*dFSCmodel_0.143_masked'] = '---'
if data.data.unmasked.d_fsc_model_0143 is not None:
    dictOverall2['dFSCmodel_0.143_unmasked'] = "%.1f" % (round(data.data.unmasked.d_fsc_model_0143, 1))
    dictOverall2['*dFSCmodel_0.143_unmasked'] = "%.2f" % (round(data.data.unmasked.d_fsc_model_0143, 2))
else:
    dictOverall2['dFSCmodel_0.143_unmasked'] = '---'
    dictOverall2['*dFSCmodel_0.143_unmasked'] = '---'
if data.data.masked.d_fsc_model_05 is not None:
    dictOverall2['dFSCmodel_0.5_masked'] = "%.1f" % (round(data.data.masked.d_fsc_model_05, 1))
    dictOverall2['*dFSCmodel_0.5_masked'] = "%.2f" % (round(data.data.masked.d_fsc_model_05, 2))
else:
    dictOverall2['dFSCmodel_0.5_masked'] = '---'
    dictOverall2['*dFSCmodel_0.5_masked'] = '---'
if data.data.unmasked.d_fsc_model_05 is not None:
    dictOverall2['dFSCmodel_0.5_unmasked'] = "%.1f" % (round(data.data.unmasked.d_fsc_model_05, 1))
    dictOverall2['*dFSCmodel_0.5_unmasked'] = "%.2f" % (round(data.data.unmasked.d_fsc_model_05, 2))
else:
    dictOverall2['dFSCmodel_0.5_unmasked'] ='---'
    dictOverall2['*dFSCmodel_0.5_unmasked'] = '---'
if data.data.masked.b_iso_overall is not None:
    dictOverall2["overall_b_iso_masked"] = "%.2f" % (round(data.data.masked.b_iso_overall, 2))
else:
    dictOverall2["overall_b_iso_masked"] = '---'
if data.data.unmasked.b_iso_overall is not None:
    dictOverall2["overall_b_iso_unmasked"] = "%.2f" % (round(data.data.unmasked.b_iso_overall, 2))
else:
    dictOverall2["overall_b_iso_unmasked"] = '---'
if data.data.masked.radius_smooth is not None:
    dictOverall2['mask_smoothing_radius'] = "%.2f" % (round(data.data.masked.radius_smooth, 2))
else:
    dictOverall2['mask_smoothing_radius'] = '---'
          
# Map statistics
dictOverall2['Map_min'] = "%.2f" % (round(data.data.counts.min_max_mean[0], 2))
dictOverall2['Map_max'] = "%.2f" % (round(data.data.counts.min_max_mean[1], 2))
dictOverall2['Map_mean'] = "%.2f" % (round(data.data.counts.min_max_mean[2], 2))
dictOverall2['Map_origin'] = data.data.counts.origin
dictOverall2['Map_last'] = data.data.counts.last
dictOverall2['Map_focus'] = data.data.counts.focus
dictOverall2['Map_all'] = data.data.counts.all
min = round(data.data.counts.min_max_mean[0], 3)
max = round(data.data.counts.min_max_mean[1], 3)
mean = round(data.data.counts.min_max_mean[2], 3)
dictOverall2['Map_min_max_mean'] = ("%.3f, %.3f, %.3f") % (min, max, mean)
Map_value = []
for item in data.data.histograms.h_map.slot_centers():
    Map_value.append(item)
dictOverall2['Map_value'] = Map_value   
Map_count = []
for item in data.data.histograms.h_map.slots():
    Map_count.append(item)
dictOverall2['Map_count'] = Map_count
dictOverall2['Map_slots'] = data.data.histograms.h_map.as_str() 
if data.data.histograms.h_half_map_1 is not None:
    HalfMap1_count = []
    for item in data.data.histograms.h_half_map_1.slots():
        HalfMap1_count.append(item)
    dictOverall2['HalfMap1_count'] = HalfMap1_count
if data.data.histograms.h_half_map_2 is not None:    
    HalfMap2_count = []
    for item in data.data.histograms.h_half_map_2.slots():
        HalfMap2_count.append(item)
    dictOverall2['HalfMap2_count'] = HalfMap2_count
if data.data.histograms.half_map_histogram_cc is not None:
    dictOverall2["HalfMapCC"] = data.data.histograms.half_map_histogram_cc

# FSC (Half-maps)
if data.data.masked.fsc_curve is not None:
    fsc_masked = []
    d_inv_masked = []
    for item in data.data.masked.fsc_curve.fsc.fsc:
        fsc_masked.append(item)
    dictOverall2['FSC_Masked'] = fsc_masked
    for item in data.data.masked.fsc_curve.fsc.d_inv:
        d_inv_masked.append(item)
    dictOverall2['d_inv_Masked'] = d_inv_masked
if data.data.unmasked.fsc_curve is not None:
    fsc_unmasked = []
    d_inv_unmasked = []
    for item in data.data.unmasked.fsc_curve.fsc.fsc:
        fsc_unmasked.append(item)
    dictOverall2['FSC_Unmasked'] = fsc_unmasked
    for item in data.data.unmasked.fsc_curve.fsc.d_inv:
        d_inv_unmasked.append(item)
    dictOverall2['d_inv_Unmasked'] = d_inv_unmasked
    
# FSC (Model-map)
if data.data.masked.fsc_curve_model.fsc is not None:
    fsc_model_map_masked = []
    d_inv_model_map_masked = []
    for item in data.data.masked.fsc_curve_model.fsc:
        fsc_model_map_masked.append(item)
    dictOverall2['FSC_Model_Map_Masked'] = fsc_model_map_masked
    for item in data.data.masked.fsc_curve_model.d_inv:
        d_inv_model_map_masked.append(item)
    dictOverall2['d_inv_Model_Map_Masked'] = d_inv_model_map_masked
if data.data.unmasked.fsc_curve_model.fsc is not None:
    fsc_model_map_unmasked = []
    d_inv_model_map_unmasked = []
    for item in data.data.unmasked.fsc_curve_model.fsc:
        fsc_model_map_unmasked.append(item)
    dictOverall2['FSC_Model_Map_Unmasked'] = fsc_model_map_unmasked
    for item in data.data.unmasked.fsc_curve_model.d_inv:
        d_inv_model_map_unmasked.append(item)
    dictOverall2['d_inv_Model_Map_Unmasked'] = d_inv_model_map_unmasked
               
# Model vs. Data
dictOverall2['CC_mask'] = round(data.model_vs_data.cc.cc_mask, 2)
dictOverall2['CC_box'] = round(data.model_vs_data.cc.cc_box, 2)
dictOverall2['CC_peaks'] = round(data.model_vs_data.cc.cc_peaks, 2)
dictOverall2['CC_volume'] = round(data.model_vs_data.cc.cc_volume, 2)
dictOverall2['CC_main_chain'] = round(data.model_vs_data.cc.cc_main_chain.cc, 2)
try:  # skip if chain has no sidechains
    dictOverall2['CC_side_chain']= round(data.model_vs_data.cc.cc_side_chain.cc, 2)
except:
   pass

# Model structure
Chain_list = []
chain_names = []
for item in data.model_vs_data.cc.cc_per_chain:
    if item.chain_id not in chain_names:
        chain_names.append(item.chain_id)
        Chain_list.append((item.chain_id, item.cc))
        dictOverall2['Chain_list'] = Chain_list
Residue_list = []
for item in data.model_vs_data.cc.cc_per_residue:
    Residue_list.append((item.model_id, item.chain_id, item.resseq, item.resname, item.cc))
    dictOverall2['Residue_list'] = Residue_list
IUPAC_aminoacids_list = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS",
"LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
IUPAC_dna_list = ["DA", "DC", "DG", "DT"]
IUPAC_rna_list = ["A", "C", "G", "U"]
Ligand_CC = []
for item in data.model_vs_data.cc.cc_per_residue:
    if item.resname.strip() not in IUPAC_aminoacids_list:
        if item.resname.strip() not in IUPAC_dna_list:
            if item.resname.strip() not in IUPAC_rna_list:
                Ligand_CC.append(item.cc)
if len(Ligand_CC) > 0:
    dictOverall2['Ligand_CC'] = "%.2f" % round(sum(Ligand_CC)/len(Ligand_CC), 2)
else:
    dictOverall2['Ligand_CC'] = "---"
""".format(VALIDATIONCRYOEMPKLFILENAME=self.VALIDATIONCRYOEMPKLFILENAME)
        if Plugin.getPhenixVersion() >= PHENIXVERSION18:
            command += """
# Rama-Z (Ramachandran plot Z-score, RMSD))
if data.model.geometry.rama_z is not None:
    dictOverall2['Rama_Z_whole_n'] = "%d" % (data.model.geometry.rama_z.whole.n)
    if data.model.geometry.rama_z.whole.n != 0:
        dictOverall2['Rama_Z_whole_value'] = "%.2f" % abs(round(data.model.geometry.rama_z.whole.value, 2))
        dictOverall2['Rama_Z_whole_std'] = "%.2f" % (round(data.model.geometry.rama_z.whole.std, 2))
    else:
        dictOverall2['Rama_Z_whole_value'] = '---'
        dictOverall2['Rama_Z_whole_std'] = '---'
    dictOverall2['Rama_Z_helix_n'] = "%d" % (data.model.geometry.rama_z.helix.n)
    if data.model.geometry.rama_z.helix.n != 0:
        dictOverall2['Rama_Z_helix_value'] = "%.2f" % abs(round(data.model.geometry.rama_z.helix.value, 2))
        dictOverall2['Rama_Z_helix_std'] = "%.2f" % (round(data.model.geometry.rama_z.helix.std, 2))
    else:
        dictOverall2['Rama_Z_sheet_n'] = '---'
        dictOverall2['Rama_Z_sheet_value'] = '---'
    dictOverall2['Rama_Z_sheet_n'] = "%d" % (data.model.geometry.rama_z.sheet.n)
    if data.model.geometry.rama_z.sheet.n != 0:
        dictOverall2['Rama_Z_sheet_value'] = "%.2f" % abs(round(data.model.geometry.rama_z.sheet.value, 2))
        dictOverall2['Rama_Z_sheet_std'] = "%.2f" % (round(data.model.geometry.rama_z.sheet.std, 2))
    else:
        dictOverall2['Rama_Z_sheet_value'] = '---'
        dictOverall2['Rama_Z_sheet_std'] = '---'
    dictOverall2['Rama_Z_loop_n'] = "%d" % (data.model.geometry.rama_z.loop.n)
    if data.model.geometry.rama_z.loop.n != 0:
        dictOverall2['Rama_Z_loop_value'] = "%.2f" % abs(round(data.model.geometry.rama_z.loop.value, 2))
        dictOverall2['Rama_Z_loop_std'] = "%.2f" % (round(data.model.geometry.rama_z.loop.std, 2))
    else:
        dictOverall2['Rama_Z_loop_value'] = '---'
        dictOverall2['Rama_Z_loop_value'] = '---'
"""
        command += """with open('%s',"w") as f:
    f.write(json.dumps(dictOverall2))
""" % (VALIDATIONTMPFILENAME)

        pythonFileName = VALIDATIONTMPFILENAME.replace('.txt', '.py')
        # write script file
        with open(pythonFileName, "w") as f:
            f.write(command)

        # execute file with phenix.python
        Plugin.runPhenixProgram("", pythonFileName)

        # read file in scipion python
        with open(VALIDATIONTMPFILENAME, "r") as f:
            self.dictOverall2 = f.read()
            # self.dataDict = json.loads(f.read())

        self._store()








