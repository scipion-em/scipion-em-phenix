=====================
Phenix scipion plugin
=====================

This plugin allows to use programs from the *PHENIX* software suite within the Scipion framework. **You need to install the Phenix suite before installing the plugin**, see section "Binary Files" for details.

Phenix is a software suite that allows model building of macromolecule structures obtained by X-ray crystallography, and that has been extended to other techniques like cryo-EM (see `Phenix home page <https://www.phenix-online.org/>`_ for details).

Current programs implemented:

  * dock in map
  * emringer
  * real space refine
  * molprobity
  * superpose pdbs
  * validation cryoem
  * search fit

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-em-phenix
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-phenix.git

2. Install:

.. code-block::

           scipion installp -p path_to_scipion-em-phenix --devel
 
 
- **Binary files** 

*PHENIX* binaries will *NOT* be installed automatically with the plugin. The independent installation of PHENIX software suite by the user is required before running the programs. Default installation path assumed is */usr/local/phenix-1.13-2998*; this path or any other of your preference has to be set in *PHENIX_HOME* in *scipion.conf* file. We recommend to install PHENIX version 1.13-2998.

  the plug-in also requires imagemagick package:  sudo apt-get install imagemagick

- **Tests**

Tested with PHENIX version: 1.13-2998.

To check the installation, simply run the following Scipion tests: 

  * scipion test --grep phenix --run 


- **Supported versions of PHENIX**

Tested with  `Phenix-1.13-2998, Phenix-1.16-3549, Phenix- 1.17.1, Phenix 1.18.2, phenix-1.19.2 and phenix 1.20.1`




=========
Protocols
=========

* emringer: Validates the agreement between the initial map and the derived low-resolution atomic structure. This program samples the density around Chi1 angles of protein sidechains. Electronic density and appropriate rotameric angles must coincide for each residue if the atomic structure backbone has been perfectly fitted to the map.
* molprobity: Validates the geometry of an atomic structure inferred from an electron density map.
* real_space_refine: Designed for extensive real-space refinement of an atomic structure against the map provided. The map can be derived from X-ray or neutron crystallography, or cryoEM. The program obtains a model that fits the map as well as possible having appropriate geometry. The model should not show validation outliers, such as Ramachandran plot or rotamer outliers.
* superpose_pdbs: Superposes two atomic structures so that they optimally match.
* validation_cryoem: generalization of molprobity implemented by Phenix package.
* search_fit: given a chain of n alanines, a 3D map and a sequence search for the subsequence of n aminoacids that better fits in the density. Only works if the atomic structure has a single chain.
* rebuild_docked_predicted_alphafold2_model: Rebuild predicted model morphs and rebuilds a model produced by AlphaFold,
     RoseTTAFold and other prediction software into a cryo EM map, using a set
     of docked domains from the predicted model as a template.
* protocol_dock_in_map: Docking of a PDB (one or several copies) into a map
* dock_and_rebuild_alphafold_model: Rebuild predicted model morphs and rebuilds a model produced by AlphaFold, RoseTTAFold and other prediction software into a cryo EM map, using a set of docked domains from the predicted model as a template.
* protocol_process_predicted_alphafold2_model: Replace values in b-factor field with estimated B values. Optionally remove low-confidence residues and split into domains.
* dock_predicted_alphafold2_modeldocks the domains from a model produced by AlphaFold, RoseTTAFold and other prediction software into a cryo EM map. It uses the connectivity of the model as a restraint in the docking process so that the docked domains normally are in a reasonable arrangement. It can take map symmetry into account.


========
Examples
========

See `Model Building Tutorial <https://scipion-em.github.io/docs/release-3.0.0/docs/user/user-documentation.html#model-building>`_

  
===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/phenix_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/phenix_prod.svg

  
