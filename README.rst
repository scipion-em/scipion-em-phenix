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

  * scipion test phenix.tests.test_protocol_dock_in_map
  * scipion test phenix.tests.test_protocol_emringer
  * scipion test phenix.tests.test_protocol_real_space_refine
  * scipion test phenix.tests.test_protocol_molprobity
  * scipion test phenix.tests.test_protocol_superpose_pdbs
  * scipion test phenix.tests.test_protocol_validation_cryoem


- **Supported versions of PHENIX**

Tested with  `Phenix-1.13-2998, Phenix-1.16-3549 and Phenix- 1.17.1`




=========
Protocols
=========

* emringer: Validates the agreement between the initial map and the derived low-resolution atomic structure. This program samples the density around Chi1 angles of protein sidechains. Electronic density and appropriate rotameric angles must coincide for each residue if the atomic structure backbone has been perfectly fitted to the map.
* molprobity: Validates the geometry of an atomic structure inferred from an electron density map.
* real_space_refine: Designed for extensive real-space refinement of an atomic structure against the map provided. The map can be derived from X-ray or neutron crystallography, or cryoEM. The program obtains a model that fits the map as well as possible having appropriate geometry. The model should not show validation outliers, such as Ramachandran plot or rotamer outliers.
* superpose_pdbs: Superposes two atomic structures so that they optimally match.
* validation_cryoem: generalization of molprobity implemented by Phenix package.


========
Examples
========

See `Model Building Tutorial <https://github.com/I2PC/scipion/wiki/tutorials/tutorial_model_building_basic.pdf>`_

  
===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/phenix_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/phenix_prod.svg

  
