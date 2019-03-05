=================
Scipion-em-phenix
=================

Wrapper to use Phenix in Scipion. Supports `Phenix-1.13 <https://www.phenix-online.org/download/nightly_builds.cgi>`_ . So far we have implemented:

  * emringer
  * real_space_refinement
  * molprobity
  * superpose pdbs

Installation
============

1. Install `Phenix-1.13 <https://www.phenix-online.org/download/nightly_builds.cgi>`_

2. Add the path to ~/.config/scipion/scipion.conf (yours might look different, please check your Phenix path).
   ``PHENIX_HOME=/usr/local/phenix-1.13-2998``
   
3. Install scipion-em-phenix
  a) In user mode: 
  ``scipion installp -p scipion-em-phenix``

  b) In developer mode: 
  ``scipion installp -p /path/to/scipion-em-phenix --devel``
  
  TESTS:
  
  * scipion test phenix.tests.test_protocol_emringer
  * scipion test phenix.tests.test_protocol_real_space_refine
  * scipion test phenix.tests.test_protocol_molprobity
  * scipion test phenix.tests.test_protocol_superpose_pdbs

  
![alt text](http://arquimedes.cnb.csic.es:9980/badges/phenix_devel.svg)

Builbot_URL: http://arquimedes.cnb.csic.es:9980/#/builders/25
  
  
