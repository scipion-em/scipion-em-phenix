=====================
Phenix scipion plugin
=====================

Wrapper to use Phenix in Scipion. Supports `Phenix-1.13 <https://www.phenix-online.org/download/nightly_builds.cgi>`_ . So far we have implemented:

  - emringer
  - real_space_refinement
  - molprobity
  - superpose pdbs


============
INSTALLATION
============

- **Install Phenix-1.13** 

    <https://www.phenix-online.org/download/nightly_builds.cgi>

- **Add the path to ~/.config/scipion/scipion.conf** (yours might look different, please check your Phenix path).
    
.. code-block::

    PHENIX_HOME=/usr/local/phenix-1.13-2998
   
- **Install this plugin:**

.. code-block::

    scipion installp -p scipion-em-powerfit

OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

Alternatively, in devel mode:

.. code-block::

    scipion installp -p local/path/to/scipion-em-powerfit --devel

- **TESTS:**
  
  - scipion test phenix.tests.test_protocol_emringer
  - scipion test phenix.tests.test_protocol_real_space_refine
  - scipion test phenix.tests.test_protocol_molprobity
  - scipion test phenix.tests.test_protocol_superpose_pdbs

  
![alt text](http://arquimedes.cnb.csic.es:9980/badges/phenix_devel.svg)


Builbot_URL: http://arquimedes.cnb.csic.es:9980/#/builders/25
  
  
