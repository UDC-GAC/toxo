============
Toxo
============

Toxo is a MATLAB library which provides a powerful, object-oriented implementation of penetrance table computations for epistasis interaction models.

Requirements
------------------

* MATLAB R2018a (correct behaviour is not guaranteed in previous or later versions, although it is likely to work).

Usage
------------------

1) Download the latest Toxo release from here_.
2) Unzip the contents of the file.
3) Add the ``src/`` folder into MATLAB path:

   .. code-block:: matlab

      addpath("path/to/src/folder");

4) *(Optional)* Import library contents for ease of use: 

   .. code-block:: matlab
   
      import toxo.*;

.. _here: https://github.com/chponte/toxo/releases/latest

Troubleshooting
------------------

If you are having trouble using Toxo, encounter any error or would like to see some additional functionality implemented, feel free to open an Issue_.

.. _Issue: https://github.com/chponte/toxo/issues

References
------------------

.. [1] Marchini, Jonathan, Peter Donnelly, and Lon R. Cardon. 2005. ‘Genome-Wide Strategies for Detecting Multiple Loci That Influence Complex Diseases’. Nature Genetics 37 (4): 413–17. https://doi.org/10.1038/ng1537.
