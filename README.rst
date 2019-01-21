============
Toxo
============

Toxo is a MATLAB library which provides a powerful, object-oriented implementation of penetrance table calculations for epistasis interaction models. Toxo calculates the penetrance table for a given epistasis model that maximizes the heritability for a given prevalence and Minor Allele Frequency (MAF). Maximizing the prevalence is an analogous process.


Requirements
------------------

* MATLAB R2018a (correct behaviour is not guaranteed in previous or later versions, although it is likely to work).


Installation
------------------

1) Download the latest Toxo release from `here <https://github.com/chponte/toxo/releases/latest>`__.
2) Unzip the contents of the file.
3) Add the ``src/`` folder into your MATLAB environment or script:

   .. code:: matlab

      addpath('path/to/src/folder');


Usage
------------------

Using Toxo starts with defining a model in CSV-like format. Not all models are supported by Toxo, make sure the desired model complies with the two conditions. Calculating the variable values is then performed, and the resulting penetrance table is written into a text file. A complete working example can be found `here <https://github.com/chponte/toxo/blob/master/generate_models.m>`__.

Model requirements
^^^^^^^^^^^^^^^^^^
In order for Toxo to calculate the values of the two variables for which the prevalence (or heritability) is maximum, two requirements must be met:

1) Penetrance expressions are non-decreasing monotonic polynomials in the real positive number space. Polynomials that meet this criteria have positive partial derivatives for all real positive values of both variables.
2) Polynomials can be sorted unequivocally in the real positive number space. This can be demonstrated analytically for all models by comparing the polynomials in pairs. As an example, the demonstration that |e1| is greater than |e2| for real positive numbers would be:

.. |e1| image:: https://latex.codecogs.com/gif.latex?x%281&plus;y%29%5E4
   :align: bottom

.. |e2| image:: https://latex.codecogs.com/gif.latex?x%281&plus;y%29%5E3
   :align: bottom

.. figure:: https://latex.codecogs.com/gif.latex?x%281&plus;y%29%5E4%20%26%5Cge%20x%281&plus;y%29%5E3
   :align: center

.. figure:: https://latex.codecogs.com/gif.latex?%281&plus;y%29%5E4%20%26%5Cge%20%281&plus;y%29%5E3
   :align: center

.. figure:: https://latex.codecogs.com/gif.latex?1&plus;y%20%26%5Cge%201
   :align: center

.. figure:: https://latex.codecogs.com/gif.latex?y%20%26%5Cge%200
   :align: center


Model description
^^^^^^^^^^^^^^^^^
Models read by Toxo are formatted in CSV-like style, where each row represents a genotype combination and, separated by a comma, its associated penetrance. The order in which the rows appear does not matter. Empty lines or lines starting with '#' (comments) are ignored.

Genotypes are represented as two characters, each one corresponding to each of the alleles from genotype. Alleles of the same genotype use the same alphabetic letter, and the difference in capitalization encodes the minor (lowercase) and major (uppercase) allele. There is no limit on the genotype combination size.

Penetrance expressions are functions of two variables. Variables can take any alphabetic name, but for simplicity we will name them x and y.

An example of model would be:

.. code:: text
   
   # Model name
   
   AABB, x
   AABb, x
   AAbb, x
   AaBB, x
   AaBb, x*(1+y)
   Aabb, x*(1+y)
   aaBB, x
   aaBb, x*(1+y)
   aabb, x*(1+y)

Penetrance table calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once the model is defined, obtaining the penetrance table is very straightforward. First the model is imported into MATLAB using the Model class. Then, the penetrance class is obtained. The last step is to save the calculated table into a file. The following code snippet exemplifies the process:


.. code:: matlab
   
   model = 'sample_model.csv';
   output = 'penetrance_table.csv';
   maf = 0.25;
   prevalence = 0.1;
   m = toxo.Model(model);
   p = m.find_max_heritability(maf, prevalence);
   p.write(output, toxo.PTable.format_csv);


Classes in Toxo
------------------
Toxo implements two main classes, Model_ and PTable_, which encapsulate all the functionality:

Model
^^^^^^^^^^^^
Model is a symbolic representation of a epistasis model. It is responsible for reading the model, parsing the text file and converting the penetrance strings to symbolic expressions and calculating the penetrance table 

Attributes
""""""""""
name : ``String``
  Name of the model.
order : ``Integer``
  Number of loci involved in the epistatic model.
penetrances : ``Array of symbolic``
  Array of symbolic expressions, representing the epistatic model.
variables : ``Array of symbolic``
  List of all variables contained in all symbolic expressions

Methods
""""""""""
Model(path)
  Construct an instance of this class from the given model.
  
  - ``path`` : ``String`` - Path to the model CSV file.
find_max_prevalence(maf, h)
  Calculate the penetrance table(s) of the model with the maximum admissible prevalence given its MAF and heritability.
  
  - ``maf`` : ``Double`` - MAF of the resulting penetrance table.
  - ``h``: ``Double`` - Heritability of the resulting penetrance table.
  - ``output`` : ``toxo.PTable`` - Resulting penetrance table.
find_max_heritability(maf, p)
  Calculate the penetrance table(s) of the model with the maximum admissible heritability given its MAF and prevalence.
  
  - ``maf``: ``Double`` - MAF of the resulting penetrance table.
  - ``p``: ``Double`` - Prevalence of the resulting penetrance table.
  - ``output`` : ``toxo.PTable`` - Resulting penetrance table.
  
PTable
^^^^^^^^^^^^
Static constants
""""""""""""""""""""
format_csv : ``Integer``
  Represents the CSV output format, taken as a parameter in the write method.
format_gametes: ``Integer``
  Represents the GAMETES output format, taken as a parameter in the write method.

Attributes
""""""""""
order : ``Integer``
  Number of loci involved in the penetrance table.
maf : ``Double``
  Common MAF of all locis involved in the interaction.
vars : ``Map``
  Values of the variables present in the original model.
gp : ``Array of symbolic``
  Genotype probabilities table array.
pt : ``Array of symbolic``
  Penetrances table array.

Methods
""""""""""
PTable(model, maf, values)
  Create a penetrance table from a given Model, using the MAF and variable values desired.
  
  - ``model``: ``toxo.Model`` - Model from which the table is constructed.
  - ``maf``: ``Double`` - MAF of the penetrance table.
  - ``values``: ``Array of double`` - Values of the variables in Model.
prevalence( )
  Calculate the prevalence of the penetrance table.
  
  - ``output`` : ``Double`` - Prevalence of the table.
heritability( )
  Calculate the heritability of the penetrance table.
  
  - ``output`` : ``Double`` - Heritability of the table.
write(path, format)
  Write the penetrance table into a text file using a specific output format.
  
  - ``path``: ``String`` - File path in which the table should be written into.
  - ``format``: ``Integer`` - Format to use for the output.

Troubleshooting
------------------

If you are having trouble using Toxo, encounter any error or would like to see some additional functionality implemented, feel free to open an `Issue <https://github.com/chponte/toxo/issues>`_.

References
------------------

.. [1] Marchini, Jonathan, Peter Donnelly, and Lon R. Cardon. 2005. "Genome-Wide Strategies for Detecting Multiple Loci That Influence Complex Diseases". Nature Genetics 37 (4): 413. https://doi.org/10.1038/ng1537.
