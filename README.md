# Toxo

Toxo is an object-oriented MATLAB library for calculating penetrance tables of any bivariate epistasis model. It is solely dedicated to bivariate models, that is, models where the penetrance is an expression of two intervening variables. The user specifies the model, its desired heritability (or prevalence) and Minor Allele Frequency (MAF) and the library maximizes the resulting table's prevalence (or heritability).

Toxo does not generate population samples, and is intended to be used together with other already existing simulation software such as [GAMETES](https://sourceforge.net/projects/gametes/). In order to acomplish that, Toxo includes multiple output formats to provide compatibility with external simulators, and can be easily extended to support additional programs.

Toxo includes, as examples, two MATLAB scripts:

1. [calculate_tables.m](http): calculates the penetrance tables of the additive, multiplicative and threshold third and fourth order models, for a number of MAF and heritability combinations.

2. [usage_example.m](http): calculates a penetrance table of the 2-way additive model using MAF = 0.25 and h² = 0.2, saves the table using GAMETES format, and calls GAMETES with the previous table to generate population samples.

## Table of contents

1. [Requirements.](https://github.com/chponte/toxo#requirements)

2. [Contents of this repository.](https://github.com/chponte/toxo#contents-of-this-repository)

3. [Usage.](https://github.com/chponte/toxo#usage)

4. [Class documentation.](https://github.com/chponte/toxo#class-documentation)

5. [Troubleshooting.](https://github.com/chponte/toxo#troubleshooting)

6. [References.](https://github.com/chponte/toxo#references)

## 1. Requirements

- MATLAB (tested against version R2018a, it is likely to work on many others).

## 2. Contents of this repository

- `models/` Sample models that can be used together with Toxo.

- `src/` Source code of the library, to be included in MATLAB' s path.

- `examples/` Example scripts showcasing the usage of Toxo.

- `LICENSE` MIT license.

- `README.md` This file.

## 3. Usage

### Installation

1. Download the latest Toxo release from `here <https://github.com/chponte/toxo/releases/latest>`__.

2. Unzip the contents of the file.

3. Add the `src/` folder into your MATLAB environment.: .. code:: matlab addpath('path/to/src/folder');

### Using Toxo

Using Toxo starts with defining a model in CSV-like format. Not all models are supported by Toxo, make sure the desired model complies with the two requirements expressed below. Model variable values are then calculated, and the resulting penetrance table is written into a text file.

#### Model requirements

In order for Toxo to calculate the values of the two variables for which the prevalence (or heritability) is maximum, two requirements must be met:

1) Penetrance expressions are non-decreasing monotonic polynomials in the real positive number space. Polynomials that meet this criteria have positive partial derivatives for all real positive values of both variables. For example, the polynomial ![x(1+y)^2](https://latex.codecogs.com/gif.latex?x%281+y%29%5E2) is monotonically non-decreasing because ![d/dx(x(1+y)^2) = (1+y)^2](https://latex.codecogs.com/gif.latex?%5Ctfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%5Cbig%28x%281&plus;y%29%5E2%5Cbig%29%20%3D%20%281&plus;y%29%5E2) and ![d/dy(x(1+y)^2) = x(2y+2)](https://latex.codecogs.com/gif.latex?%5Ctfrac%7B%5Cpartial%7D%7B%5Cpartial%20y%7D%5Cbig%28x%281&plus;y%29%5E2%5Cbig%29%20%3D%20x%282y%20&plus;%202%29) are positive for x > 0 and y > 0.
2) Polynomials can be sorted unequivocally in the real positive number space. This can be demonstrated analytically for all models by comparing the polynomials in pairs. Example:
   
   ![x(1+y)^4 >= x(1+y)^3](https://latex.codecogs.com/gif.latex?x%281&plus;y%29%5E4%20%26%5Cge%20x%281&plus;y%29%5E3)
   
   ![1+y >= 1](https://latex.codecogs.com/gif.latex?1&plus;y%20%26%5Cge%201)
   
   ![y >= 0](https://latex.codecogs.com/gif.latex?y%20%26%5Cge%200)

#### Model description

Models read by Toxo are formatted in CSV-like style, where each row represents a genotype combination and its associated penetrance (separated by a comma). The order in which the rows appear does not matter. Empty lines or lines starting with '#' (comments) are ignored.

Genotypes are represented as two characters, each one corresponding to each of the alleles from genotype. Alleles of the same genotype use the same alphabetic letter, and the difference in capitalization encodes the minor (lowercase) and major (uppercase) allele. There is no limit on the genotype combination size.

Penetrance expressions are functions of two variables. Variables can take any alphabetic name, but for simplicity we will name them x and y.

An example of model would be:

```matlab
# 2-way additive model
AABB, x
AABb, x
AAbb, x
AaBB, x
AaBb, x*(1+y)
Aabb, x*(1+y)
aaBB, x
aaBb, x*(1+y)
aabb, x*(1+y)   
```

#### Penetrance table calculation

Once the model is defined, obtaining the penetrance table is very straightforward. First the model is imported into MATLAB using the Model class. Then, the penetrance class is obtained. The last step is to save the calculated table into a file. The following code snippet exemplifies the process:

```matlab
model = 'sample_model.csv';
output = 'penetrance_table.csv';
mafs = [0.25, 0.25];
heritability = 0.1;
m = toxo.Model(model);
p = m.find_max_prevalence(mafs, heritability);
p.write(output, toxo.PTable.format_csv);
```

## 4. Classes documentation

Toxo implements two main classes, Model_ and PTable_, which encapsulate all the functionality:

Model

Model is a symbolic representation of an epistasis model. It is responsible for reading the model, parsing the text file and converting the penetrance strings to symbolic expressions. It offers two methods to calculate penetrance tables which maximize the associated penetrance or heritability under certain constraints.

<div>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#ccc;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#ccc;color:#333;background-color:#fff;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#ccc;color:#333;background-color:#f0f0f0;}
.tg .tg-s268{text-align:left}
</style>
<table class="tg">
  <tr>
    <th class="tg-s268">Model<br></th>
  </tr>
  <tr>
    <td class="tg-s268">+name<br>+order<br>+penetrances<br>+variables<br></td>
  </tr>
  <tr>
    <td class="tg-s268">metodos</td>
  </tr>
</table>
</div>

Attributes

name : ``String``
  Name of the model.
order : ``Integer``
  Number of loci involved in the epistatic model.
penetrances : ``Array of symbolic``
  Array of symbolic expressions, representing the epistatic model.
variables : ``Array of symbolic``
  List of all variables contained in all symbolic expressions

Methods

Model(path)
  Construct an instance of this class from the given model.

- ``path`` : ``String`` - Path to the model CSV file.
  
  find_max_prevalence(maf, h)
  
  Calculate the penetrance table(s) of the model with the maximum admissible prevalence given its MAF and heritability.

- ``maf`` : ``Double`` - MAF of the resulting penetrance table.

- ``h`` : ``Double`` - Heritability of the resulting penetrance table.

- ``output`` : ``toxo.PTable`` - Resulting penetrance table.
  
  find_max_heritability(maf, p)
  
  Calculate the penetrance table(s) of the model with the maximum admissible heritability given its MAF and prevalence.

- ``maf`` : ``Double`` - MAF of the resulting penetrance table.

- ``p`` : ``Double`` - Prevalence of the resulting penetrance table.

- ``output`` : ``toxo.PTable`` - Resulting penetrance table.

PTable

Numeric representation of a penetrance table. This class provides methods to calculate several metrics, as well as a method to write the table to a file in several formats.

Static constants
"""""""""""""""""""""""""""""""""""""
format_csv : ``Integer``
  Represents the CSV output format, taken as a parameter in the write method.
format_gametes: ``Integer``
  Represents the GAMETES output format, taken as a parameter in the write method.

Attributes
"""""""""""""""""""""""""""""""""""""
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
"""""""""""""""""""""""""""""""""""""
PTable(model, maf, values)
  Create a penetrance table from a given Model, using the MAF and variable values desired.

- ``model`` : ``toxo.Model`` - Model from which the table is constructed.

- ``maf`` : ``Double`` - MAF of the penetrance table.

- ``values`` : ``Array of double`` - Values of the variables in Model.
  
  prevalence( )
  
  Calculate the prevalence of the penetrance table.

- ``output`` : ``Double`` - Prevalence of the table.
  
  heritability( )
  
  Calculate the heritability of the penetrance table.

- ``output`` : ``Double`` - Heritability of the table.
  
  write(path, format)
  
  Write the penetrance table into a text file using a specific output format.

- ``path`` : ``String`` - File path in which the table should be written into.

- ``format`` : ``Integer`` - Format to use for the output.

Troubleshooting
-------------------------------------

If you are having trouble using Toxo, encounter any error or would like to see some additional functionality implemented, feel free to open an `Issue <https://github.com/chponte/toxo/issues>`_.

References
-------------------------------------

.. [1] Marchini, Jonathan, Peter Donnelly, and Lon R. Cardon. 2005. "Genome-Wide Strategies for Detecting Multiple Loci That Influence Complex Diseases". Nature Genetics 37 (4): 413. https://doi.org/10.1038/ng1537.
