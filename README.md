# Toxo

Toxo is an object-oriented MATLAB library for calculating penetrance tables of any bivariate epistasis model. It is solely dedicated to bivariate models, that is, models where the penetrance is an expression of two intervening variables. The user specifies the model, its desired heritability (or prevalence) and Minor Allele Frequency (MAF) and the library maximizes the resulting table's prevalence (or heritability).

Toxo does not generate population samples, and is intended to be used together with other already existing simulation software such as [GAMETES](https://sourceforge.net/projects/gametes/). In order to acomplish that, Toxo includes multiple output formats to provide compatibility with external simulators, and can be easily extended to support additional programs.

Toxo includes, as examples, two MATLAB scripts:

1. [calculate_tables.m](examples/calculate_tables.m): calculates the penetrance tables of the additive, multiplicative and threshold third and fourth order models, for a number of MAF and heritability combinations.

2. [usage_example.m](examples/usage_example.m): calculates a penetrance table of the 2-way additive model using MAF = 0.25 and h² = 0.2, saves the table using GAMETES format, and calls GAMETES with the previous table to generate population samples.

---

## Table of contents

1. [Requirements.](#1-requirements)

2. [Contents of this repository.](#2-contents-of-this-repository)

3. [Usage.](#3-usage)

4. [Class documentation.](#4-classes-documentation)

5. [Troubleshooting.](#5-troubleshooting)

6. [Acknowledgements.](#6-acknowledgements)

---

## 1. Requirements

- MATLAB (tested against version R2018a, it is likely to work on many others).

## 2. Contents of this repository

- `models/` [Marchini's models](https://doi.org/10.1038/ng1537) (generalized to 3rd-8th order) that can be used together with Toxo.

- `src/` Source code of the library, to be included in MATLAB' s path.

- `examples/` Example scripts showcasing the usage of Toxo.

- `LICENSE` MIT license.

- `README.md` This file.

## 3. Usage

### Installation

1. Download the latest Toxo release from [here](https://github.com/chponte/toxo/releases/latest).

2. Unzip the contents of the file.

3. Add the `src/` folder into your MATLAB environment:
   
   ```matlab
   addpath('path/to/src/folder');
   ```

### Using Toxo

Using Toxo starts with defining a model in CSV-like format. Not all models are supported by Toxo, make sure the desired model complies with the two requirements expressed below. Variable values are then calculated so that the resulting table fits some desired parameter, and the table is written into a file.

#### Model requirements

In order for Toxo to calculate the values of the two variables for which the prevalence (or heritability) is maximum, two requirements must be met:

1) Penetrance expressions are non-decreasing monotonic polynomials in the real positive number space. Polynomials that meet this criteria have positive partial derivatives for all real positive values of both variables. For example, the polynomial ![x(1+y)^2](https://latex.codecogs.com/gif.latex?x%281+y%29%5E2) is monotonically non-decreasing because ![d/dx(x(1+y)^2) = (1+y)^2](https://latex.codecogs.com/gif.latex?%5Ctfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%5Cbig%28x%281&plus;y%29%5E2%5Cbig%29%20%3D%20%281&plus;y%29%5E2) and ![d/dy(x(1+y)^2) = x(2y+2)](https://latex.codecogs.com/gif.latex?%5Ctfrac%7B%5Cpartial%7D%7B%5Cpartial%20y%7D%5Cbig%28x%281&plus;y%29%5E2%5Cbig%29%20%3D%20x%282y%20&plus;%202%29) are positive for x > 0 and y > 0.
2) Polynomials can be sorted unequivocally in the real positive number space. This can be demonstrated analytically for all models by comparing the polynomials in pairs. Example:
   
   ![x(1+y)^4 >= x(1+y)^3](https://latex.codecogs.com/gif.latex?x%281&plus;y%29%5E4%20%26%5Cge%20x%281&plus;y%29%5E3)
   
   ![1+y >= 1](https://latex.codecogs.com/gif.latex?1&plus;y%20%26%5Cge%201)
   
   ![y >= 0](https://latex.codecogs.com/gif.latex?y%20%26%5Cge%200)

#### Model description

Models read by Toxo are formatted in a CSV-like style, where each row represents a genotype combination and its associated penetrance (separated by a comma). The order in which the rows appear does not matter. Empty lines or lines starting with '#' (comments) are ignored.

Genotypes are represented as two characters, each one corresponding to each of the alleles from genotype. Alleles of the same genotype use the same alphabetic letter, and the difference in capitalization encodes the minor (lowercase) and major (uppercase) allele. There is no limit on the genotype combination size.

Penetrance expressions are functions of two variables. Variables can take any alphabetic name, but for simplicity we will name them x and y.

An example of a model would be:

```matlab
# 2-way additive model
AABB, x
AABb, x*(1+y)
AAbb, x*(1+y)^2
AaBB, x*(1+y)
AaBb, x*(1+y)^2
Aabb, x*(1+y)^3
aaBB, x*(1+y)^2
aaBb, x*(1+y)^3
aabb, x*(1+y)^4
```

#### Penetrance table calculation

Once the model is defined, obtaining the penetrance table is very straightforward. First the model is imported into MATLAB using the Model class. Then, the penetrance class is obtained. The last step is to save the calculated table into a file. The following code snippet exemplifies the process:

```matlab
model = 'sample_model.csv';
output = 'penetrance_table.csv';
mafs = [0.25, 0.25];
heritability = 0.1;
m = toxo.Model(model);
pt = m.find_max_prevalence(mafs, heritability);
pt.write(output, toxo.PTable.format_csv);
```

## 4. Classes documentation

Toxo implements two main classes, **Model** and **PTable**, which encapsulate all the functionality of the library:

| Model                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **+name**: *char*<br/>The name of the model as appears in its file name.<br/>**+order**: *double*<br/>    The number of locus invoved in the epistasis relationship.<br/>**+penetrances**: *sym*<br/>    Sorted array of symbolic penetrance expressions.<br/>**+variables**: *sym*<br/>    Array containing the symbolic variables used throughout the table.                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| **«constructor»+Model**(path: *char*)<br/>    Class constructor, responsible for reading the epistasis model from the given path.<br/>**-max_penetrance**(): *sym*<br/>    Returns the largest polynomial from all penetrance expressions, for any real and possitive value of the two variables.<br/>**-solve**(varargin: *sym*): *sym*<br/>    Solve the equation system made of the provided equations.<br/>**+find_max_prevalence**(mafs: *double*, h: *double*): *PTable*<br/>    Finds the table whose prevalence is maximum for the given MAFs and heritability, and returns it as a PTable object.<br/>**+find_max_heritability**(mafs: *double*, p: *double*): *PTable*<br/>    Finds the table whose heritability is maximum for the given MAFs and prevalence, and returns it as a PTable object. |

| PTable                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <u>**+format_csv**: *double* = 0</u><br/>    Static constant, storing the CSV format code for the write method.<br/><u>**+format_gametes**: *double* = 0</u><br/>    Static constant, storing the GAMETES format code for the write method.<br/>**+order**: *double*<br/>    The number of locus invoved in the epistasis relationship.<br/>**+vars**: *sym*<br/>    Struct containing the original variable names and its value.<br/>**+pt**: *sym*<br/>    Array of symbolic values representing the prenetrance table.                                                                                                                                                                                                                                                                                                              |
| **«constructor»+PTable**(model: *Model*, values: *double*)<br/>    Class constructor, creates a penetrance table from a Model and its variables values.<br/>**-to_gametes**(fmask: *char*, mafs: *double*): *char*<br/>    Returns a char array containing the table description, using the GAMETES format.<br/>**+prevalence**(mafs: *double*): *double*<br/>    Returns the prevalence of the table.<br/>**+heritability**(mafs: *double*): *double*<br/>    Returns the heritability of the table.<br/>**+marginal_penetrances**(mafs: *double*): *double*<br/>    Returns the marginal penetrance of each allele from every locus.<br/>**+write**(path: *char*, format: *double*, varargin)<br/>    Writes the table to the specified file using a specific format. Depending on the format, additional arguments may be provided. |

## 5. Troubleshooting

> Why do I get "*There is no solution to the problem defined*" when finding a penetrance table?

The MAF and prevalence (or heritability) combination that you have specified leads to an incompatible system. Make sure that that prevalence or heritability level can be achieved for the specified MAF. Plotting your model's prevalence or heritability surface can help to understand it.

> Why do I get "*Could not find a solution to the problem defined*" when finding a penetrance table?

The model that you are specifying may be too complex for the MATLAB's solver to process. This may be the case if, for example, your model contains big exponents (>32).

If you are having trouble using Toxo, encounter any error or would like to see some additional functionality implemented, feel free to open an [Issue](https://github.com/chponte/toxo/issues).

## 6. Acknowledgements

Thanks to María J. Martín and Jorge González-Domínguez for their supervision of the project and writting of the article.
