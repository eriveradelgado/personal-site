---
title: "Chemically Aware Tables in R"
author: "Edgardo Rivera-Delgado"
date: '2021-06-30'
slug: chemistry-in-rstudio-with-RDKit-and-reticulate
categories: 
 - chemistry
 - library-exploration
tags: [reticulate, RDKit, flextable]
---

<link href="{{< blogdown/postref >}}index_files/tabwid/tabwid.css" rel="stylesheet" />
<link href="{{< blogdown/postref >}}index_files/tabwid/tabwid.css" rel="stylesheet" />

# Motivation

Often in pharmaceutical analysis it helps to have the structure of the molecule
alongside key properties one is interested in evaluating. A table displaying
the structure with some molecule-derived information is often the preferred format.
In this post I explore how to use R to achieve this goal. At the same time,
we’ll explore how Python and R can be used together, as promised in the [last post](https://riveradelgado.com/post/2020/04/18/chemistry-in-r/) where I
examined how R can be used for input display and writing to file of two common
chemical formats with packages that are native to R.

This time I’ll make use of the `RDKit` library, a Python library that has been of great
interest for me to explore for a while. The library `reticulate` will helps us
call Python libaries and use them within R without having to leave RStudio.

# Reticulate Setup

To load reticulate some setup is necessary. Systems can have multiple versions
of Python installed which can make selection a little difficult. From the
different ways possible, I found that using reticulate with mini-conda the most
straightforward way of getting started. When installing reticulate the package
will recommend to create a mini-conda installation. I suggest you do this
and use the py_install functions in `reticulate` to install `RDKit`.

Remember to make sure you have reticulate installed by calling the install
function.

    install.packages("reticulate")

Once `reticulate` is installed, the necessary packages need to be loaded and
the Python environment needs to be set up to be used for the reticulate session.

``` r
library(reticulate)
library(tidyverse)
library(flextable)
library(magrittr)
library(data.table)
library(magrittr)
reticulate::use_condaenv("r-reticulate")
```

Installing RDKit can be achieved by the `py_install()` function as shown below.

    reticulate::py_install("rdkit")

With setup out of the way reticulate is ready to run Python code in RStudio and
we are on our way to building tables with chemical structures in them.

## Testing Reticulate Works

There are two ways to speak Python with reticulate. Python code can be written
directly into the notebook or Python functions can be imported into R and used
as regular R functions.

Here’s my attempt at writing Python straight from R using a matplotlib example.
In order for the code chunk to run we start a new chunk as we would usually do
and embrace it in curly parenthesis with {python} instead of {r}. Python code
can then be written as it would normally be written in a jupyter notebook or
any other Python editor. For `matplotlib` the results will be shown inline as
in any other chunk written in R.

``` python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
x = np.linspace(1,10,100)
y = x**2
fig = plt.figure(figsize=(14,8))
plt.plot(x,y)
plt.show()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-1-1.png" width="1344" />

Other png outputs from Python like the PIL image format may need to be saved to
file first and then read into R via `knitr::include_graphics()` or via Rmarkdown
image referencing mechanism.

# Reticulated Chemistry

In the spirit of getting an output as fast as possible let’s create a display
function using `RDKit` that can show molecules inline. The function and logic
are modified from the nice implementation developed by [Lars Richter](https://github.com/richterlars/larspack).

``` r
# Assign needed RDKIT python libraries to r functions via reticuale::import()
Chem <- reticulate::import("rdkit.Chem") 
Draw <- reticulate::import("rdkit.Chem.Draw")
Plt <- reticulate::import("matplotlib.pyplot")
AllChem <- reticulate::import("rdkit.Chem.AllChem")
# Use caffeine as an example for molecular display
caffeine <- "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

display_structure <- function(mol){
    molecule <- Chem$MolFromSmiles(mol)
    filepath <- fs::path(mol, ext = "png")
    Draw$MolToFile(molecule, filepath) 
    knitr::include_graphics(filepath)
}

display_structure(mol = caffeine)
```

<img src="CN1C=NC2=C1C(=O)N(C(=O)N2C)C.png" width="150" />

Notice above how the function `Chem <- reticulate::import("rdkit.Chem")` is
replicating the `import pandas as pd` role in the python code chunk. This
allows to store the functions of that particular RDKit library in an R object.
Each one of Chem’s functions can then be accessed with the dollar `$` operator.
In the previous chunk, I am importing several of the libraries from RDKit that
are used for displaying chemical structures.

## Displaying Mutiple Molecules

To display multiple molecules let’s first create a vector containing our
gastronomy inspired set of example molecules.

``` r
multiple_smiles <- c(
  paste0("O=C(N(C(C1=C(C(OC([H])([H])[H])=C(O[H])C(=C1[H])[H])[H])",
         "([H])[H])[H])C(C(C(C(C(=C(C(C([H])([H])[H])(C([H])([H])[H])[H])",
         "[H])[H])([H])[H])([H])[H])([H])[H])([H])[H]"
         ), # Chili 
         "CCO", # Rum
         "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" # Coffee
                 ) 
chemical <-   c("Capsaicin", "Ethanol", "Caffeine")
source   <-   c("Chili", "Rum", "Coffee")
```

Displaying a molecule grid with RDKIT in Rmarkdown can similarly be achieved
with the following function.

``` r
display_grid <- function(multiple_smiles,
                         size = 800, 
                         nrow = length(multiple_smiles)
                         ){
  
# mapping with python functions
mol_list <- purrr::map(multiple_smiles, Chem$MolFromSmiles)

filepath_grid <- fs::path("plot_grid", ext = "png")

img <- Draw$MolsToGridImage(mol_list,
                            molsPerRow = as.integer(nrow),
                            # Set the size of each molecule in the grid
                            # as.integer() is a necessary part for python to interpret
                            # correctly the input
                            subImgSize = c(as.integer(size), as.integer(size))
                            )
  
    img$save(filepath_grid)
 knitr::include_graphics(filepath_grid)
 
}
display_grid(multiple_smiles)
```

<img src="plot_grid.png" width="1200" />

Notice how we can use `purrr::map` with a Python function that has been imported
to R as we would with any native R function. Since the the RDKit libraries are
loaded in the global environment they are accessible to all the
functions in this document. For a more robust implementation, assignment from
the Python libraries should probably be done inside of each function.

There are some peculiarities with working with Python code from R. For example,
when calling `Draw$MolsToGridImage` the `molsPerRow` argument requires its input
to be a strict integer and does not coerce a double input to integer
automatically. Therefore we need to convert its input to integer for the
function to work properly.

For more details on reticulate’s capacity to intercovert between R and python’s
structures I found reticulate’s [cheat sheet](https://ugoproto.github.io/ugo_r_doc/pdf/reticulate.pdf)
to be a good reference.

# Chemical Tables

To achieve our intended goal let’s use RDKit to create a function that converts
SMILES and chemical names into chemical images saved to a directory. Then we
can read the images from the directory into the [flextable](https://davidgohel.github.io/flextable/)
library to build our chemically aware table.

``` r
save_structure <- function(smiles, chemical, path){

Chem <- reticulate::import("`RDKit`.Chem")
Draw <- reticulate::import("`RDKit`.Chem.Draw")
Plt <-  reticulate::import("matplotlib.pyplot")
AllChem <- reticulate::import("`RDKit`.Chem.AllChem")

mol <- Chem$MolFromSmiles(smiles)

filepath <- fs::path(path, chemical, ext = "png")
# Unlike the functions before this function just saves to file without calling
# knitr::include_graphics()
Draw$MolToFile(mol = mol, filename = filepath) 
}
```

We’ll loop over the SMILES and Chemical names with the help of purrr’s `walk2`
function which allows for two vectors to be fed into a function simultaneously.

``` r
fs::dir_create("molecules")

table_chemicals <- 
  # Allocate vectors to data frame
  data.frame(
  SMILES = multiple_smiles, 
  Chemical = chemical,
  Source = source, 
  stringsAsFactors = FALSE) %>% 
  # Loop over the smiles and chemical names to save to file
  # Use possibly to catch any errors if there are any.
  mutate(
    savemol = purrr::walk2(
    SMILES, Chemical,
    purrr::possibly(~save_structure(
      smiles = .x,
      chemical = .y,
      path = "./molecules"),
      otherwise = "oh-no!"
      )
    )
    ) %>% 
  mutate(
    # Record the path where images are stored. Convert to character the path
    # to allow flextable to handle the paths.
    structure_path = as.character(fs::path("molecules", Chemical, ext = "png")),
    # Leave one column empty for the eventual location of the image
    Structure = "",
    SMILES = stringr::str_wrap(SMILES, width = 34)
  ) 
```

## Build the Flextable

To build the flextable I will be following along the gallery example
for [tennis players](https://ardata-fr.github.io/flextable-gallery/gallery/2021-03-26-tennis-players/).

The first step is to initialize the flextable by using the `flextable()`
function. The columns to be displayed are provided as a character vector to the
argument `col_keys`.

``` r
table_chemicals_ft <- 
  flextable::flextable(
    table_chemicals,
    col_keys = c(
       "SMILES", "Chemical", "Source", "Structure"
      )
    )
table_chemicals_ft
```

<template id="646b0ffc-711a-4f4c-881c-52eb506e76b0"><style>
.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-cf2bb2fe{}.cl-cf206b60{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-cf209374{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-cf20d960{width:54pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf20d96a{width:54pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf20d974{width:54pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-cf2bb2fe'>
<thead><tr style="overflow-wrap:break-word;"><td class="cl-cf20d974"><p class="cl-cf209374"><span class="cl-cf206b60">SMILES</span></p></td><td class="cl-cf20d974"><p class="cl-cf209374"><span class="cl-cf206b60">Chemical</span></p></td><td class="cl-cf20d974"><p class="cl-cf209374"><span class="cl-cf206b60">Source</span></p></td><td class="cl-cf20d974"><p class="cl-cf209374"><span class="cl-cf206b60">Structure</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60">O=C(N(C(C1=C(C(OC([H])([H])<br>[H])=C(O[H])C(=C1[H])[H])[H])([H])<br>[H])[H])C(C(C(C(C(=C(C(C([H])([H])<br>[H])(C([H])([H])[H])[H])[H])[H])<br>([H])[H])([H])[H])([H])[H])([H])<br>[H]</span></p></td><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60">Capsaicin</span></p></td><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60">Chili</span></p></td><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60">CCO</span></p></td><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60">Ethanol</span></p></td><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60">Rum</span></p></td><td class="cl-cf20d960"><p class="cl-cf209374"><span class="cl-cf206b60"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-cf20d96a"><p class="cl-cf209374"><span class="cl-cf206b60">CN1C=NC2=C1C(=O)N(C(=O)N2C)C</span></p></td><td class="cl-cf20d96a"><p class="cl-cf209374"><span class="cl-cf206b60">Caffeine</span></p></td><td class="cl-cf20d96a"><p class="cl-cf209374"><span class="cl-cf206b60">Coffee</span></p></td><td class="cl-cf20d96a"><p class="cl-cf209374"><span class="cl-cf206b60"></span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="65518071-d6b8-4ff7-9a96-a01f80ec1e15"></div>
<script>
var dest = document.getElementById("65518071-d6b8-4ff7-9a96-a01f80ec1e15");
var template = document.getElementById("646b0ffc-711a-4f4c-881c-52eb506e76b0");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

The object `table_chemicals_ft` still retains information from the other columns
even when they are not displayed. These columns can still be referenced from the
saved object. The function `flextable::compose()` can then be used to add images
to the `flextable` object. The column to be used to store the image is supplied
to the `j` argument as a character and the `value` argument is supplied with a
nested combination of `as_paragraph()` and `as_image()`, the last one which will
take the path where the images are located via its `src` argument. The `width`
and `height` arguments inside of `as_image()` can be used to control the final
size of the image in the table.

``` r
table_chemicals_ft <- flextable::compose(table_chemicals_ft, 
    # Specify the column which will contain the composed information
    j = "Structure",
    value = as_paragraph(as_image(src = structure_path, width = 1.5, height = 1.5))
    )
```

Finally, we can add some light weight formatting to the table for final display.

``` r
table_chemicals_ft <- flextable::theme_vanilla(table_chemicals_ft) %>% 
      flextable::autofit()
table_chemicals_ft
```

<template id="e4fd0ed5-9a40-4090-9977-62fc37b351af"><style>
.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-cf7426f6{}.cl-cf6ba60c{font-family:'Helvetica';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-cf6ba616{font-family:'Helvetica';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-cf6bb430{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-cf6bd9c4{width:69.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9ce{width:826.8pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9d8{width:57.8pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9d9{width:128.5pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9da{width:69.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9e2{width:826.8pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9e3{width:57.8pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9e4{width:128.5pt;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9ec{width:69.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9ed{width:826.8pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9f6{width:57.8pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9f7{width:128.5pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bd9f8{width:69.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bda00{width:826.8pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bda01{width:57.8pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-cf6bda02{width:128.5pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-cf7426f6'>
<thead><tr style="overflow-wrap:break-word;"><td class="cl-cf6bda00"><p class="cl-cf6bb430"><span class="cl-cf6ba60c">SMILES</span></p></td><td class="cl-cf6bd9f8"><p class="cl-cf6bb430"><span class="cl-cf6ba60c">Chemical</span></p></td><td class="cl-cf6bda01"><p class="cl-cf6bb430"><span class="cl-cf6ba60c">Source</span></p></td><td class="cl-cf6bda02"><p class="cl-cf6bb430"><span class="cl-cf6ba60c">Structure</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-cf6bd9ce"><p class="cl-cf6bb430"><span class="cl-cf6ba616">O=C(N(C(C1=C(C(OC([H])([H])<br>[H])=C(O[H])C(=C1[H])[H])[H])([H])<br>[H])[H])C(C(C(C(C(=C(C(C([H])([H])<br>[H])(C([H])([H])[H])[H])[H])[H])<br>([H])[H])([H])[H])([H])[H])([H])<br>[H]</span></p></td><td class="cl-cf6bd9c4"><p class="cl-cf6bb430"><span class="cl-cf6ba616">Capsaicin</span></p></td><td class="cl-cf6bd9d8"><p class="cl-cf6bb430"><span class="cl-cf6ba616">Chili</span></p></td><td class="cl-cf6bd9d9"><p class="cl-cf6bb430"><img style="vertical-align:middle;width:108pt;height:108pt;" src="data:image/png;base64,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" /></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-cf6bd9e2"><p class="cl-cf6bb430"><span class="cl-cf6ba616">CCO</span></p></td><td class="cl-cf6bd9da"><p class="cl-cf6bb430"><span class="cl-cf6ba616">Ethanol</span></p></td><td class="cl-cf6bd9e3"><p class="cl-cf6bb430"><span class="cl-cf6ba616">Rum</span></p></td><td class="cl-cf6bd9e4"><p class="cl-cf6bb430"><img style="vertical-align:middle;width:108pt;height:108pt;" src="data:image/png;base64,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" /></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-cf6bd9ed"><p class="cl-cf6bb430"><span class="cl-cf6ba616">CN1C=NC2=C1C(=O)N(C(=O)N2C)C</span></p></td><td class="cl-cf6bd9ec"><p class="cl-cf6bb430"><span class="cl-cf6ba616">Caffeine</span></p></td><td class="cl-cf6bd9f6"><p class="cl-cf6bb430"><span class="cl-cf6ba616">Coffee</span></p></td><td class="cl-cf6bd9f7"><p class="cl-cf6bb430"><img style="vertical-align:middle;width:108pt;height:108pt;" src="data:image/png;base64,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" /></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="63e7c76b-d09a-46e8-863c-760775045a9f"></div>
<script>
var dest = document.getElementById("63e7c76b-d09a-46e8-863c-760775045a9f");
var template = document.getElementById("e4fd0ed5-9a40-4090-9977-62fc37b351af");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

`flextable` provides a broad range of outputs. It can save files to html,
word and powerpoint making it very versatile when creating reports and
presentations. The [R community](https://rfortherestofus.com/2019/11/how-to-make-beautiful-tables-in-r/)
has a variety of tabular formats that would also be worth exploring to handle
chemical information.

There are efforts to bring some of the functionality of RDKit into the R
ecosystem, such as the nicely looking repository [tidychem](https://github.com/nanxstats/tidychem/). The author from tidychem outlines a vision for the integration of RDKit with R. If
you are interested in learning more about RDKit, good resources are the YouTube
videos by [Jan Jensen](https://youtu.be/ERvUf_lNopo) and the website
[Is life worth living?](https://iwatobipen.wordpress.com/).

In future posts, I will discuss how we can source SMILES, SDF’s and chemical related
information using R.
