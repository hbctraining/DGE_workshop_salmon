---
title: "R refresher"
author: "Meeta Mistry, Radhika Khetani, Mary Piper"
---

Approximate time: 45 minutes

## Learning Objectives

* Describe the various data types and data structures (including tibbles) used by R
* Use functions in R and describe how to get help with arguments
* Describe how to install and use packages in R
* Use the pipe (%>%) from the dplyr package
* Describe the syntax used by ggplot2 for making plots

## Refresher exercises!

1. Create a new project called `R_refresher`
1. Create a new R script called `reviewing_R.R`
1. Load the tidyverse and ggplot2 libraries
1. Create the following folders in the project directory - `data`, `figures`
1. Download a counts file to the `data` folder by [right-clicking here](https://github.com/hbctraining/DGE_workshop_salmon/blob/master/data/raw_counts_mouseKO.csv?raw=true)
1. Use `read.csv()` to read in the downloaded file and save it in the object/variable `counts`
1. What is the data structure of `counts`?
1. Create a data frame called `meta` that looks like the one shown below:  <img src="../img/refresher_meta1.png">
1. Summarize the contents of the `meta` object, how many data types are represented?
1. Check that the row names in the `meta` data frame are identical to the column names in `counts` (content and order)
1. Convert the existing `replicate` column into a factor data type
1. Create a list called `project1` with the `meta` and `counts` objects, as well as a new vector with all the sample names extracted from one of the 2 data frames.
1. Use `%>%` to get the average of the age_in_days column from `meta`
1. Display only the metadata for those samples which have an age of 17 days or older
1. Using `%>%` create a tibble of the `meta` object and call it `meta_tb` (make sure you don't lose the rownames!)
1. Add the `meta_tb` object to the list `project1`
1. Plot a boxplot of the mean_expression of the KO and WT samples, give it a title that is centered and give the axes new names.


