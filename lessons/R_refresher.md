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
1. Download a counts file to the `data` folder by [right-clicking here](https://www.dropbox.com/s/b7n7fs453pgsy20/raw_counts_mouseKO.csv?dl=0)
1. Use `read.csv()` to read in the downloaded file and save it in the object/variable `counts`
1. What is the data structure of `counts`?
1. Create a data frame called `meta` that looks like the one shown below:  <img src="../img/refresher_meta1.png">
1. Summarize the contents of the `meta` object, how many data types are represented?
1. Check that the row names in the `meta` data frame are identical to the column names in `counts` (content and order)
1. Convert the existing `replicate` column into a factor data type
1. Create a list called `project1` with the `meta` and `counts` objects, as well as a vector with all the samplenames extracted from one of the 2 data frames.
1. Use `%>%` to selectively display only columns " " and " " from the metadata data frame
1. Display only those rows in `meta` that correspond to "..."
1. Create a tibble of the `meta` object and call it `meta_tb` (make sure you don't lose the rownames!)
1. Add the `meta_tb` object to the list `project1`
1. Plot a histogram of the counts in "sample1", give it a title that is centered and give the axes sensible names.


