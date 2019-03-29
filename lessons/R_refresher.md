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

## Setting up

Let’s create a new project directory for our "R refresher lesson:
  
  - Create a new project called `R_refresher`
  - Create a new R script called `reviewing_R.R`
  - Create the following folders in the project directory - `data`, `figures`
  - Download a counts file to the `data` folder by [right-clicking here](https://github.com/hbctraining/DGE_workshop_salmon/blob/master/data/raw_counts_mouseKO.csv?raw=true)

Now that we have our directory structure set-up, let's load our libraries and read in our data:

  - Load the tidyverse library
  - Use `read.csv()` to read in the downloaded file and save it in the object/variable `counts`
  - What is the data structure of `counts`?

## Creating vectors/factors and dataframes

We are performing RNA-Seq on cancer samples being with genotypes of p53 WT and KO. You have 8 samples total, with 4 replicates per genotype. Write the R code you would use to construct your metadata table as described below.  
     - Create the vectors/factors for each column (Hint: you can type out each vector/factor, or if you want the process go faster try exploring the `rep()` function).
     - Put them together into a dataframe called `meta`.
     - Use the `rownames()` function to assign row names to the dataframe (Hint: you can type out the row names as a vector, or if you want the process go faster try exploring the `paste0()` function).
     
Your finished metadata table should have information for the variables `sex`, `stage`, `genotype`, and `myc` levels: 

     | |sex	| stage	| genotype	| myc |
     |:--:|:--: | :--:	| :------:	| :--: |
     |KO1 |	M	|I	|KO	|23|
     |KO2|	F	|II	|KO	|4|
     |KO3	|M	|II	|KO	|45|
     |KO4	|F	|I	|KO	|90|
     |WT1|	M	|II	|WT	|34|
     |WT2|	F|	II|	WT|	35|
     |WT3|	M|	I|	WT|	9|
     |WT4|	F|	II|	WT|	10|

## Subsetting vectors/factors and dataframes

2. Using the `meta` data frame from question #1, write out the R code you would use to perform the following operations (questions **DO NOT** build upon each other):

     - return only the `genotype` and `sex` columns using `[]`:
     - return the `genotype` values for samples 5, 7, 9, and 10 using `[]`:
     - use `filter()` to return all data for those samples receiving treatment `P`:
     - use `filter()`/`select()`to return only the `stage` and `genotype` columns for those samples with `myc` > 50:
     - add a column called `pre_treatment` to the beginning of the dataframe with the values T, F, T, F, T, F, T, F 
     - change the names of the columns to: "A", "B", "C", "D":

     
## Refresher exercises!



1. Create a data frame called `meta` that looks like the one shown below:  <img src="../img/refresher_meta1.png">
1. Summarize the contents of the `meta` object, how many data types are represented?
1. Check that the row names in the `meta` data frame are identical to the column names in `counts` (content and order)
1. Convert the existing `replicate` column into a factor data type
1. Create a list called `project1` with the `meta` and `counts` objects, as well as a new vector with all the sample names extracted from one of the 2 data frames.
1. Use `%>%` to get the average of the age_in_days column from `meta`
1. Display only the metadata for those samples which have an age of 17 days or older
1. Using `%>%` create a tibble of the `meta` object and call it `meta_tb` (make sure you don't lose the rownames!)
1. Add the `meta_tb` object to the list `project1`
1. Plot a boxplot of the mean_expression of the KO and WT samples using `theme_minimal()` and give the plot new axes names and a centered title.


