# ringGradients
Tools to the visualization of within-tree wood property gradients, implemented in R.

The research data and methods used in: 
Pehkonen, M., Holopainen, M., Kukko, A., Hyyppä, J. and Pyörälä, J., 2022. How tree morphology reflects ring properties in Norway spruce?. 10.31219/osf.io/984tk

Installation

Clone the repository to your computer:

git clone https://github.com/jpyorala/ringGradients.git

Install the following dependencies:
nlme
rsq
reshape2
ggplot2
matrix

Usage

The repository contains four R-scripts: 
mainDriver.R
helperFunctions.R
knotGradients.R
ringGradients.R

mainDriver.R contains the pipeline to read data, source and call functions, and write out results. Outline:
- Source scripts and define functions for stem and wood property models, and read in data
- Loop through target groups within the data (e.g. stands), and the target knot variables
- Fuse the knottiness data with the ring observations, and calculate features that describe the size and position of the crown with respect to the ring observations

helperFunctions.R contains three functions: organizeLogs, rings2Logs, and knots2Rings. These are the data-fusion functions used to match and combine the log-tomographic data with the ring data.
- rings2Logs ...
- organizeLogs takes as inputs the log-specific data: log dimensions and height in tree, and any wood quality- or knottiness variable. Outputs a    reorganized dataframe that can be inputted directly to the R-function lm().
- knots2Rings takes as inputs...

knotGradients.R takes as inputs parameterized stem taper and knot models, and a dataframe containing the basic information of trees within the target group: tree height (H), and diameter-at-breast-height (DBH). The function used the models and the input data to visualize the within-group mean gradients of the target feature as functions of height and stem diameter.

ringGradient.R is called from the mainDriver.R, and uses ggplot2 to visualize the model predictions as within-tree gradients. 
Input: model data, and stem data + define grouping factor (default=Stand) -> Runs through all the groups specified for the data. (better solution will be added later)

Data

An example data set is provided with the code:

ring_data.txt: contains ring-to-ring measurements sampled from 52 Norway spruces (Picea abies), at variable heights, coupled with tree-morphological features

stem_data.txt: contains log-to-log dimensions used to calculate stem taper models to each stand in ringGradient.R (better solution will be added later)
