# ringGradients
Tool to the visualization of within-tree wood property gradients.

Partial research data and method implementation in R used in: 
Pehkonen, M., Holopainen, M., Kukko, A., Hyyppä, J. and Pyörälä, J., 2022. How tree morphology reflects ring properties in Norway spruce?. 10.31219/osf.io/984tk

Installation

Clone the repository to your computer:

git clone https://github.com/jpyorala/ringGradients.git

Install the following dependencies:
nlme
rsq
reshape2
ggplot2

Usage

The repository contains two R-scripts: 
mainDriver.R
ringGradients.R

mainDriver.R is used to read in the ring property data, and linear mixed models for each wood property.

ringGradient.R is called from the mainDriver.R, and uses ggplot2 to visualize the model predictions as within-tree gradients. 
Input: model data, and stem data + define grouping factor (default=Stand) -> Runs through all the groups specified for the data. (better solution will be added later)

Data

An example data set is provided with the code:
ring_data.txt: contains ring-to-ring measurements sampled from 52 Norway spruces (Picea abies), at variable heights, coupled with tree-morphological features
stem_data.txt: contains log-to-log dimensions used to calculate stem taper models to each stand in ringGradient.R (better solution will be added later)
