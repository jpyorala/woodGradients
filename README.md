# woodGradients
Tools to the visualization of within-tree wood property gradients, implemented in R.

The research data and methods used in: 
Pehkonen, M., Holopainen, M., Kukko, A., Hyyppä, J. and Pyörälä, J., (Manuscript). How tree morphology reflects ring properties in Norway spruce?. 10.31219/osf.io/984tk

Installation

Clone the repository to your computer:

git clone https://github.com/jpyorala/ringGradients.git

Install the following dependencies:
nlme
rsq
reshape2
ggplot2
matrix

# Usage

Methods

The repository contains five R-scripts: 
helperFunctions.R
knotGradients.R
mainDriver.R
mixedModels.R
ringGradients.R

mainDriver.R contains the pipeline used in Pehkonen et al. (Manuscript) to read data, source and call functions, and write out results. Outline:
- Source scripts and define functions for stem and wood property models, and read in data
- Loop through target groups within the data (e.g. stands), and the target knot variables
- Fuse the knottiness data with the ring observations, and calculate features that describe the size and position of the crown with respect to the ring observations

helperFunctions.R contains three functions: organizeLogs, rings2Logs, and knots2Rings. These are the data-fusion functions used to match and combine the log-tomographic data with the ring data.
- rings2Logs ...
- organizeLogs takes as inputs the log-specific data: log dimensions and height in tree, and any wood quality- or knottiness variable. Outputs a    reorganized dataframe that can be inputted directly to the R-function lm().
- knots2Rings takes as inputs...

mixedModels.R contains two functions: one for testing all simple linear regressions with stand-level random intercepts, and another for testing different combinations of multiple linear regressions with stand-level random intercepts and slopes. The final model structures used in the manuscript are retained in the file.

woodGradients1.R and woodGradients2.R take as inputs the parameterized stem taper and knot models, and the basic information of trees within the target group: tree height (H), and diameter-at-breast-height (DBH). The functions use these models and the input data to visualize the within-tree mean gradients of the target features in the specified group (e.g. stand). The functions use ggplot2 for the visualization. woodGradients1.R works only with positional information: stem radius, and height (both in absolute and relative terms). 
woodGradients2.R also takes as inputs existing gradients, to enable downstream simulations of wood properties that are interconnected. In Pehkonen et al. (Manuscript), cambial age (CA) is expressed as a function of position (height and stem radius). Ring area (RA) is then the function of position and CA, latewood percentage (LWP) the function of position, CA and RA, and ring density (RD) the function of position, CA, RA and LWP. 
Please note that the tool is only intended for visualization purposes, and the use of predicted values in the downstream simulations are not accurate estimates, but incorporate the propagation of errors upstream. The actual prediction outcomes at the sampled locations are to be checked in mixedModels.R.


Data

An example data set is provided with the code:

stem_data.txt: contains log-to-log dimensions and knottiness variables measured at the Honkalahti sawmill (by the courtesy of Stora Enso Finland oyj), coupled with field data (1131 logs from 450 Norway spruces):
- Stand
- TreeID
- StemSection
- ForestType
- DBH
- H
- buttD
- topD
- knotIndex
- knotSizeIndex
- whorlDistance

ring_data.txt: contains ring-to-ring measurements at variable heights (stump height + log tops), coupled with field data (250 ring sample discs from 52 Norway spruces):
- Stand
- TreeID
- StemSection
- Year
- Cambial age (CA)
- Ring area (RA)
- Latewood percentage (LWP)
- Ring density (RD)
