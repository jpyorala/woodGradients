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
mainDriver.R
woodGradient.R
gradientFigure.R
mixedModels.R
helperFunctions.R

mainDriver.R contains the pipeline used in Pehkonen et al. (Manuscript) to read data, source and call functions, and write out results. Outline:
- Source scripts and define functions for stem and wood property models, and read in data
- Loop through target groups within the data (e.g. stands), and the target knot variables
- Fuse the knottiness data with the ring observations, and calculate features that describe the size and position of the crown with respect to the ring observations
- Call simple and multiple mixed models and save results
- Visualize the ring property gradients

woodGradient.R contains one function with the same name. The function takes as inputs the parameterized stem taper and knot models, and the basic information of trees within the target group: tree height (H), and diameter-at-breast-height (DBH). The function uses these models and the input data to extrapolate the targeted wood properties within stems, and outputs a 2-D matrix (stem diameter x height) containing the values at specific locations. Default resolution is at 1 decimeter for the height, and at 1 mm for the diameter.
Optionally, the function takes as inputs pre-calculated gradient matrices, to enable downstream simulations of wood properties that are interconnected. In Pehkonen et al. (Manuscript), cambial age (CA) is expressed as a function of position (height and stem radius). Ring area (RA) is then the function of position and CA, latewood percentage (LWP) the function of position, CA and RA, and ring density (RD) the function of position, CA, RA and LWP. See mainDriver.R, lines ... Please note that when the function is used in this way, the outcome can only be considered as a visualization, or a demonstration: the use of predicted values as inputs in the downstream simulations will propagate all errors upstream. The actual prediction outcomes and their statistical properties at the sampled locations are to be checked in mixedModels.R.

gradientFigure.R visualizes the within-tree mean gradients of the target features in the specified group (e.g. stand), using ggplot2. Takes as inputs a matrix containing the values to be plotted, variable name, color scale, and optional parameters regarding the viewing and saving of the plotted figure.

mixedModels.R contains two functions: one for testing all simple linear regressions with stand-level random intercepts, and another for testing different combinations of multiple linear regressions with stand-level random intercepts and slopes. The final model structures used in the manuscript are retained in the file.

helperFunctions.R contains three functions: organizeLogs, rings2Logs, and knots2Rings. These are the data-fusion functions used to match and combine the log-tomographic data with the ring data in the example study.
- rings2Logs matches the rings in the ring data to corresponding sawmill data, and field measurements.
- organizeLogs takes as inputs the log-specific data: log dimensions and height in tree, and any wood quality- or knottiness variable. Outputs a    reorganized dataframe that can be inputted directly to the R-function lm().
- knots2Rings takes as inputs the matched ring and stem data, and the calculated knot gradients: calculates tree morphological features to each ring observation.


Data

An example data set is provided with the code:

stem_data.txt: contains log-to-log dimensions and knottiness variables measured at the Honkalahti sawmill (by the courtesy of Stora Enso Finland oyj), coupled with field data (1131 logs from 450 Norway spruces):
- Stand
- TreeID
- StemSection: by height, 1=butt log etc...
- ForestType: according to the Finnish forest types. MT=Myrtillus, OMT=Oxalis-Myrtillus
- DBH: Diameter at breast height (mm)
- H: Tree height (dm)
- buttD: Diameter of the log's butt-end (mm)
- topD: Diameter of the log's top-end (mm)
- knotIndex: ratio of knot volume to log volume
- knotSizeIndex: mean whorlm volume (cm3)
- whorlDistance: mean distance between adjacent whorls within log (mm)

ring_data.txt: contains ring-to-ring measurements at variable heights (stump height + log tops), coupled with field data (250 ring sample discs from 52 Norway spruces):
- Stand
- TreeID
- StemSection
- Year
- Cambial age (CA)
- Ring area (RA)
- Latewood percentage (LWP)
- Ring density (RD)
