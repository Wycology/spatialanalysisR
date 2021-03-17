# Created on 11th March 2021
# Last edited 17th March 2021
# Created by Wyclife Agumba Oluoch
# Task: Calculating standard deviation on sdm model outputs
# Relevance: This is super helpful with appreciating uncertainty around model
# outputs...like where do the models disagree most and where to they agree more

# This script is not standalone it continues from prediction stage of sdm

# In running predict function, we generated several layers and stored as .img 
# files as specified in the finename = "blabla.img".

# Here we can also save the same image as .tif so that we can easily postanalyze 
# them in any other QGIS programme. 

predicted_layers <- predict(model, predictors, 
                            filename = 'predicted_bands.tif') # Note the .tif

# This is saved as a single layer with several bands. Each band represent model 
# id (methods*replications*reps). So number 1 can be SVM under bootsrapping for 
# first rep.

# To read this into R, use stack function in R

stacked_layers <- stack('predicted_bands.tif') #Reads all the bands. If you use 
# raster, only one band will be read at a time.

stacked_layers # This shows that nlayers are 30 (equal to my model outputs)

# Now we can calculate median, mean, sd et al using calc function in raster

stacked_sd <- calc(stacked_layers, fun = sd)

# Now we have our uncertainty layer to report with the precited layers

plot(stacked_sd)

# Great...
# Remember, since I saved them as tif, other beautification can be done in QGIS 
# et al.
