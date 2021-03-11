# Created on 11th March, 2021
# Created by Wyclife Agumba Oluoch
# Task: Calculating standard deviation on sdm model outputs

# Here I am going to show how to run standard deviation on model outputs. 
# This can help reveal the regions of maximum uncertainty within the study area

# In running predict function, we generated several layers and stored as .img 
# files as specified in the finename = "blabla.img".

# Here we can also save the same image as .tif so that we can easily postanalyze it

predicted_layers <- predict(model, predictors, 
                            filename = 'predicted_bands.tif') # Note the .tif

# This is saved as a single layer with several bands (methods*replications*reps)

# To read this into R, use stack function in R

stacked_layers <- stack('predicted_bands.tif') #Reads all the bands. If you use raster,
# only one band will be read at a time.

stacked_layers # This shows that nlayers are 30 (equal to the modl outputs)

# Now we can calculate median, mean, sd et al using calc function in raster

stacked_sd <- calc(stacked_layers, fun = sd)

# Now we have our uncertainty layer to report with the precited layers

plot(stacked_sd)

# Great...
# Remember, since I saved them as tif, other beautification can be done in QGIS et al.
