#################################################
# Created on 13th May 2021                      #
# Last edited on 25th Oc 2021              #
# Author: Wyclife Agumba Oluoch                 #
# Last edited by the author on 21st October 2021#
# Projection into the future                    #
#################################################

library(raster)

getData('CMIP5', var='tmin', res=10, rcp=85, model='AC', year=70)

# A more complete one will look like:

getData(name = 'CMIP5',   # others would be worldclim, GADM, alt, SRTM(plus lon,lat)
        var = 'tmin',     # others include tmax, prec, and bio
        model = 'AC',     # others "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO"
        rcp = 45,         # others 26, 60 and 85
        year = 50,        # other 70
        res = 2.5,        # others 2.5 and 5
        lon = 35,         # as need be
        lat = -0.23,      # as need be
        download = T,     # in case you want to download
        path = 'D:/JUNK') # where to store downloads, default is current working directory