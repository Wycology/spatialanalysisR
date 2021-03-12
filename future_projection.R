# Making a projection into the future

# This is super simple, the beginning is to download the future climate scenario
# Here is a sample code for downloading the future scenario 

getData('CMIP5', var='tmin', res=10, rcp=85, model='AC', year=70)

# A more complete one will look like:

getData(name = 'CMIP5', # others would be worldclim, GADM, alt, SRTM(plus lon,lat)
        var = 'tmin',   # others include tmax, prec, and bio
        model = 'AC',   # others "BC", "CC", "CE", "CN", "GF", "GD", "GS", "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO"
        rcp = 45,       # others 26, 60 and 85
        year = 50,      # other 70
        res = 10,       # others 2.5 and 5
        lon = 35,       # as need be
        lat = -0.23,    # as need be
        download = T,   # in case you want to download
        path = 'D:/JUNK') # where to store downloads, default is working directory

function (var, model, rcp, year, res, lon, lat, path, download = TRUE)
  
# 'model' should be one of "AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", 
# "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO".
# 'rcp' should be one of 26, 45, 60, 80----
# 'year' should be 50 or 70

