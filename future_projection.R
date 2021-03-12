# Making a projection into the future

# This is super simple, the beginning is to download the future climate scenario
# Here is a sample code for downloading the future scenario 

getData('CMIP5', var='tmin', res=10, rcp=85, model='AC', year=70)

# A more complete one will look 
function (var, model, rcp, year, res, lon, lat, path, download = TRUE)
  
# 'model' should be one of "AC", "BC", "CC", "CE", "CN", "GF", "GD", "GS", 
# "HD", "HG", "HE", "IN", "IP", "MI", "MR", "MC", "MP", "MG", or "NO".
# 'rcp' should be one of 26, 45, 60, 80----
# 'year' should be 50 or 70

