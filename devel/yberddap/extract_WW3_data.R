#########################################################################
Extract wave height and wave period from WW3 for give location and date
erik.franklin@hawaii.edu
March 23, 2020
#########################################################################

#prelim
install.packages("ncdf4", dependencies = TRUE) 
install.packages("parsedate", dependencies = TRUE)
install.packages("plotdap", dependencies = TRUE)
install.packages("rerddap", dependencies = TRUE)
install.packages("sp", dependencies = TRUE)
install.packages("rerddapXtracto", dependencies = TRUE)

# if CRAN doesn't work...
#install.packages("devtools")
#devtools::install_github("rmendels/rerddapXtracto")

require("rerddap")
require("rerddapXtracto")

# data
# longitude
xpos = c(-152.269,
         -159.643,
         -176.347,
         179.1517,
         -177.778,
         -154.661,
         -156.714)
# since crosses date line
xpos.1 = 360 + xpos
xpos.1[4] = xpos.1[4] - 360

#latitude
ypos = c(18.86819,
         23.9083,
         28.7772,
         30.5724,
         29.91861,
         21.0801,
         21.58615)

# date
tpos = c('2010-08-18',
         '2010-08-26',
         '2010-09-01',
         '2010-09-03',
         '2010-09-04',
         '2010-09-19',
         '2010-09-20')

# APDRC ERDDAP server
urlBase = "http://apdrc.soest.hawaii.edu/erddap/"

# wave height
parameter = "htsgwsfc"
ww3Info = rerddap::info('hawaii_soest_98bb_253a_eb1c', url = urlBase)
ww3_swh = rxtracto(ww3Info, parameter = parameter, xcoord = xpos, ycoord = ypos, tcoord = tpos, progress_bar = TRUE)
ww3_swh

# wave period
parameter = "perpwsfc"
ww3Info = rerddap::info('hawaii_soest_98bb_253a_eb1c', url = urlBase)
ww3_per = rxtracto(ww3Info, parameter = parameter, xcoord = xpos, ycoord = ypos, tcoord = tpos, progress_bar = TRUE)
ww3_per

# table of output
output = cbind(xpos, ypos, tpos, ww3_swh$`mean htsgwsfc`, ww3_per$`mean perpwsfc`)
colnames(output) = c("lat", "lon", "date", "wave_ht_m", "wave_per_s")
write.csv(output, "wave_ht_and_per.csv")


