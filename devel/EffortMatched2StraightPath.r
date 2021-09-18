#Acoustic Effort Matched to Straight Segments of Trackline (using the straighPath output)
#Yvonne Barkley
# April 2, 2020


### Files needed:
straightPath_1303.csv
Effort_1303.csv
SWlocations_1303.csv
###


# Determine when and where the acoustics team was ON effort along straight segments.

library(lubridate) 
library(tidyverse)
library(RcppRoll)
library(data.table)


#load straight path data and spermie location data
mydir = 'C:\\Users\\yvers\\Documents\\CHP 3\\SpermWhales\\data\\'

survey = 1303
str8segs = read.csv(paste0(mydir, 'straightPath_', survey, '.csv'))
str8SW = read.csv(paste0(mydir, 'SWlocations_', survey, '.csv'))


#load effort data
eff <- read.csv(paste0(mydir, 'Effort_', survey, '.csv'))
colnames(eff) <- c("UTC", "PCTime", "cruise_number", "ac_observer", "ac_effort", "vis_effort")

eff$UTC <- mdy_hm(as.character(eff$UTC))
eff$PCTime <- mdy_hm(as.character(eff$PCTime))


#####
# Match up the effort to the straight path times using data.table functions
df1 <- as.data.frame(as.POSIXct(str8segs$UTC ))   #gps dataframe after using straightPath
colnames(df1) <- "UTC"

df2 <- data.frame(as.POSIXct(eff$PCTime), eff$ac_effort)   #effort dataframe
colnames(df2) <- c("UTC", "aeffort")

setDT(df1)
setDT(df2)

df1[, date := UTC] #creates duplicate of UTC
setkey(df1, date) #set column to perform the join on
setkey(df2, UTC)

join = df2[df1, roll = Inf]  

#make df adding the acoustic effort column to the straightPath output
str8eff <- data.frame(str8segs, join$aeffort) %>% drop_na(join.aeffort)


#Plot results of effort with sperm whales colored by the acoustic effort

(plotSperm <- ggplot(str8eff, aes(x=Longitude, y=Latitude)) +
    geom_point(aes(group=timeGroup, col=join.aeffort), size = 1) +
    scale_color_manual(values = c('red', 'darkgreen')) +
    
    geom_point(data=str8SW, mapping = aes(lon, lat), shape = 15))

