library(Distance)

source("add_df_covar_line.R")

data(minke)
ddf1 <- ds(minke, truncation=1.5, formula=~Region.Label)

# plot without the dots
plot(ddf1, showpoints=FALSE)

# add South covariate line
add_df_covar_line(ddf1, data.frame(Region.Label="South"), lty=2)

# add North
add_df_covar_line(ddf1, data.frame(Region.Label="North"), lty=3)

legend(1, 1, c("Average", "South", "North"), lty=1:3)
