require(sp)
require(tidyverse)
require(fs)
require(lubridate)

#########################################################################
# Response
#########################################################################

# # empirical spatial means
# spatial_avg <- suburb_df_long %>% group_by(suburb, tx, ty, longitude, latitude) %>% summarise(mu_emp = mean(area.prop))
#  
# ggplot(spatial_avg) +
#   geom_point(aes(latitude, mu_emp), colour='darkgreen', alpha=0.5) +
#   xlab("Latitude") + ylab('Tree canopy %') + theme_bw()
# 
# ggplot(spatial_avg) +
#   geom_point(aes(longitude, mu_emp), colour='darkgreen', alpha=0.5) +
#   xlab("Longitude") + ylab('Tree canopy %') + theme_bw()
# 


#########################################################################
# Temporal
#########################################################################
