bbs_eco = rdataretriever::fetch("breed-bird-survey") 

Years = (bbs_eco$counts$Year)
bbs = bbs_eco$breed_bird_survey_counts
bbs$Year = as.numeric(bbs$year)
bbs$stateroute = bbs$statenum*1000 + bbs$route

Years = (bbs_eco$breed_bird_survey_counts)
# bbs_eco$breed_bird_survey_counts = as.numeric(bbs_eco$breed_bird_survey_counts$Year)
Years$stateroute = Years$statenum*1000 + Years$route

#### BBS prep to replace dataset 1 #####
# Get subset of stateroutes that have been surveyed every year from 2001-2015
good_rtes =  Years %>% 
  dplyr::filter(year > 2000, year < 2016) %>% 
  dplyr::select(year, stateroute) %>%
  unique() %>%    
  dplyr::count(stateroute) %>% 
  filter(n == 15)  # have to stay at 15 to keep # of years consistent

# Calculate occupancy for all species at subset of stateroutes above
bbs_sub1 = bbs %>% 
  filter(year >= 2001 & year < 2016, stateroute %in% good_rtes$stateroute) %>% 
  dplyr::select(stateroute, Year, aou, speciestotal)# %>% filter(stateroute != 7008)

# bbs_sub1 = read.csv("data/BBS/bbs_2000_2014.csv", header = TRUE)
bbs_w_aou = bbs_sub1 %>% filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010) # %>%
  # filter(sporder != "Accipitriformes", 
  #       sporder != "Falconiformes", 
  #       sporder != "Anseriformes",
  #       sporder != "Cathartiformes")
# write.csv(bbs_w_aou, "Data/bbs_2001_2015.csv", row.names = FALSE)
# write.csv(bbs_w_aou, "Data/bbs_2015_on.csv", row.names = FALSE)
