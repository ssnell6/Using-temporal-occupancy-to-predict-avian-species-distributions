# bbs_eco = rdataretriever::fetch("breed-bird-survey") 

bbs = bbs_eco$breed_bird_survey_counts
bbs$Year = as.numeric(bbs$year)
bbs$stateroute = bbs$statenum*1000 + bbs$route

weather <- bbs_eco$breed_bird_survey_weather
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))

Years <- filter(bbs, stateroute %in% RT1$stateroute & year > 2000 & year < 2016)
# bbs_eco$breed_bird_survey_counts = as.numeric(bbs_eco$breed_bird_survey_counts$Year)
Years$stateroute = Years$statenum*1000 + Years$route

#### BBS prep to replace dataset 1 #####
# Get subset of stateroutes that have been surveyed every year from 2001-2015
good_rtes =  Years %>% 
  dplyr::filter(rpid == 101) %>% 
  dplyr::select(Year, stateroute) %>%
  unique() %>%    
  dplyr::count(stateroute) %>% 
  filter(n == 15)  # have to stay at 15 to keep # of years consistent

bad_rtes =  Years %>% 
  dplyr::filter(rpid == 101) %>% 
  dplyr::select(Year, stateroute) %>%
  unique() %>%    
  dplyr::count(stateroute) %>% 
  filter(n != 15)

# Calculate occupancy for all species at subset of stateroutes above
bbs_sub1 = Years %>% 
  filter(year >= 2001 & year < 2016, stateroute %in% good_rtes$stateroute, rpid == 101) %>% 
  dplyr::select(stateroute, year, aou, speciestotal)# %>% filter(stateroute != 7008)

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



Years <- filter(bbs, stateroute %in% RT1$stateroute & year == 2016)
# bbs_eco$breed_bird_survey_counts = as.numeric(bbs_eco$breed_bird_survey_counts$Year)
Years$stateroute = Years$statenum*1000 + Years$route

# Calculate occupancy for all species at subset of stateroutes above
bbs_sub1 = Years %>% 
  filter(year == 2016, stateroute %in% good_rtes$stateroute, rpid == 101) %>% 
  dplyr::select(stateroute, year, aou, speciestotal)# %>% filter(stateroute != 7008)

# bbs_sub1 = read.csv("data/BBS/bbs_2000_2014.csv", header = TRUE)
bbs_w_aou = bbs_sub1 %>% filter(aou > 2880) %>%
  filter(aou < 3650 | aou > 3810) %>%
  filter(aou < 3900 | aou > 3910) %>%
  filter(aou < 4160 | aou > 4210) %>%
  filter(aou != 7010) # %>%
# write.csv(bbs_w_aou, "Data/bbs_2016.csv", row.names = FALSE)
