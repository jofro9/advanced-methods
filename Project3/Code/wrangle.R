library(janitor)

data = read.csv("../data_raw/frmgham2.csv") %>% clean_names()

# ten year risk! look at data before 10 years
# hoping to do an analysis on people who have never had a stroke
# check how many people have had a stroke at baseline
# we are studying to the first stroke ("survival")
# time till the first stroke is the survival
