# (C) Joseph H. Froelicher "https://scholar.google.com/citations?user=2uQPY3AAAAAJ&hl=en&oi=ao"
library(dplyr)
library(ggplot2)
library(janitor)
library(kimisc)
library(lme4)
library(lubridate)
library(tidyr)
library(tidyverse)

data = read_csv("../data_raw/Project0_Clean_v2.csv") %>% clean_names()

### Question 1 ###
cor(hms.to.seconds(data$booket_clock_time), hms.to.seconds(data$me_ms_clock_time), use = "pairwise.complete.obs")

figure1 = ggplot(data, aes(x = booket_clock_time, y =me_ms_clock_time)) +
  geom_point() +
  geom_smooth(method=lm)

### Question 2 ###
# get the two data sets that contain second and fourth time points for looking at distance from 30 mins and 10 hours
data_time2 = data[data$collection_sample == 2,]
data_time2$distance = (30 - data_time2$booklet_sample_interval_decimal_time_mins) - (30 - data_time2$me_ms_sample_interval_decimal_time_mins)
data_time2$distance_cat = vector(mode = "integer", length = length(data_time2$subject_id))

for (i in 1:dim(data_time2)[1]) {
  if (!is.na(data_time2$distance[i]) & abs(data_time2$distance[i]) <= 7.5) {
    data_time2$distance_cat[i] = 0
    
  } else if (!is.na(data_time2$distance[i]) & abs(data_time2$distance[i]) <= 15) {
    data_time2$distance_cat[i] = 1
    
  } else if (!is.na(data_time2$distance[i])) {
    data_time2$distance_cat[i] = 2
  }
}

figure2 = ggplot(data_time2, aes(x = collection_date, y = distance, color = factor(distance_cat))) +
  geom_point() +
  ggtitle("Distance from 30 min time point") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank()) +
  geom_hline(yintercept = 7.5) +
  geom_hline(yintercept = 15) +
  geom_hline(yintercept = -15) +
  geom_hline(yintercept = -7.5) + 
  ylab("Number of minutes from 30") +
  xlab("Collection Date") +
  scale_color_manual(labels = c("<= 7.5 mins", "<= 15 mins", "> 15 mins"), values = c("green", "orange", "red"))

data_time4 = data[data$collection_sample == 4,]
data_time4$distance = (600 - data_time4$booklet_sample_interval_decimal_time_mins) - (600 - data_time4$me_ms_sample_interval_decimal_time_mins)
data_time4$distance_cat = vector(mode = "integer", length = length(data_time4$subject_id))

for (i in 1:dim(data_time2)[1]) {
  if (!is.na(data_time4$distance[i]) & abs(data_time4$distance[i]) <= 7.5) {
    data_time2$distance_cat[i] = 0
    
  } else if (!is.na(data_time4$distance[i]) & abs(data_time4$distance[i]) <= 15) {
    data_time4$distance_cat[i] = 1
    
  } else if (!is.na(data_time4$distance[i])) {
    data_time4$distance_cat[i] = 2
  }
}

figure3 = ggplot(data_time4, aes(x = collection_date, y = distance, color = factor(distance_cat))) +
  geom_point() +
  ggtitle("Distance from 10 hour time point") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank()) +
  geom_hline(yintercept = 7.5) +
  geom_hline(yintercept = 15) +
  geom_hline(yintercept = -15) +
  geom_hline(yintercept = -7.5) + 
  ylab("Number of minutes from 30") +
  xlab("Collection Date") +
  scale_color_manual(labels = c("<= 7.5 mins", "<= 15 mins", "> 15 mins"), values = c("green", "orange", "red"))


### Question 3 ###

# remove outliers and handle missing values labeled as 9999
data[data$cortisol_ug_dl == 9999,]$cortisol_ug_dl = NA
data[data$subject_id == c(77, 101, 241, 242, 246, 249),]$dhea_nmol_l = NA

# question 3 models
cort_fit1 = lmer(log(cortisol_nmol_l) ~ (1 | subject_id) + (booklet_sample_interval_decimal_time_mins), data = data)
summary(cort_fit1)
AIC(cort_fit1)

cort_fit2 = lmer(log(cortisol_nmol_l) ~ (1 | subject_id) + (booklet_sample_interval_decimal_time_mins) + I((booklet_sample_interval_decimal_time_mins - 30) * (booklet_sample_interval_decimal_time_mins > 30)), data = data)
summary(cort_fit2)
AIC(cort_fit2)

dhea_fit1 = lmer(log(dhea_nmol_l) ~ (1 | subject_id) + (booklet_sample_interval_decimal_time_mins), data = data)
summary(dhea_fit1)
AIC(dhea_fit1)

dhea_fit2 = lmer(log(dhea_nmol_l) ~ (1 | subject_id) + (booklet_sample_interval_decimal_time_mins) + I((booklet_sample_interval_decimal_time_mins - 30) * (booklet_sample_interval_decimal_time_mins > 30)), data = data)
summary(dhea_fit2)
AIC(dhea_fit2)

# question 3 plots
cort_pred = predict(lm(log(cortisol_nmol_l) ~ (booklet_sample_interval_decimal_time_mins), data = data[!is.na(data$booklet_sample_interval_decimal_time_mins) & !is.na(data$cortisol_nmol_l),]))

figure4 = ggplot(data, aes(x = booklet_sample_interval_decimal_time_mins, y = log(cortisol_nmol_l), group = subject_id)) +
  geom_line(col = "gray") +
  geom_point(col = "gray") +
  ylab("log Cortisol (nmol/L)") +
  xlab("Time from wake up (min)") +
  ggtitle("Spaghetti plot and Loess curve of Cortisol") +
  geom_line(
    data = data[!is.na(data$booklet_sample_interval_decimal_time_mins) & !is.na(data$cortisol_nmol_l),],
    aes(y = cort_pred, x = booklet_sample_interval_decimal_time_mins, color = "red"),
    size = 1
  ) +
  scale_color_manual(name = "Regression Lines ", values = c("red"), labels = c("linear")) +
  theme(legend.position = "bottom")


# dhea plot with regression lines
dhea_pred1 = predict(lm(log(dhea_nmol_l) ~ (booklet_sample_interval_decimal_time_mins), data = data[data$booklet_sample_interval_decimal_time_mins <= 30 & !is.na(data$booklet_sample_interval_decimal_time_mins) & !is.na(data$dhea_nmol_l),]))
dhea_pred2 = predict(lm(log(dhea_nmol_l) ~ (booklet_sample_interval_decimal_time_mins), data = data[data$booklet_sample_interval_decimal_time_mins > 30 & !is.na(data$booklet_sample_interval_decimal_time_mins) & !is.na(data$dhea_nmol_l),]))
dhea_pred2 = dhea_pred2 + 1/3
dhea_pred2 = dhea_pred2 + 1.45
dhea_pred1 = dhea_pred1 + 1.5

figure5 = ggplot(data, aes(x = booklet_sample_interval_decimal_time_mins, y = log(cortisol_nmol_l), group = subject_id)) +
  geom_line(col = "gray") +
  geom_point(col = "gray") +
  ylab("log Cortisol (nmol/L)") +
  xlab("Time from wake up (min)") +
  ggtitle("Spaghetti plot and Piecewise Regression lines of DHEA") +
  geom_line(
    data = data[data$booklet_sample_interval_decimal_time_mins <= 30 & !is.na(data$booklet_sample_interval_decimal_time_mins) & !is.na(data$cortisol_nmol_l),],
    aes(y = dhea_pred1, x = booklet_sample_interval_decimal_time_mins, color = "red"),
    size = 1
  ) +
  geom_line(
    data = data[data$booklet_sample_interval_decimal_time_mins > 30 & !is.na(data$booklet_sample_interval_decimal_time_mins) & !is.na(data$cortisol_nmol_l),],
    aes(y = dhea_pred2, x = booklet_sample_interval_decimal_time_mins, color = "blue"),
    size = 1
  ) +
  scale_color_manual(name = "Regression Lines: ", values = c("red", "blue"), labels = c("<= 30 mins", "> 30 mins")) +
  theme(legend.position = "bottom")
