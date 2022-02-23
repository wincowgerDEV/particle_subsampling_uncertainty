#Libraries ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv)
library(dtplyr)
library(data.table)

#Options ----
options(scipen = 999)
#template df ----
data <- tibble(
  class_num = numeric(),
  subsample_size = numeric(),
  sample_count = numeric(),
  error = numeric()
)


#functions ----
boot_mean_error_data_pub <- function(sample_subsample_size, sample_count){
  #Disallow sample count to if(sample_count < sample_classes)
  #Disallow sample round(sample_subsample_size*sample_count, 0) < 1)
  #prob = c(0.5, 0.2, rep(0.3/8, 8)) if I add this to the first sample argument, we get a less conservative estimate, the errors go down. I think it is because you are less likely to mess up a skewed distribution than a uniform one. 
  
  b = 1000 #Changing b doesn't do much, probably more accurate with larger number but this will save time while in dev mode.
  error <- numeric(length = b)
  
  set.seed(37)
  
  for(n in 1:b){
    sample_classes <- 2    #Sample classes must be greater than 1 but less than 20
    particle_categories <- 1:sample_classes

    #Simulation
    particles <- sample(particle_categories, size = sample_count, replace = T) #could add weights to this so that 1 or two of them always have a big sway. But we dont know that for sure yet.
    
    subsetparticles <- sample(particles, size = sample_subsample_size)
    
    #Difference Metric
    error[n] <- as.data.frame(table(particles)/length(particles)) %>%
      dplyr::rename(subsetparticles = particles) %>%
      left_join(as.data.frame(table(subsetparticles)/length(subsetparticles)) %>%
                  dplyr::rename(Freq2 = Freq)) %>%
      dplyr::mutate(difference = Freq - Freq2) %>%
      dplyr::mutate(difference = ifelse(is.na(difference), Freq, difference)) %>% #This sets any classes which weren't accounted for from the original group to be completely unaccounted for. 
      pull(difference) %>%
      abs() %>%
      mean() #this is the MAE or mean absolute error. Could also do mean here, might be a useful metric. Or RMSE or some other more commonly used metric.
    
  }
  
  quantile(error, probs = c(0.95))
}

microplastic_pubs <- read.csv("G:/My Drive/GrayLab/Projects/Plastics/ActiveProjects/UniWaterloo/Web of Science MP articles .csv") %>%
  filter(DOI != "")

yearly_number <- microplastic_pubs %>%
  group_by(Publication.Year) %>%
  summarize(count = n()) %>%
  filter(count > 100)

sample_sizes <- c(yearly_number$count)
subsample_size <- rep.int(100, times = length(sample_sizes))

vec_boot_mean_data_pub <- Vectorize(boot_mean_error_data_pub)

error <- vec_boot_mean_data_pub(sample_count = sample_sizes, sample_subsample_size = subsample_size)


