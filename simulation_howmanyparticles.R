#Libraries ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

#Options ----
options(scipen = 99)
population_sizes <- 10^(1:10)
groups <- 1:5
max_error <- c(10^(c(-5:-1)), 0.05)
subsample_size <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
count <- c(10, 100, 10^3, 10^4, 10^5, 10^6) #Probably decrease the size here. 

#template df ----
data <- tibble(
    class_num = numeric(),
    subsample_size = numeric(),
    sample_count = numeric(),
    error = numeric()
)

#functions ----

finite_corrected_eq <- function(pop_size, uncorrected){
  round(pop_size * uncorrected / (uncorrected + pop_size - 1), 0)
} 

uncorrected_eq <- function(confidence = 0.95, groups = 1, expected_p = 0.5, max_error = 0.05){
  critical_val = qnorm(p=(1-confidence^(1/groups))/2) # The group root is a new thing here. 
  round(critical_val^2*expected_p*(1-expected_p)/max_error^2, 0)
}

boot_mean_error_class_num <- function(sample_subsample_size, sample_count, classes){
    #Disallow sample count to if(sample_count < sample_classes)
    #Disallow sample round(sample_subsample_size*sample_count, 0) < 1)
    #prob = c(0.5, 0.2, rep(0.3/8, 8)) if I add this to the first sample argument, we get a less conservative estimate, the errors go down. I think it is because you are less likely to mess up a skewed distribution than a uniform one. 
    
    b = 100 #Changing b doesn't do much, probably more accurate with larger number but this will save time while in dev mode.
    error <- numeric(length = b)
    
    for(n in 1:b){
        sample_classes <- classes    #Sample classes must be greater than 1 but less than 20
        particle_categories <- 1:sample_classes
        
        #Add gausian distribution to classes
        #values <- rnorm(n = length(particle_categories), mean = 10, sd = 1)
        #Or can swap for a poison distribution
        values <- rpois(n = length(particle_categories), lambda = 1) + 1
        weights <- values/sum(values)
        
        #Simulation
        particles <- sample(particle_categories, size = sample_count, prob = weights, replace = T) #could add weights to this so that 1 or two of them always have a big sway. But we dont know that for sure yet.
        
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
            quantile(., c(0.95)) #this is the MAE or mean absolute error. Could also do mean here, might be a useful metric. Or RMSE or some other more commonly used metric.
        
    }
    
    error
}

boot_mean <- function(x, b = 100){
  
  value <- numeric(length = b)
  
  for(n in 1:b){
    value[n] <- mean(sample(x, b, replace = T), na.rm = T)
  }
  
  value
}

boot_mean_error <- function(sample_subsample_size, sample_count){
    #Disallow sample count to if(sample_count < sample_classes)
    #Disallow sample round(sample_subsample_size*sample_count, 0) < 1)
    #prob = c(0.5, 0.2, rep(0.3/8, 8)) if I add this to the first sample argument, we get a less conservative estimate, the errors go down. I think it is because you are less likely to mess up a skewed distribution than a uniform one. 
    
    b = 100 #Changing b doesn't do much, probably more accurate with larger number but this will save time while in dev mode.
    error <- numeric(length = b)
    
    for(n in 1:b){
        sample_classes <- sample(2:10, 1)    #Sample classes must be greater than 1 but less than 20
        particle_categories <- 1:sample_classes
        
        #Add gausian distribution to classes
        #values <- rnorm(n = length(particle_categories), mean = 10, sd = 1)
        #Or can swap for a poison distribution
        values <- rpois(n = length(particle_categories), lambda = 1) + 1
        weights <- values/sum(values)
        
        #Simulation
        particles <- sample(particle_categories, size = sample_count, prob = weights, replace = T) #could add weights to this so that 1 or two of them always have a big sway. But we dont know that for sure yet.
        
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
                    quantile(., c(0.95)) #this is the MAE or mean absolute error. Could also do mean here, might be a useful metric. Or RMSE or some other more commonly used metric.
                
    }
    
    error
}


#test single scenario----
errors <- boot_mean_error(sample_count = 10, sample_subsample_size = 1) #Just say we will set this to 20, that is really what can be reliably characterized, can do a quick lit search to figure out how people are reporting it.    

hist(errors)
quantile(errors, c(0.025, 0.5, 0.975))
mean(errors)

#Simulate all scenarios----
test_df <- expand.grid(subsample_size, count) %>%
    filter(Var1 * Var2 >= 1) %>% 
    mutate(num_particles = Var1*Var2) %>%
        rename(sample_count = Var2 , subsample_proportion = Var1) %>%
    as_tibble()

#Just polymer types
for(n in 1:nrow(test_df)){
    set.seed(37)
    test_df[n, "median_error"] <- quantile(
        boot_mean_error(sample_count = unlist(test_df[n,"sample_count"]), 
                        sample_subsample_size = unlist(test_df[n,"num_particles"])),
        c(0.95)
    )
}

#Polymers, colors, morphologies, and sizes
for(n in 1:nrow(test_df)){
    set.seed(38)
    print(n)
    test_df[n, "median_error_multiple"] <- quantile(
        c(
        quantile(
        boot_mean_error(sample_count = unlist(test_df[n,"sample_count"]), 
                        sample_subsample_size = unlist(test_df[n,"num_particles"])),
        c(0.95)
    ), 
    quantile(
        boot_mean_error(sample_count = unlist(test_df[n,"sample_count"]), 
                        sample_subsample_size = unlist(test_df[n,"num_particles"])),
        c(0.95)
    ), 
    quantile(
        boot_mean_error(sample_count = unlist(test_df[n,"sample_count"]), 
                        sample_subsample_size = unlist(test_df[n,"num_particles"])),
        c(0.95)
        ),
    quantile(
        boot_mean_error(sample_count = unlist(test_df[n,"sample_count"]), 
                        sample_subsample_size = unlist(test_df[n,"num_particles"])),
        c(0.95)
    )
    ),
    c(0.95)
    )
}

# Real Data Example ----

lake_data <- fread("datasets_SubsamplingDatasets/2_MPs_features.csv")

env_data_test <- function(lake_data, subsample_size) {
 
  full_results_shape <- lake_data |>
    group_by(sample_lake, shape) |>
    summarise(count = n()) |>
    ungroup() |>
    group_by(sample_lake) |>
    mutate(proportion = count/sum(count)) |>
    ungroup()
  
  full_results_color <- lake_data |>
    group_by(sample_lake, color) |>
    summarise(count = n()) |>
    ungroup() |>
    group_by(sample_lake) |>
    mutate(proportion = count/sum(count)) |>
    ungroup()
  
  results_subsampled <-  lake_data |>
    group_by(sample_lake) |>
    mutate(sample_particles = n()) |>
    filter(sample_particles > subsample_size) |>
    sample_n(size = subsample_size, replace = F) 
  
  results_subsampled_shape <-  results_subsampled |>
    group_by(sample_lake, shape) |>
    summarise(count = n()) |>
    ungroup() |>
    group_by(sample_lake) |>
    mutate(proportion = count/sum(count)) |>
    ungroup() |>
    left_join(full_results_shape, by = c("sample_lake", "shape")) |>
    mutate(error = abs(proportion.y - proportion.x)) |> 
    group_by(sample_lake) |>
    summarise(max_error = quantile(error, 0.95),
              sample_size = sum(count.y), 
              subsample_size = subsample_size, 
              max_prop = max(proportion.y)
    ) |>
    ungroup() 
  
  results_subsampled_color <-  results_subsampled |>
    group_by(sample_lake, color) |>
    summarise(count = n()) |>
    ungroup() |>
    group_by(sample_lake) |>
    mutate(proportion = count/sum(count)) |>
    ungroup() |>
    left_join(full_results_color, by = c("sample_lake", "color")) |>
    mutate(error = abs(proportion.y - proportion.x)) |> 
    group_by(sample_lake) |>
    summarise(max_error = quantile(error, 0.95),
              sample_size = sum(count.y), 
              subsample_size = subsample_size, 
              max_prop = max(proportion.y)
    ) |>
    ungroup() 
  
  joined_results <- bind_rows(results_subsampled_shape, results_subsampled_color) |>
    group_by(sample_lake) |>
    filter(max_error == max(max_error)) |>
    mutate(mathematic_sub_count = finite_corrected_eq(pop_size = sample_size, uncorrected = uncorrected_eq(max_error = max_error, 
                                                                                                           #expected_p = max_prop, 
                                                                                                           groups = 2)))


 bootsy <- boot_mean(joined_results$mathematic_sub_count) |>
   quantile(probs = c(0.025, 0.975)) 
 
 data.table(subsample_size = subsample_size, 
            mean = mean(joined_results$mathematic_sub_count, na.rm = T), 
            min = bootsy[1], 
            max = bootsy[2])
}

set.seed(100)
enviro_data <- lapply(seq(10, 200, by = 10), function(x){env_data_test(lake_data = lake_data, subsample_size = x)}) |>
  rbindlist()

# Equations ----

test_ungrouped <- expand.grid(population_sizes = population_sizes, max_error = max_error) %>%
  mutate(sample_size = finite_corrected_eq(pop_size = population_sizes, uncorrected = uncorrected_eq(max_error = max_error)))

test_grouped <- expand.grid(population_sizes = population_sizes, max_error = max_error) %>%
  mutate(sample_size = finite_corrected_eq(pop_size = population_sizes, uncorrected = uncorrected_eq(max_error = max_error, groups = 4)))


#Plots ----
ggplot() +
  geom_point(data = test_ungrouped %>% arrange(desc(max_error)), aes(x = population_sizes, y = sample_size, color = factor(max_error, levels = c(10^(c(-5:-1)), 0.05))), size = 10) +
  geom_point(data = test_grouped %>% arrange(desc(max_error)), aes(x = population_sizes, y = sample_size, color = factor(max_error, levels = c(10^(c(-5:-1)), 0.05))), size = 10, shape = 2) +
  scale_x_log10(breaks = unique(population_sizes)) +
  scale_y_log10(breaks = unique(population_sizes)) + 
  scale_color_viridis_d() +
  theme_dark(base_size = 10) +
  labs(color = "Error", x = "Sample Size", y = "Subsample Size") +
  coord_equal()

test_nopop_nogroup <- expand.grid(population_sizes = population_sizes, max_error = max_error) %>%
  mutate(sample_size = uncorrected_eq(max_error = max_error))

test_nopop <- expand.grid(population_sizes = population_sizes, max_error = max_error) %>%
  mutate(sample_size = uncorrected_eq(max_error = max_error, groups = 4))

ggplot() +
  geom_line(data = test_nopop %>% arrange(desc(max_error)), aes(x = max_error, y = sample_size), size = 2) +
  geom_line(data = test_nopop_nogroup %>% arrange(desc(max_error)), aes(x = max_error, y = sample_size), size = 2, color = "yellow") +
  scale_x_log10(breaks = unique(max_error)) +
  scale_y_log10(breaks = unique(population_sizes)) + 
  theme_dark(base_size = 20) +
  labs(x = "Error", y = "Sample Size") +
  coord_equal()


simulation_validation <- test_df %>%
  mutate(mathematic_sub_count = finite_corrected_eq(pop_size = sample_count, uncorrected = uncorrected_eq(max_error = median_error))) %>%
  mutate(mathematic_sub_count_grouped = finite_corrected_eq(pop_size = sample_count, uncorrected = uncorrected_eq(max_error = median_error_multiple, groups = 4)))

ggplot(simulation_validation) +
  geom_point(aes(x = num_particles, y = mathematic_sub_count), size = 4, alpha = 0.5) + 
  geom_point(aes(x = num_particles, y = mathematic_sub_count_grouped), size = 4, alpha = 0.5, color = "red") + 
  scale_x_log10(limits = c(10,150000)) + 
  scale_y_log10(limits = c(10,150000)) +
  geom_abline(slope=1, intercept = 0) +
  theme_dark(base_size = 20) +
  coord_equal() +
  labs(x = "Simulated Sample Counts", y = "Calculated Sample Counts")

#final recommendations for sample size and subsample size. ----
uncorrected_eq()
uncorrected_eq(max_error = 0.1)

#final recommendation for all group analysis at once.
uncorrected_eq(groups = 4)

#Number of papers to review 
finite_corrected_eq(pop_size = 1000, uncorrected = uncorrected_eq())

0.95^(1/3)

#Checking correspondence with simulations. 
finite_corrected_eq(pop_size = 100, uncorrected = uncorrected_eq(max_error = 0.3))
finite_corrected_eq(pop_size = 1000, uncorrected = uncorrected_eq(max_error = 0.1))
finite_corrected_eq(pop_size = 10000, uncorrected = uncorrected_eq(max_error = 0.03))
finite_corrected_eq(pop_size = 100000, uncorrected = uncorrected_eq(max_error = 0.01))
finite_corrected_eq(pop_size = 1000000, uncorrected = uncorrected_eq(max_error = 0.003))

#Checking correspondence with simulations. 
finite_corrected_eq(pop_size = 100, uncorrected = uncorrected_eq(max_error = 0.3, groups = 50))
finite_corrected_eq(pop_size = 1000, uncorrected = uncorrected_eq(max_error = 0.1, groups = 50))
finite_corrected_eq(pop_size = 10000, uncorrected = uncorrected_eq(max_error = 0.03, groups = 50))
finite_corrected_eq(pop_size = 100000, uncorrected = uncorrected_eq(max_error = 0.01, groups = 50))
finite_corrected_eq(pop_size = 1000000, uncorrected = uncorrected_eq(max_error = 0.003, groups = 50))

