#Libraries ----
library(dplyr)
library(ggplot2)
library(tidyr)
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

sample_subsample_size = 100
sample_count = 1000

#functions ----

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
            dplyr::mutate(difference = (Freq2 - Freq)/Freq * 100) %>%
            dplyr::mutate(difference = ifelse(is.na(difference), 100, difference)) %>% #This sets any classes which weren't accounted for from the original group to be completely unaccounted for. 
            pull(difference) %>%
            abs() %>%
            quantile(., c(0.95)) #this is the MAE or mean absolute error. Could also do mean here, might be a useful metric. Or RMSE or some other more commonly used metric.
        
    }
    
    error
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

#boot_mean_error_vector <- Vectorize(boot_mean_error)


#test single scenario----
errors <- boot_mean_error(sample_count = 10, sample_subsample_size = 1) #Just say we will set this to 20, that is really what can be reliably characterized, can do a quick lit search to figure out how people are reporting it.    

hist(errors)
quantile(errors, c(0.025, 0.5, 0.975))
mean(errors)

#Simulate all scenarios----

subsample_size <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
count <- c(10, 100, 10^3, 10^4, 10^5, 10^6) #Probably decrease the size here. 

test_df <- expand.grid(subsample_size, count) %>%
    filter(Var1 * Var2 >= 32) %>% #change to 1 for earlier figures. 
    #mutate(error = quantile(
    #                    boot_mean_error_vector(
    #                        sample_count = Var2, 
    #                        sample_subsample_size = Var1),
    #                c(0.95)
    #                )
    #) %>%
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

#Just polymer types with class
for(n in 1:nrow(test_df)){
    set.seed(37)
    test_df[n, "median_error_2class"] <- quantile(
        boot_mean_error_class_num(sample_count = unlist(test_df[n,"sample_count"]), 
                        sample_subsample_size = unlist(test_df[n,"num_particles"]), 
                        classes = 2),
        c(0.95)
    )
}

for(n in 1:nrow(test_df)){
    set.seed(37)
    test_df[n, "median_error_4class"] <- quantile(
        boot_mean_error_class_num(sample_count = unlist(test_df[n,"sample_count"]), 
                                  sample_subsample_size = unlist(test_df[n,"num_particles"]), 
                                  classes = 4),
        c(0.95)
    )
}

for(n in 1:nrow(test_df)){
    set.seed(37)
    test_df[n, "median_error_8class"] <- quantile(
        boot_mean_error_class_num(sample_count = unlist(test_df[n,"sample_count"]), 
                                  sample_subsample_size = unlist(test_df[n,"num_particles"]), 
                                  classes = 8),
        c(0.95)
    )
}

for(n in 1:nrow(test_df)){
    set.seed(37)
    test_df[n, "median_error_16class"] <- quantile(
        boot_mean_error_class_num(sample_count = unlist(test_df[n,"sample_count"]), 
                                  sample_subsample_size = unlist(test_df[n,"num_particles"]), 
                                  classes = 16),
        c(0.95)
    )
}

for(n in 1:nrow(test_df)){
    set.seed(37)
    test_df[n, "median_error_32class"] <- quantile(
        boot_mean_error_class_num(sample_count = unlist(test_df[n,"sample_count"]), 
                                  sample_subsample_size = unlist(test_df[n,"num_particles"]), 
                                  classes = 32),
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

#plot of tiles ----

#cut <- cut(test_df$median_error, c( 0.0001, 0.001, 0.01, 0.1, 1))
ggplot(test_df, aes(x = sample_count, y = num_particles)) +
    geom_tile(aes(fill = log10(median_error)))+ 
    geom_text(aes(label = round(median_error, 3)))+ 
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                  labels = c(1, 10, 100, 1000, 10000, 100000, 1000000))+ 
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                  labels = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
    scale_fill_viridis_c() + 
    coord_equal(ratio = 1) + 
    labs(x = "Sample Count", y = "Subsample Count") + 
    theme_classic(base_size = 20)

#cut <- cut(test_df$median_error, c( 0.0001, 0.001, 0.01, 0.1, 1))
ggplot(test_df, aes(x = sample_count, y = num_particles)) +
    geom_tile(aes(fill = log10(median_error_multiple)))+ 
    geom_text(aes(label = round(median_error_multiple, 3)))+ 
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                  labels = c(1, 10, 100, 1000, 10000, 100000, 1000000))+ 
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), 
                  labels = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
    scale_fill_viridis_c() + 
    coord_equal(ratio = 1) + 
    labs(x = "Sample Count", y = "Subsample Count") + 
    theme_classic()

#Plot of linear model ----
#ggplot(test_df, aes(x = Var1, y = median_error)) + geom_point(aes(color = Var2))+ scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")+ scale_fill_viridis_c()# + labs(x = "Proportion of Subsample", y = "Sample Count")
#ggplot(test_df, aes(x = sample_count, y = median_error)) #+ geom_point(aes(color = Var1))+ scale_y_log10()+ scale_x_log10()+ geom_smooth(method = "lm") + scale_fill_viridis_c()# + labs(x = "Proportion of Subsample", y = "Sample Count")
ggplot(test_df, aes(x = num_particles, y = median_error)) + 
    geom_hline(yintercept = 0.05) + 
    geom_point() + 
    scale_y_log10()+ 
    scale_x_log10() + 
    geom_smooth(method = "lm", se = F) + 
    #coord_equal() + 
    labs(x = "Number of Particles Subsampled", y = "High Absolute Error (decimal proportion)") + 
    theme_classic(base_size = 20) + 
    guides(color=guide_legend(title="Proportion Subsampled"))

ggplot(test_df) + 
    geom_hline(yintercept = 0.05) + 
    geom_point(aes(x = num_particles, y = median_error_2class), color = "blue") + 
    geom_point(aes(x = num_particles, y = median_error_4class), color = "red") + 
    geom_point(aes(x = num_particles, y = median_error_8class), color = "black") + 
    geom_point(aes(x = num_particles, y = median_error_16class), color = "green") + 
    geom_point(aes(x = num_particles, y = median_error_32class), color = "pink") + 
    scale_y_log10()+ 
    scale_x_log10() + 
    geom_smooth(aes(x = num_particles, y = median_error_2class), method = "lm", se = F, color = "blue") + 
    geom_smooth(aes(x = num_particles, y = median_error_4class), method = "lm", se = F, color = "red") + 
    geom_smooth(aes(x = num_particles, y = median_error_8class), method = "lm", se = F, color = "black") + 
    geom_smooth(aes(x = num_particles, y = median_error_16class), method = "lm", se = F, color = "green") + 
    geom_smooth(aes(x = num_particles, y = median_error_32class), method = "lm", se = F, color = "pink") + 
    # geom_smooth(method = "lm", se = F) + 
    #coord_equal() + 
    labs(x = "Number of Particles Subsampled", y = "High Relative Error (decimal proportion)") + 
    theme_classic(base_size = 20) + 
    guides(color=guide_legend(title="Proportion Subsampled"))


ggplot(test_df, aes(x = num_particles, y = median_error_multiple)) + 
    geom_hline(yintercept = 0.05) + 
    geom_point() + 
    scale_y_log10()+ 
    scale_x_log10() + 
    geom_smooth(method = "lm") + 
    labs(x = "Number of Particles Subsampled", y = "Median Uncertainty (decimal proportion)") + 
    theme_classic(base_size = 20) + 
    guides(color=guide_legend(title="Proportion Subsampled"))

model_df <- test_df %>%
    filter(num_particles != sample_count)

#Model development ----
#Wow, almost an rsq of 1, I am sure there is some basic math I have missed and that is why I didn't just do that, but this is promising. Maybe we can bring in category and improve the model some more. 
model_single <- lm(log(model_df$num_particles)~  log(model_df$median_error)) #could try ordernorm on this or something else. 
summary(model_single)             

model_multiple <- lm(log(model_df$num_particles)~  log(model_df$median_error_multiple)) #could try ordernorm on this or something else. 
summary(model_multiple)  

#Global plastic sampling. ----
#How large should a sample be? 121 particles. If so, subsample all of them. 
exp(model_single$coefficients[2] * log(0.05) + model_single$coefficients[1])
exp(model_multiple$coefficients[2] * log(0.05) + model_multiple$coefficients[1])

#How many particles do we need to sample from the whole world?
exp(model_multiple$coefficients[2] * log(0.01) + model_multiple$coefficients[1])
exp(model$coefficients[2] * log(10^-8) + model$coefficients[1])

#exact equations
exp(-2*log(0.05)-0.4)
exp(-2.1*log(0.01)-0.3)

