library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv)


#give classes
#randomly select some Pecent
#compare that subsample to the full sample.
#sum up all differences in absolute value space. 
data <- tibble(
    class_num = numeric(),
    subsample_size = numeric(),
    sample_count = numeric(),
    error = numeric()
)

boot_mean_error <- function(sample_subsample_size, sample_count){
    #Disallow sample count to if(sample_count < sample_classes)
    #Disallow sample round(sample_subsample_size*sample_count, 0) < 1)
    #, prob = c(0.5, 0.2, rep(0.3/8, 8)) if I add this to the first sample argument, we get a less conservative estimate, the errors go down. I think it is because you are less likely to mess up a skewed distribution than a uniform one. 
    
    b = 100 #Changing b doesn't do much, probably more accurate with larger number but this will save time while in dev mode.
    error <- numeric(length = b)
    
    set.seed(37)
    
    for(n in 1:b){
        sample_classes <- sample(2:10, 1)    #Sample classes must be greater than 1 but less than 20
        particle_categories <- 1:sample_classes
        
        #Simulation
        particles <- sample(particle_categories, size = sample_count, replace = T) #could add weights to this so that 1 or two of them always have a big sway. But we dont know that for sure yet.
        
        subsetparticles <- sample(particles, round(sample_subsample_size*length(particles), 0))
        
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
    
    error
}

errors <- boot_mean_error(sample_count = 10, sample_subsample_size = 0.1) #Just say we will set this to 20, that is really what can be reliably characterized, can do a quick lit search to figure out how people are reporting it.    

hist(errors)
quantile(errors, c(0.025, 0.5, 0.975))
mean(errors)

#table(sample(particle_categories, size = sample_count, replace = T))

subsample_size <- c(seq(0.001, 0.01, by = 0.001), seq(0.01, 0.1, by = 0.01), seq(0.1, 0.9, by = 0.1))
count <- c(seq(10, 100, by = 10), seq(100, 1000, by = 100), seq(1000, 10000, by = 1000), seq(10000, 100000, by = 10000))

test_df <- expand.grid(subsample_size, count) %>%
    filter(Var1 * Var2 >= 1)

for(n in 1:nrow(test_df)){
    test_df[n, "median_error"] <- quantile(
        boot_mean_error(sample_count = unlist(test_df[n,"Var2"]), 
                        sample_subsample_size = unlist(test_df[n,"Var1"])),
        c(0.95)
    )
}

test_df <- test_df %>%
    #mutate(highly_accurate = median_error < 0.05, cut = cut(median_error, breaks = c(1, 0.1, 0.01, 0.001, 0.0001))) %>%
    mutate(num_particles = Var1*Var2)

ggplot(test_df, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = cut))+ scale_y_log10()+ scale_x_log10() + scale_fill_viridis_d() + labs(x = "Proportion of Subsample", y = "Sample Count") + theme_classic()
ggplot(test_df, aes(x = Var1, y = median_error)) + geom_point(aes(color = Var2))+ scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")+ scale_fill_viridis_c()# + labs(x = "Proportion of Subsample", y = "Sample Count")
ggplot(test_df, aes(x = Var2, y = median_error)) + geom_point(aes(color = Var1))+ scale_y_log10()+ scale_x_log10()+ geom_smooth(method = "lm") + scale_fill_viridis_c()# + labs(x = "Proportion of Subsample", y = "Sample Count")
ggplot(test_df, aes(x = num_particles, y = median_error)) + 
    geom_hline(yintercept = 0.05) + 
    geom_point(aes(color = Var1)) + 
    scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1), labels = c(0.0001, 0.001, 0.01, 0.1, 1), limits = c(0.0001, 1))+ 
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = c(1, 10, 100, 1000, 10000, 100000)) + 
    labs(x = "Number of Particles Subsampled", y = "Median Uncertainty (decimal proportion)") + 
    theme_classic(base_size = 20) + 
    guides(color=guide_legend(title="Proportion Subsampled"))

#Wow, almost an rsq of 1, I am sure there is some basic math I have missed and that is why I didn't just do that, but this is promising. Maybe we can bring in category and improve the model some more. 
model <- lm(log(test_df$median_error) ~ log(test_df$Var2) + log(test_df$num_particles)) #could try ordernorm on this or something else. 
summary(model)             

#Parameters
count_mod = 10000
#error_mod = 0.05
subsample_mod = 0.01

#Var1 is subsample_size, Var2 is count, This seems to be corresponding well to the model
10^(model$coefficients[2] * log10(subsample_mod) + model$coefficients[3] * log10(count_mod) + model$coefficients[1])

subsample_size_mod <- seq(0.001, 0.999, by = 0.001)
count_mod <- 10:100000
figure_mod_df <- expand.grid(subsample_size_mod, count_mod) %>%
    mutate(error = 10^(model$coefficients[2] * log10(Var1) + model$coefficients[3] * log10(Var2) + model$coefficients[1]))
ggplot(figure_mod_df, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = error))+ scale_y_log10() + scale_fill_viridis_c() + labs(x = "Proportion of Subsample", y = "Sample Count")

#e^((log(Uncertainty) - 0.100639 * log(proportion sampled) + 0.967713) /  -0.609322 )
exp((log(0.05) - 0.10 * log(10^20) + 0.97) /  -0.61)

#Flip this around to predict the count

#I don't think I calculated this right. The result doesn't correspond to my expectations. 
10^((log10(error_mod)-model$coefficients[1]) / (model$coefficients[3] * log10(count_mod) * model$coefficients[2]))*count_mod

test_df$logerror <- log10(test_df$median_error)
modelgam <- gam(logerror ~ s(log10(Var1)) + s(log10(Var2)), data = test_df)
