finite_corrected_eq <- function(pop_size, uncorrected){
  round(pop_size * uncorrected / (uncorrected + pop_size - 1), 0)
} 

uncorrected_eq <- function(confidence = 0.95, groups = 1, expected_p = 0.5, max_error = 0.05){
  critical_val = qnorm(p=(1-confidence^(1/groups))/2, lower.tail=FALSE) # The group root is a new thing here. 
  round(critical_val^2*expected_p*(1-expected_p)/max_error^2, 0)
}

round(critical_val^2*expected_p*(1-expected_p)/max_error^2, 0)
1.96^2 * 0.5 * (0.5)
uncorrected_eq()
1/0.05^2
0.05^-2
100^(-1/2)

library(dplyr)
library(ggplot2)

options(scipen = 99)
population_sizes <- 10^(1:10)
groups <- 1:5
max_error <- c(10^(c(-5:-1)), 0.05)

test_ungrouped <- expand.grid(population_sizes = population_sizes, max_error = max_error) %>%
  mutate(sample_size = finite_corrected_eq(pop_size = population_sizes, uncorrected = uncorrected_eq(max_error = max_error)))

test_grouped <- expand.grid(population_sizes = population_sizes, max_error = max_error) %>%
  mutate(sample_size = finite_corrected_eq(pop_size = population_sizes, uncorrected = uncorrected_eq(max_error = max_error, groups = 4)))

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
  geom_abline(slope=1, intercept= 0) +
  theme_dark(base_size = 20) +
  coord_equal() +
  labs(x = "Simulated Sample Counts", y = "Calculated Sample Counts")



#final recommendations for sample size and subsample size.
uncorrected_eq()
uncorrected_eq(max_error = 0.1)

#final recommendation for all group analysis at once.
uncorrected_eq(groups = 4)



finite_corrected_eq(pop_size = 100, uncorrected = )

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
