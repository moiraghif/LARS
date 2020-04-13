rm(list = ls())
load("simulation.RData")

library(tidyverse)


get_data <- function(ridge, lars,
              names = c(map_chr(lambdas, ~ paste0("λ=", .x)), "lars"))
    tibble(
        value = c(ridge$bias,     lars$"bias^2",
                  ridge$variance, lars$variance),
        index_name = rep(c("bias", "variance"), each = length(names)),
        model_name = factor(rep(names, times = 2), levels = names))


barplot_simulation <- function(ridge, lars)
    ggplot(get_data(ridge, lars),
           aes(x = model_name, y = value, fill = index_name)) +
        geom_bar(stat = "identity", position = "stack") +
        xlab("Model") + ylab("Value") +
        ggtitle("Reducible Error scomposition") +
        guides(fill = guide_legend(title = "")) +
        theme_minimal()


barplot_simulation(ridge_fixed, lars_fixed) %>%
    ggsave(filename = "fixed_barplot.png",
           height = 9, width = 16, dpi = 500)

barplot_simulation(ridge_random, lars_random) %>%
    ggsave(filename = "random_barplot.png",
           height = 9, width = 16, dpi = 500)


lineplot_simulation <- function(ridge, lars)
    ggplot(get_data(ridge, lars),
           aes(x = model_name, y = value,
               group = index_name, colour = index_name)) +
        geom_line() +
        geom_point()
