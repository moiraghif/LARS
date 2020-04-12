library(tidyverse)



names <- c(map_chr(lambdas, ~ paste0("λ=", .x)), "lars")


bias_variance <- tibble(
   value = c(ridge_simulations.bias,     lars_simulation$"bias^2",
             ridge_simulations.variance, lars_simulation$variance),
   index_name = rep(c("bias", "variance"), each = length(names)),
   model_name = factor(rep(names, times = 2), levels = names))


ggplot(bias_variance, aes(x = model_name, fill = index_name)) +
    geom_bar(aes(y = value), stat = "identity", position = "stack") +
    xlab("") + ylab("") + ggtitle("Reducible Error scomposition") +
    guides(fill = guide_legend(title = "")) +
    scale_y_sqrt()
