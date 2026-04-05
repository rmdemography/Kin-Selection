#------------------------------------------------------------------------------#
# Paper:  Kin selection and population dynamics
# Title:  Visualisation
# Data:   HMD & HFD
# Author: Rahul Mondal
# Date:   04/04/2026
#------------------------------------------------------------------------------#
pacman::p_load(tidyvers, ggplot2, RColorBrewer)
load("sens_lro.RData")
# Stochastic Growth Rate

# LRO
sens_mean$beta <- factor(
  sens_mean$beta,
  levels = c("gm", "os", "age_m"),
  labels = c("beta[gm]", "beta[os]", "beta[ma]")
)
display.brewer.all(colorblindFriendly = TRUE)
col1 <- c(
  brewer.pal(9, "YlGn")[9],
  brewer.pal(9, "Blues")[9],
  brewer.pal(9, "PuRd")[9]
)
col1 <- c(
  brewer.pal(9, "YlOrBr")[6],
  brewer.pal(9, "OrRd")[8],
  brewer.pal(9, "Reds")[9]
)
ggplot(sens_mean, aes(x = country, y = sensitivity, fill = country, color = country)) +
  geom_bar(stat = "identity", linewidth = 0.7, width = 0.6) +
  facet_wrap(~ beta, labeller = label_parsed) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray40") +
  scale_x_discrete(labels = c("FRATNP" = "France", "JPN" = "Japan", "SWE" = "Sweden")) +
  scale_fill_manual(
    values = c("FRATNP" = alpha(col1[1], 0.6), "JPN" = alpha(col1[2], 0.6), "SWE" = alpha(col1[3], 0.6)),
    labels = c("FRATNP" = "France", "JPN" = "Japan", "SWE" = "Sweden")
  ) +
  scale_color_manual(
    values = c("FRATNP" = col1[1], "JPN" = col1[2], "SWE" = col1[3]),
    labels = c("FRATNP" = "France", "JPN" = "Japan", "SWE" = "Sweden")
  ) +
  labs(
    x     = NULL,
    fill  = "Country",
    color = "Country",
    y     = expression("Sensitivity" ~ (d*tilde(rho)[1] / d*beta[k]))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    panel.border = element_rect(fill = NA, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(1,1,1,1,"cm"),
    plot.background = element_blank(),
    legend.background = element_blank(),
    panel.background = element_blank()
  )
ggsave("fig2.png", width = 20, height = 14, units = "cm", dpi = 300)
