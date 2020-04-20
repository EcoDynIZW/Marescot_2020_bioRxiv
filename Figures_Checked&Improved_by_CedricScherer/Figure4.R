library(tidyverse)
library(ggtext)
library(patchwork)
library(showtext)

font_add("Times New Roman", 
         regular = "C:/Windows/Fonts/times.ttf",
         bold = "C:/Windows/Fonts/timesbd.ttf",
         italic = "C:/Windows/Fonts/timesi.ttf")


results <- read.table(here::here("data", "ImmunityV2.txt"), header = TRUE) %>% 
  mutate(type = "A", step = row_number())
resultsNO <- read.table(here::here("data", "NoImmunityV2.txt"), header = TRUE) %>% 
  mutate(type = "B", step = row_number())


results_all <- rbind(results, resultsNO) %>% 
  dplyr::select(step, type, abundance, susjuv, susad, infecjuv, infecad, recovjuv, recovad) %>% 
  pivot_longer(
    cols = susjuv:recovad,
    names_to = "group",
    values_to = "inds"
  ) %>% 
  mutate(SIR = case_when(
    str_detect(group, "sus") ~ "S",
    str_detect(group, "infec") ~ "I",
    str_detect(group, "recov") ~ "R"
  ))


####################################################################################################################
## PLOT ############################################################################################################
####################################################################################################################

theme_lucile <-
  theme_classic(base_size = 18, base_family = "Times New Roman") +
  theme(axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_markdown(hjust = .4, margin = margin(r = 15)),
        axis.line = element_line(colour = 'black', size = .8),
        axis.ticks = element_line(size = .8), 
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(hjust = .5, face = "bold"),
        legend.position = "none")


## Susceptible
lab <-
  results_all %>% 
  filter(SIR == "S", type == "A", step == 365) %>% 
  mutate(
    step = if_else(group == "susjuv", step + 0, step + 10),
    inds = if_else(group == "susjuv", inds + 20, inds - 10),
    label = if_else(group == "susjuv", "Juveniles", "Adults"),
    label_x = step + 50,
    label_y = if_else(group == "susjuv", inds + 170, inds - 170)
  )

plot_a_proto <- 
  results_all %>% 
  filter(SIR == "S", type == "A") %>% 
  ggplot(aes(step, inds)) +
    geom_line(
      aes(
        color = group, 
        linetype = group
      ),
      size = 1.3
    ) +
    geom_segment(
      data = lab,
      aes(
        x = step,
        xend = label_x + .1,
        y = inds,
        yend = label_y,
        color = group
      ),
      size = .8
    ) +
    scale_x_continuous(
      expand = c(.015, .015),
      breaks = seq(0, 520, by = 52)
    ) +
    scale_y_continuous(
      expand = c(.03, .03),
      limits = c(0, 850), 
      breaks = seq(0, 800, by = 200)
    ) +
    scale_color_manual(values = c(colorspace::darken("#0090f7", .4), "#0090f7")) +
    scale_fill_manual(values = c(colorspace::darken("#0090f7", .4), "#0090f7")) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    theme_lucile +
    labs(
      x = "Time (weeks)", 
      y = "&#35; of susceptibles&nbsp;&nbsp;*S*", 
      title = "ORC",
      tag = expression(paste("(", italic("a"), ")"))
    )

plot_a <- 
  plot_a_proto +
  geom_label(
    data = lab,
    aes(
      label_x, label_y, 
      label = label,
      fill = group
    ),
    family = "Times New Roman",
    fontface = "bold",
    color = "white",
    size = 4.5,
    hjust = 0,
    label.r = unit(.6, "lines"),
    label.padding = unit(9, "pt"),
    label.size = 0
  )


plot_a2 <- 
  plot_a_proto +
  geom_label(
    data = lab,
    aes(
      label_x, label_y, 
      label = label,
      color = group
    ),
    family = "Times New Roman",
    fontface = "bold",
    fill = "white",
    size = 4.5,
    hjust = 0,
    label.r = unit(.4, "lines"),
    label.padding = unit(9, "pt"),
    label.size = 1
  ) +
  geom_label(
    data = lab,
    aes(
      label_x, label_y, 
      label = label
    ),
    family = "Times New Roman",
    fontface = "bold",
    color = "black",
    fill = "transparent",
    size = 4.5,
    hjust = 0,
    label.r = unit(.4, "lines"),
    label.padding = unit(9, "pt"),
    label.size = 0
  ) 

plot_b <-
  results_all %>% 
  filter(SIR == "S", type == "B") %>% 
  ggplot(
    aes(
      step, inds, 
      color = group, 
      linetype = group
      )
    ) +
    geom_line(size = 1.3) +
    scale_x_continuous(
      expand = c(.015, .015),
      breaks = seq(0, 520, by = 52)
    ) +
    scale_y_continuous(
      expand = c(.03, .03),
      limits = c(0, 850), 
      breaks = seq(0, 800, by = 200)
    ) +
    scale_color_manual(values = c(colorspace::darken("#0090f7", .4), "#0090f7")) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    theme_lucile +
    labs(
      x = "Time (weeks)", 
      y = "", 
      title = "Baseline scenario",
      tag = expression(paste("(", italic("b"), ")"))
    )


## Infected
plot_c <- 
  results_all %>% 
  filter(SIR == "I", type == "A") %>% 
  ggplot(
    aes(
      step, inds, 
      color = group, 
      linetype = group
    )
  ) +
  geom_line(size = 1.3) +
  scale_x_continuous(
    expand = c(.015, .015),
    breaks = seq(0, 520, by = 52)
  ) +
  scale_y_continuous(
    expand = c(.03, .03),
    limits = c(0, 420), 
    breaks = seq(0, 400, by = 100)
  ) +
  scale_color_manual(values = c(colorspace::darken("#f20000", .4), "#f20000")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  theme_lucile +
  labs(
    x = "Time (weeks)", 
    y = "&#35; of infected&nbsp;&nbsp;*I*", 
    title = NULL,
    tag = expression(paste("(", italic("c"), ")"))
  )


plot_d <-
  results_all %>% 
  filter(SIR == "I", type == "B") %>% 
  ggplot(
    aes(
      step, inds, 
      color = group, 
      linetype = group
    )
  ) +
  geom_line(size = 1.3) +
  scale_x_continuous(
    expand = c(.015, .015),
    breaks = seq(0, 520, by = 52)
  ) +
  scale_y_continuous(
    expand = c(.03, .03),
    limits = c(0, 420), 
    breaks = seq(0, 400, by = 100)
  ) +
  scale_color_manual(values = c(colorspace::darken("#f20000", .4), "#f20000")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  theme_lucile +
  labs(
    x = "Time (weeks)", 
    y = "", 
    title = NULL,
    tag = expression(paste("(", italic("d"), ")"))
  )


## Recovered
plot_e <- 
  results_all %>% 
  filter(SIR == "R", type == "A") %>% 
  ggplot(
    aes(
      step, inds, 
      color = group, 
      linetype = group
    )
  ) +
  geom_line(size = 1.3) +
  scale_x_continuous(
    expand = c(.015, .015),
    breaks = seq(0, 520, by = 52)
  ) +
  scale_y_continuous(
    expand = c(.03, .03),
    limits = c(0, 850), 
    breaks = seq(0, 800, by = 200)
  ) +
  scale_color_manual(values = c(colorspace::darken("#00de6f", .4), "#00de6f")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  theme_lucile +
  labs(
    x = "Time (weeks)", 
    y = "&#35; of recovered&nbsp;&nbsp;*R*", 
    title = NULL,
    tag = expression(paste("(", italic("e"), ")"))
  )


plot_f <-
  results_all %>% 
  filter(SIR == "R", type == "B") %>% 
  ggplot(
    aes(
      step, inds, 
      color = group, 
      linetype = group
    )
  ) +
  geom_line(size = 1.3) +
  scale_x_continuous(
    expand = c(.015, .015),
    breaks = seq(0, 520, by = 52)
  ) +
  scale_y_continuous(
    expand = c(.03, .03),
    limits = c(0, 850), 
    breaks = seq(0, 800, by = 200)
  ) +
  scale_color_manual(values = c(colorspace::darken("#00de6f", .4), "#00de6f")) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  theme_lucile +
  labs(
    x = "Time (weeks)", 
    y = "", 
    title = NULL,
    tag = expression(paste("(", italic("f"), ")"))
  )


## FULL PANEL
(plot_a + plot_b) / (plot_c + plot_d) / (plot_e + plot_f) + plot_layout(heights = c(1, .97, .97))
ggsave(here::here("figures", "fig_4.pdf"), width = 16, height = 11, device = cairo_pdf)

(plot_a2 + plot_b) / (plot_c + plot_d) / (plot_e + plot_f) + plot_layout(heights = c(1, .97, .97))
ggsave(here::here("figures", "fig_4_alt.pdf"), width = 16, height = 11, device = cairo_pdf)
