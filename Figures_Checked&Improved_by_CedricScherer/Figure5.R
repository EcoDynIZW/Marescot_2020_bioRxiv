library(tidyverse)
library(ggtext)
library(patchwork)
library(showtext)

font_add("Times New Roman", 
         regular = "C:/Windows/Fonts/times.ttf",
         bold = "C:/Windows/Fonts/timesbd.ttf",
         italic = "C:/Windows/Fonts/timesi.ttf")


results <- read.table(here::here("data", "ImmunityV2.txt"), header = TRUE)
resultsNO <- read.table(here::here("data", "NoImmunityV2.txt"), header = TRUE)


theme_lucile <-
  theme_classic(base_size = 18, base_family = "Times New Roman") +
  theme(axis.title.x = element_text(size = 19, margin = margin(t = 15)),
        axis.title.y = element_text(size = 19, margin = margin(r = 8)),
        axis.text = element_text(size = 16),
        axis.ticks = element_line(size = .8), 
        axis.ticks.length = unit(.2, "cm"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = .8),
        plot.title = element_text(face = "bold", size = 24, margin = margin(0, 0, 10, 0), hjust = .5),
        plot.tag = element_text(face = "italic", size = 24),
        legend.position = "none")


results$rest <- results$sealed - results$protected - results$locked
results$open <- 20 - results$sealed
vec<-matrix(0, nrow=(length(results$rest)*4), ncol=2)

r <- 1
for(i in 1:length(results$rest)) {
  vec[r:(r+3),1] <- rep(i,4)
  vec[r:(r+3),2] <- t(results[i,c("locked", "protected", "rest", "open")])
  r <- r + 4
}  


vecbis_a <- NULL
for(i in 1:200){
  vec_a <- rep(1, results[i,"sealed"])
  vec_a2 <- rep(2, results[i,"protected"])
  vec_a3 <- rep(3, (results[i,"open"]+ results[i,"rest"]))
  
  
  bound <- c(vec_a, vec_a2,vec_a3)
  vecbis_a <- rbind(vecbis_a, cbind(rep(i,length(bound)), bound))
}  


df_vecbis_a <- as_tibble(vecbis_a) %>% 
  mutate(bound = factor(bound, levels = c("3", "2", "1")))



resultsNO$rest <- resultsNO$sealed - resultsNO$protected - resultsNO$locked
resultsNO$open <- 20 - resultsNO$sealed
vec<-matrix(0, nrow=(length(resultsNO$rest)*4), ncol=2)

r <- 1
for(i in 1:length(resultsNO$rest)){
  vec[r:(r+3),1] <- rep(i,4)
  vec[r:(r+3),2] <- t(resultsNO[i,c("locked", "protected", "rest", "open")])
  r <- r + 4
}  

vecbis_b<-NULL
for(i in 1:200){
  vec_b <- rep(1, resultsNO[i,"sealed"])
  vec_b2 <- rep(2, resultsNO[i,"protected"])
  vec_b3 <- rep(3, (resultsNO[i,"open"]+ resultsNO[i,"rest"]))

  bound <- c(vec_b, vec_b2,vec_b3)
  vecbis_b <- rbind(vecbis_b, cbind(rep(i,length(bound)), bound))
}  


df_vecbis_b <- as_tibble(vecbis_b) %>% 
  mutate(bound = factor(bound, levels = c("3", "2", "1")))


a <- 
  ggplot(tibble(x = 1:10, y = 1:10),
         aes(x, y)) +
    labs(tag = expression(paste("(", italic("a"), ")"))) +
    theme_void() +
    theme(plot.tag = element_text(family = "Times New Roman",
                                  face = "italic", size = 24),)
    

labs_b <-
  tibble(
    x = c(53, 53, 125),
    y = c(.23, .7, .5),
    label = c("**sealed groups**<br>(recovered adults and<br>infected juveniles)",
              "**open groups**<br>(susceptible and/or<br>infected adults)",
              "**sealed groups**<br>(recovered adults and<br>susceptible juveniles)")
  )

b <- 
  ggplot(df_vecbis_a, aes(V1, ..count..)) +
    geom_density(aes(fill = bound),
                 position = "fill", color = NA) +
    geom_richtext(data = labs_b, 
                  aes(x, y, label = label),
                  family = "Times New Roman",
                  size = 4.5,
                  fill = NA, 
                  label.color = NA) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(0, 200)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1)) +
    scale_fill_manual(values = c("#83c2ca", "#ebce77", "#b96777")) +
    #scale_fill_manual(values = c("grey85", "grey65", "grey45")) +
    labs(x = NULL, 
         y = "Proportion of groups", 
         title = "ORC", 
         tag = expression(paste("(", italic("b"), ")"))) +
    theme_lucile

labs_c <-
  tibble(
    x = c(47, 47, 125),
    y = c(.22, .62, .77),
    label = c("**open groups**<br>(recovered adults and<br>infected juveniles)",
              "**open groups**<br>(susceptible and/or<br>infected adults)",
              "**open groups**<br>(recovered adults and<br>susceptible juveniles)")
  )

c <- 
  ggplot(df_vecbis_b, aes(V1, ..count.., fill = bound)) +
    geom_density(position = "fill", color = NA) +
    geom_richtext(data = labs_c, 
                  aes(x, y, label = label),
                  family = "Times New Roman",
                  size = 4.5,
                  fill = NA, 
                  label.color = NA) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(0, 200)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1)) +
    scale_fill_manual(values = c("#83c2ca", "#ebce77", "#b96777")) +
    #scale_fill_manual(values = c("grey85", "grey65", "grey45")) +
    labs(x = "Weeks", 
         y = "Proportion of groups", 
         title = "Baseline scenario", 
         tag = expression(paste("(", italic("c"), ")"))) +
    theme_lucile


a / b / c

ggsave(here::here("figures", "fig_5.pdf"), width = 10, height = 13.5, device = cairo_pdf)  
