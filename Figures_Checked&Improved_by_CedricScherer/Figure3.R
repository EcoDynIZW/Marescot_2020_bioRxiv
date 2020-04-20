library(tidyverse)
library(rpart)
library(ggdendro)
library(ggtext)
library(showtext)

# added by SB:
library(rcartocolor)

font_add("Times New Roman", 
         regular = "C:/Windows/Fonts/times.ttf",
         bold = "C:/Windows/Fonts/timesbd.ttf",
         italic = "C:/Windows/Fonts/timesi.ttf")

results <- read.table(here::here("data", "FinalresultsFrequency.txt"), header = TRUE)
resultsDdp <- read.table(here::here("data", "FinalresultsDDP.txt"), header = TRUE)

results$transmission<-"Frequency"
resultsDdp$transmission<-"Density"

totsimu<-rbind(results, resultsDdp)

output<-as.data.frame(totsimu)
output$modularity<-0.1/output$within.group.contact

select <-
  data.frame(
    Ext = output$ext,
    `immunity.barrier` = output$immunity.barrier,
    `modulrity` = output$modularity,
    `between-group contact` = output$between.group.contact,
    `infection-length` = output$infection.length,
    virulence = output$virulence,
    `group-size` = output$group.size,
    transmission = output$transmission    
  )



gridim <- data.frame(NULL)

tabpar <-
  expand.grid(
    unique(output[, "modularity"]),
    unique(output[, "infection.length"]),
    unique(output[, "virulence"]),
    unique(output[, "group.size"]),
    unique(output[, "fertility"]),
    unique(output[,"transmission"]),
    unique(output[, "survival"])
  )
colnames(tabpar) <-
  c(
    "modularity",
    "infection.length",
    "virulence",
    "group.size",
    "fertility",
    "transmission",
    "survival"
  )

for (r in 1:nrow(tabpar)) {
  select <-
    subset(
      output,
      modularity == tabpar[r, "modularity"] &
        group.size == tabpar[r, "group.size"] &
        virulence == tabpar[r, "virulence"] &
        infection.length == tabpar[r, "infection.length"] &
        fertility == tabpar[r, "fertility"] & 
        transmission == tabpar[r, "transmission"] &
        survival == tabpar[r, "survival"]
    )
  
  extinction <-
    aggregate(select[, "ext"], list(rnot = select$rnot, age =  select[, "immunity.barrier"]), sum)
  outof <-
    aggregate(select[, "ext"], list(rnot = select$rnot, age =  select[, "immunity.barrier"]), length)
  #param<-matrix(data=rep(tabpar[r,], nrow(extinction)), byrow=FALSE, nrow = nrow(extinction), ncol= ncol(tabpar))
  gridim <-
    rbind(
      cbind(
        extinction,
        outof$x,
        rep(tabpar[r, "modularity"], nrow(extinction)),
        rep(tabpar[r, "infection.length"], nrow(extinction)),
        rep(tabpar[r, "virulence"], nrow(extinction)),
        rep(tabpar[r, "group.size"], nrow(extinction)), rep(
          tabpar[r, "fertility"], nrow(extinction)),  rep(
            tabpar[r, "transmission"], nrow(extinction)), rep(tabpar[r, "survival"], nrow(extinction))
      ), gridim)
}

colnames(gridim) = c(
  "R0",
  "immunity.barrier",
  "Extinction",
  "samplesize",
  "association",
  "infection.length",
  "virulence",
  "group.size",
  "fertility",
  "transmission",
  "survival"
)


mod <-
  glm(cbind(Extinction, samplesize) ~  immunity.barrier * (
    scale(within.group.contact)  + scale(infection.length) + scale(virulence) + scale(group.size) +  scale(survival) + scale(fertility)
  ) ^ 2, family=binomial(link='logit'),
  data = gridim)


immunity<-subset(gridim, immunity.barrier==52)
noimmunity<-subset(gridim, immunity.barrier==0)
immunity$diff<-(immunity$Extinction - noimmunity$Extinction) / 30
hist(immunity$diff)
immunity$contact.ratio<-immunity$association
immunity$trans<-as.numeric(immunity$transmission)

library(rpart)       # performing regression trees

m1 <- rpart(
  formula = diff ~ R0 + contact.ratio + infection.length + virulence + survival + fertility + group.size + transmission, data=immunity, method = "anova")


#immunity<-as.data.frame(immunity)
#tree.fit <- tree(diff ~ R0 + contact.ratio + infection.length + virulence + survival + fertility + group.size + transmission, data=immunity, method = "recursive.partition")

tree_data <- dendro_data(model = m1, type = "proportional")

## uniform labels
tree_data$labels$label <- c("7.5 < R<sub>0</sub> \u2265 7.5", 
                            "3.5 < R<sub>0</sub> \u2265 3.5", 
                            "0.0375 < fertility \u2265 0.0375", 
                            "0.0175 < fertility \u2265 0.0175", 
                            "0.025 < virulence \u2265 0.025", 
                            "15 < infection length \u2265 15", 
                            "30 < group size \u2265 30", 
                            "0.075 < virulence \u2265 0.075", 
                            "frequency = transmission = density", 
                            "0.6 < contact ratio \u2265 0.6", 
                            "0.0375 < fertility \u2265 0.0375",
                            "0.075 < virulence \u2265 0.075")

## plain-bold labels # SB changes here only
tree_data$labels$label <- c("7.5 < **R<sub>0</sub>** \u2265 7.5", # unchanged
                            "3.5 < **R<sub>0</sub>** \u2265 3.5", # unchanged
                            "0.0375 < **fertility** \u2265 0.0375", # unchanged
                            "0.0175 < **fertility** \u2265 0.0175", # unchanged
                            "0.025 > **virulence**  \u2264 0.025",   # changed, previous: "0.025 < **virulence** \u2265 0.025","
                            "15 > **infection length** \u2264 15", # changed, previous: ""15 < **infection length** \u2265 15","
                            "30 > **group size** \u2264 30", # changed, previous: "30 < **group size** \u2265 30",
                            "0.075 > **virulence** \u2264 0.075", # changed, previous: " "0.075 < **virulence** \u2265 0.075","
                            "frequency = **transmission** = density", #unchanged
                            "0.6 < **contact ratio** \u2265 0.6", #  unchanged
                            "0.0375 < **fertility** \u2265 0.0375", # unchanged
                            "0.075 > **virulence** \u2264 0.075") # changed, previous: ""0.075 < **virulence** \u2265 0.075")"

# end SB changes 

## two-colored labels
tree_data$labels$label <- c("<span style='color:#9e9e9e;'>7.5 < </span>R<sub>0</sub></b><span style='color:#9e9e9e;'> \u2265 7.5</span>", 
                            "<span style='color:#9e9e9e;'>3.5 < </span>R<sub>0</sub><span style='color:#9e9e9e;'> \u2265 3.5</span>", 
                            "<span style='color:#9e9e9e;'>0.0375 < </span>fertility<span style='color:#9e9e9e;'> \u2265 0.0375</span>", 
                            "<span style='color:#9e9e9e;'>0.0175 < </span>fertility<span style='color:#9e9e9e;'> \u2265 0.0175</span>", 
                            "<span style='color:#9e9e9e;'>0.025 < </span>virulence<span style='color:#9e9e9e;'> \u2265 0.025</span>", 
                            "<span style='color:#9e9e9e;'>15 < </span>infection length<span style='color:#9e9e9e;'> \u2265 15</span>", 
                            "<span style='color:#9e9e9e;'>30 < </span>group size<span style='color:#9e9e9e;'> \u2265 30</span>",
                            "<span style='color:#9e9e9e;'>0.075 < </span>virulence<span style='color:#9e9e9e;'> \u2265 0.075</span>", 
                            "<span style='color:#9e9e9e;'>frequency = </span>transmission<span style='color:#9e9e9e;'> = density</span>", 
                            "<span style='color:#9e9e9e;'>0.6 < </span>contact ratio<span style='color:#9e9e9e;'> \u2265 0.6</span>", 
                            "<span style='color:#9e9e9e;'>0.0375 < </span>fertility<span style='color:#9e9e9e;'> \u2265 0.0375</span>",
                            "<span style='color:#9e9e9e;'>0.075 < </span>virulence<span style='color:#9e9e9e;'> \u2265 0.075</span>")

tree_data$leaf_labels$labels <- c("<b style='font-size:10pt;'>0.0092</b><br>n=4374 (50%)",
                                  "<b style='font-size:10pt;'>0.0053</b><br>n=972 (11%)",
                                  "<b style='font-size:10pt;'>0.16</b><br>n=486 (6%)",
                                  "<b style='font-size:10pt;'>0.022</b><br>n=648 (7%)",
                                  "<b style='font-size:10pt;'>0.18</b><br>n=324 (4%)",
                                  "<b style='font-size:10pt;'>0.12</b><br>n=864 (10%)",
                                  "<b style='font-size:10pt;'>0.055</b><br>n=144 (2%)",
                                  "<b style='font-size:10pt;'>0.47</b><br>n=288 (3%)",
                                  "<b style='font-size:10pt;'>0.13</b><br>n=216 (2%)",
                                  "<b style='font-size:10pt;'>0.37</b><br>n=108 (1%)",
                                  "<b style='font-size:10pt;'>0.14</b><br>n=54 (1%)",
                                  "<b style='font-size:10pt;'>0.51</b><br>n=108 (1%)",
                                  "<b style='font-size:10pt;'>0.6</b><br>n=162 (2%)")

ggplot(segment(tree_data)) +
  geom_segment(aes(x = x, 
                   y = y, 
                   xend = xend, 
                   yend = yend), #color = n), ## for fitting with the {tree} package
               color = "grey30",
               size = 0.8) +
  geom_richtext(data = label(tree_data), 
                aes(x = x, 
                    y = y, 
                    label = label), 
                size = 3.4,
                family = "Times New Roman",
                fontface = "bold", ## comment for bold-plain labels
                hjust = .5, 
                vjust = -0.2,
                label.padding = unit(c(0.4, 0.3, 0.25, 0.3), "lines")) +
  geom_rect(data = leaf_label(tree_data), 
            aes(xmin = x - .45, xmax = x + .45, 
                ymin = y - .014, ymax = y, 
                fill = as.numeric(as.character(label)))
  ) +
  geom_richtext(data = leaf_label(tree_data), 
                aes(x = x, 
                    y = y, 
                    label = labels),
                size = 3.1,
                family = "Times New Roman",
                fill = NA, 
                label.color = NA,
                vjust = 1) +
  scale_x_continuous(expand = c(.04, .04)) +
  scale_y_continuous(expand = c(.02, .02)) +
  rcartocolor::scale_fill_carto_c(palette = "Emrld",
                                   name = "Predicted response value",
                                   limits = c(0, 1),
                                   breaks = seq(0, 1, by = .1)) +
  guides(fill = guide_colorbar(direction = "horizontal",
                                title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom",
                                label.hjust = .5,
                                barwidth = unit(15, "lines"),
                                barheight = unit(.3, "lines"))) +
  theme_dendro() +
  theme(legend.position = c(.78, .8),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 11),
        legend.text = element_text(family = "Times New Roman", size = 10),
        plot.margin = margin(0, 0, -10, 0))

ggsave(here::here("figures", "fig_3_face_SB.pdf"), width = 12, height = 8.5, device = cairo_pdf)

