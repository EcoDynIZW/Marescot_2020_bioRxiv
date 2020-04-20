library(ggplot2)
library(showtext)

font_add("Times New Roman", 
         regular = "C:/Windows/Fonts/times.ttf",
         bold = "C:/Windows/Fonts/timesbd.ttf",
         italic = "C:/Windows/Fonts/timesi.ttf")

results <- read.table(here::here("data", "FinalresultsFrequency.txt"), header = TRUE)
resultsDdp <- read.table(here::here("data", "FinalresultsDDP.txt"), header = TRUE)

results$transmission<-"Frequency"
resultsDdp$transmission<-"Density"

totfreq<-subset(results, rnot==20)
totDdp<-subset(resultsDdp, rnot==20)
totsimu<-rbind(results, resultsDdp)

output<-as.data.frame(totsimu)

select <-
  data.frame(
    Ext = output$ext,
    `immunity.barrier` = output$immunity.barrier,
    `within-group contact` = output$within.group.contact,
    `between-group contact` = output$between.group.contact,
    `infection-length` = output$infection.length,
    virulence = output$virulence,
    `group-size` = output$group.size
  )

gridim <- data.frame(NULL)

tabpar <-
  expand.grid(
    unique(output[, "within.group.contact"]),
    unique(output[, "infection.length"]),
    unique(output[, "virulence"]),
    unique(output[, "group.size"]),
    unique(output[, "fertility"]),
    unique(output[, "survival"])
  )

colnames(tabpar) <-
  c(
    "within.group.contact",
    "infection.length",
    "virulence",
    "group.size",
    "fertility",
    "survival"
  )

for (r in 1:nrow(tabpar))
{
  select <-
    subset(
      output,
      within.group.contact == tabpar[r, "within.group.contact"] &
        group.size == tabpar[r, "group.size"] &
        virulence == tabpar[r, "virulence"] &
        infection.length == tabpar[r, "infection.length"] &
        fertility == tabpar[r, "fertility"] &
        survival == tabpar[r, "survival"]
    )
  
  extinction <-
    aggregate(select[, "ext"], list(rnot = select$rnot, age =  select[, "immunity.barrier"]), sum)
  outof <-
    aggregate(select[, "ext"], list(rnot = select$rnot, age =  select[, "immunity.barrier"]), length)
  gridim <-
    rbind(
      cbind(
        extinction,
        outof$x,
        rep(tabpar[r, "within.group.contact"], nrow(extinction)),
        rep(tabpar[r, "infection.length"], nrow(extinction)),
        rep(tabpar[r, "virulence"], nrow(extinction)),
        rep(tabpar[r, "group.size"], nrow(extinction)), rep(
          tabpar[r, "fertility"], nrow(extinction)), rep(tabpar[r, "survival"], nrow(extinction))
      ), gridim)
}

colnames(gridim) = c(
  "R0",
  "immunity.barrier",
  "Extinction",
  "samplesize",
  "within.group.contact",
  "infection.length",
  "virulence",
  "group.size",
  "fertility",
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

totsimu<-rbind(results, resultsDdp)

output<-as.data.frame(totsimu)
output[which(output$ext == "persistence" ),"ext"]<- 0
output[which(output$ext == "extinct"),"ext"]<- 1

select<-output
gridim<-aggregate(select[,"ext"], list(rnot=select$rnot, age=  select[,"immunity.barrier"]),mean)
colnames(gridim) = c("Basic reproduction number R0", "Between-group contact", "Epidemic fade-out probability")

gridim[which(gridim$`Between-group contact`==52), "Between-group contact"]<-"ORC"
gridim[which(gridim$`Between-group contact`==0), "Between-group contact"]<-"Baseline scenario" 


####################################################################################################################
## PLOT ############################################################################################################
####################################################################################################################

theme_lucile <-
  theme_classic(base_size = 16, base_family = "Times New Roman") +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        axis.line = element_line(colour = 'black', size = .8),
        axis.ticks = element_line(size = .8), 
        axis.ticks.length = unit(.3, "cm"),
        legend.position = c(.7, .85),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.key.width = unit(3, "lines"),
        legend.key.height = unit(1.3, "lines"))

ggplot(gridim, 
  aes(`Basic reproduction number R0`, 
      `Epidemic fade-out probability`, 
       color = `Between-group contact`,
       shape = `Between-group contact`)) + 
  geom_line(size = 1.5, alpha = .7) + 
  geom_point(size = 2.8, stroke = 1.7, fill = "white") +
  geom_point(size = 2.8, stroke = 1.7) +
  scale_color_manual(values = c("#7d408a", "#eb6f6c"))+
  scale_shape_manual(values = c(16, 21))+
  labs(x = expression("Basic reproduction number R"[0]),
       y = "Epidemic fade-out probability") +
  theme_lucile

ggsave(here::here("figures", "fig_1.pdf"), width = 8, height = 5, device = cairo_pdf)
