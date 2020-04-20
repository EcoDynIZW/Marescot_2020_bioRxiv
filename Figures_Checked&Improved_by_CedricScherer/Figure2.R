library(tidyverse)
library(patchwork)
library(showtext)

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
output$modularity<-0.10/output$within.group.contact


##################mean difference in epidemic fadeout between age-dependent and age-independent networks

path <- here::here("data", "fig2_immunity.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)
  
  tabpar <-
    expand.grid(
      unique(output[, "modularity"]),
      unique(output[, "infection.length"]),
      unique(output[, "virulence"]),
      unique(output[, "group.size"]),
      unique(output[, "fertility"]),
      unique(output[, "survival"]),
      unique(output[, "rnot"]),
      unique(output[, "transmission"])
    )
  colnames(tabpar) <-
    c(
      "modularity",
      "infection.length",
      "virulence",
      "group.size",
      "fertility",
      "survival",
      "rnot",
      "transmission" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
        modularity == tabpar[r, "modularity"] &
          group.size == tabpar[r, "group.size"] &
          virulence == tabpar[r, "virulence"] &
          infection.length == tabpar[r, "infection.length"] &
          fertility == tabpar[r, "fertility"] &
          survival == tabpar[r, "survival"] &
          rnot == tabpar[r, "rnot"] &
          transmission == tabpar[r, "transmission"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(immunity.barrier =  select[, "immunity.barrier"]), sum)
    outof <-
      aggregate(select[, "ext"], list(immunity.barrier =  select[, "immunity.barrier"]), length)
    persistence<-outof - extinction
    b<-matrix(c(extinction$x, persistence$x), nrow = 2, dimnames = list(immunity.barrier = extinction$immunity.barrier, Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
      gridim <-
        rbind(
          cbind(
            extinction,
            rep(Fisher$p.value, nrow(extinction)),
            rep(tabpar[r, "modularity"], nrow(extinction)),
            rep(tabpar[r, "infection.length"], nrow(extinction)),
            rep(tabpar[r, "virulence"], nrow(extinction)),
            rep(tabpar[r, "group.size"], nrow(extinction)),  rep(
              tabpar[r, "fertility"], nrow(extinction)), rep(tabpar[r, "survival"], nrow(extinction)), rep(tabpar[r, "rnot"], nrow(extinction)), rep(tabpar[r, "transmission"], nrow(extinction))), gridim)
    }
    
  }
  colnames(gridim) = c(
    "immunity",
    "Extinction",
    "significant",
    "modularity",
    "infection.length",
    "virulence",
    "group.size",
    "fertility",
    "survival",
    "rnot",
    "transmission"
  )
  
  immunity<-subset(gridim, immunity==52)
  noimmunity<-subset(gridim, immunity==0)
  immunity$difference<-(immunity$Extinction - noimmunity$Extinction) / 30
  
  saveRDS(immunity, path)
  
}else{
  immunity <- readRDS(path)
}


################## mean difference in epidemic fadeout between modular and none modular networks####################

path <- here::here("data", "fig2_modularity.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)
  
  tabpar <-
    expand.grid(
      unique(output[, "immunity.barrier"]),
      unique(output[, "infection.length"]),
      unique(output[, "virulence"]),
      unique(output[, "group.size"]),
      unique(output[, "fertility"]),
      unique(output[, "survival"]),
      unique(output[, "rnot"]),
      unique(output[, "transmission"])
    )
  colnames(tabpar) <-
    c(
      "immunity.barrier",
      "infection.length",
      "virulence",
      "group.size",
      "fertility",
      "survival",
      "rnot",
      "transmission" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
        immunity.barrier == tabpar[r, "immunity.barrier"] &
          group.size == tabpar[r, "group.size"] &
          virulence == tabpar[r, "virulence"] &
          infection.length == tabpar[r, "infection.length"] &
          fertility == tabpar[r, "fertility"] &
          survival == tabpar[r, "survival"] &
          rnot == tabpar[r, "rnot"] &
          transmission == tabpar[r, "transmission"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(modularity =  select[, "modularity"]), sum)
    outof <-
      aggregate(select[, "ext"], list(modularity =  select[, "modularity"]), length)
    persistence<-outof - extinction
    b<-matrix(c(extinction$x, persistence$x), nrow = 3, dimnames = list(Modularity = c(10, 50, 100), Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
    gridim <-
      rbind(
        cbind(
          extinction,
          rep(Fisher$p.value, nrow(extinction)),
          rep(tabpar[r, "immunity.barrier"], nrow(extinction)),
          rep(tabpar[r, "infection.length"], nrow(extinction)),
          rep(tabpar[r, "virulence"], nrow(extinction)),
          rep(tabpar[r, "group.size"], nrow(extinction)), rep(
            tabpar[r, "fertility"], nrow(extinction)), rep(tabpar[r, "survival"], nrow(extinction)), rep(tabpar[r, "rnot"], nrow(extinction)), rep(tabpar[r, "transmission"], nrow(extinction))), gridim)
    }
    
  }
  colnames(gridim) = c(
    "modularity",
    "Extinction",
    "significant",
    "immunity.barrier",
    "infection.length",
    "virulence",
    "group.size",
    "fertility",
    "survival",
    "rnot",
    "transmission"
  )
  
  modularity<-subset(gridim, modularity==0.1)
  nomodularity<-subset(gridim, modularity==1)
  modularity$diff<-(modularity$Extinction - nomodularity$Extinction) / 30

  saveRDS(modularity, path)
  
}else{
  modularity <- readRDS(path)
}


######### mean difference in epidemic fadeout between networks with high and low survival ############

path <- here::here("data", "fig2_survival.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)
  
  tabpar <-
    expand.grid(
      unique(output[, "immunity.barrier"]),
      unique(output[, "infection.length"]),
      unique(output[, "virulence"]),
      unique(output[, "group.size"]),
      unique(output[, "fertility"]),
      unique(output[, "modularity"]),
      unique(output[, "rnot"]),
      unique(output[, "transmission"])
    )
  colnames(tabpar) <-
    c(
      "immunity.barrier",
      "infection.length",
      "virulence",
      "group.size",
      "fertility",
      "modularity",
      "rnot",
      "transmission" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
        immunity.barrier == tabpar[r, "immunity.barrier"] &
          group.size == tabpar[r, "group.size"] &
          virulence == tabpar[r, "virulence"] &
          infection.length == tabpar[r, "infection.length"] &
          fertility == tabpar[r, "fertility"] &
          modularity == tabpar[r, "modularity"] &
          rnot == tabpar[r, "rnot"] &
          transmission == tabpar[r, "transmission"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(survial =  select[, "survival"]), sum)
    outof <-
      aggregate(select[, "ext"], list(survival =  select[, "survival"]), length)
    persistence<-outof - extinction
    b<-matrix(c(extinction$x, persistence$x), nrow = 3, dimnames = list(survival = extinction$survival, Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
      gridim <-
        rbind(
          cbind(
            extinction,
            rep(Fisher$p.value, nrow(extinction)),
            rep(tabpar[r, "immunity.barrier"], nrow(extinction)),
            rep(tabpar[r, "infection.length"], nrow(extinction)),
            rep(tabpar[r, "virulence"], nrow(extinction)),
            rep(tabpar[r, "group.size"], nrow(extinction)),  rep(
              tabpar[r, "fertility"], nrow(extinction)), rep(tabpar[r, "modularity"], nrow(extinction)), rep(tabpar[r, "rnot"], nrow(extinction)), rep(tabpar[r, "transmission"], nrow(extinction))), gridim)
    }
    
  }
  colnames(gridim) = c(
    "survival",
    "Extinction",
    "significant",
    "immunity.barrier",
    "infection.length",
    "virulence",
    "group.size",
    "fertility",
    "modularity",
    "rnot",
    "transmission"
  )
  
  survival<-subset(gridim, survival==0.9)
  lowsurv<-subset(gridim, survival==0.6)
  survival$diff<-(survival$Extinction - lowsurv$Extinction) / 30  
  
  saveRDS(survival, path)
  
}else{
  survival <- readRDS(path)
}


######### mean difference in epidemic fadeout between networks with large and small group size #######

path <- here::here("data", "fig2_groupsize.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)
  
  tabpar <-
    expand.grid(
      unique(output[, "immunity.barrier"]),
      unique(output[, "infection.length"]),
      unique(output[, "virulence"]),
      unique(output[, "survival"]),
      unique(output[, "fertility"]),
      unique(output[, "modularity"]),
      unique(output[, "rnot"]),
      unique(output[, "transmission"])
    )
  colnames(tabpar) <-
    c(
      "immunity.barrier",
      "infection.length",
      "virulence",
      "survival",
      "fertility",
      "modularity",
      "rnot",
      "transmission" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
          immunity.barrier == tabpar[r, "immunity.barrier"] &
          modularity == tabpar[r, "modularity"] &
          virulence == tabpar[r, "virulence"] &
          infection.length == tabpar[r, "infection.length"] &
          fertility == tabpar[r, "fertility"] &
          survival == tabpar[r, "survival"] &
          rnot == tabpar[r, "rnot"] &
          transmission == tabpar[r, "transmission"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(group.size =  select[, "group.size"]), sum)
    outof <-
      aggregate(select[, "ext"], list(group.size =  select[, "group.size"]), length)
    persistence<-outof - extinction
    b<-matrix(c(extinction$x, persistence$x), nrow = 3, dimnames = list(group.size = extinction$group.size, Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
      gridim <-
        rbind(
          cbind(
            extinction,
            rep(Fisher$p.value, nrow(extinction)),
            rep(tabpar[r, "immunity.barrier"], nrow(extinction)),
            rep(tabpar[r, "infection.length"], nrow(extinction)),
            rep(tabpar[r, "virulence"], nrow(extinction)),
            rep(tabpar[r, "modularity"], nrow(extinction)), rep(tabpar[r, "survival"], nrow(extinction)), rep(
              tabpar[r, "fertility"], nrow(extinction)), rep(tabpar[r, "rnot"], nrow(extinction)), rep(tabpar[r, "transmission"], nrow(extinction))), gridim)
      }
    
  }
  colnames(gridim) = c(
    "group.size",
    "Extinction",
    "significant",
    "immunity.barrier",
    "infection.length",
    "virulence",
    "modularity",
    "survival",
    "fertility",
    "rnot",
    "transmission"
  )
  
  gpsize<-subset(gridim, group.size==100)
  smallgp<-subset(gridim, group.size==10)
  gpsize$diff<-(gpsize$Extinction - smallgp$Extinction) / 30
  
  saveRDS(gpsize, path)

}else{
  gpsize <- readRDS(path)
}


######### Infection length ###########################################################################

path <- here::here("data", "fig2_infection.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)
  
  tabpar <-
    expand.grid(
      unique(output[, "immunity.barrier"]),
      unique(output[, "group.size"]),
      unique(output[, "virulence"]),
      unique(output[, "survival"]),
      unique(output[, "fertility"]),
      unique(output[, "modularity"]),
      unique(output[, "rnot"]),
      unique(output[, "transmission"])
    )
  colnames(tabpar) <-
    c(
      "immunity.barrier",
      "group.size",
      "virulence",
      "survival",
      "fertility",
      "modularity",
      "rnot",
      "transmission" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
        immunity.barrier == tabpar[r, "immunity.barrier"] &
          modularity == tabpar[r, "modularity"] &
          virulence == tabpar[r, "virulence"] &
          group.size == tabpar[r, "group.size"] &
          fertility == tabpar[r, "fertility"] &
          survival == tabpar[r, "survival"] &
          rnot == tabpar[r, "rnot"] &
          transmission == tabpar[r, "transmission"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(infection.length =  select[, "infection.length"]), sum)
    outof <-
      aggregate(select[, "ext"], list(infection.length =  select[, "infection.length"]), length)
    persistence<-outof - extinction
    b<-matrix(c(extinction$x, persistence$x), nrow = 3, dimnames = list(infection.length = extinction$infection.length, Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
      gridim <-
        rbind(
          cbind(
            extinction,
            rep(Fisher$p.value, nrow(extinction)),
            rep(tabpar[r, "immunity.barrier"], nrow(extinction)),
            rep(tabpar[r, "group.size"], nrow(extinction)),
            rep(tabpar[r, "virulence"], nrow(extinction)),
            rep(tabpar[r, "modularity"], nrow(extinction)), rep(tabpar[r, "survival"], nrow(extinction)), rep(
              tabpar[r, "fertility"], nrow(extinction)), rep(tabpar[r, "rnot"], nrow(extinction)), rep(tabpar[r, "transmission"], nrow(extinction))), gridim)
    }
    
  }
  colnames(gridim)[c(1,2)] = c(
    "infection.length",
    "Extinction")
  
  infection<-subset(gridim, infection.length==30)
  lowinfec<-subset(gridim, infection.length==10)
  infection$diff<-(infection$Extinction - lowinfec$Extinction) / 30
  
  saveRDS(infection, path)

}else{
  infection <- readRDS(path)
}


######### high and low virulence #####################################################################

path <- here::here("data", "fig2_virulence.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)

  tabpar <-
    expand.grid(
      unique(output[, "immunity.barrier"]),
      unique(output[, "group.size"]),
      unique(output[, "infection.length"]),
      unique(output[, "survival"]),
      unique(output[, "fertility"]),
      unique(output[, "modularity"]),
      unique(output[, "rnot"]),
      unique(output[, "transmission"])
    )
  colnames(tabpar) <-
    c(
      "immunity.barrier",
      "group.size",
      "infection.length",
      "survival",
      "fertility",
      "modularity",
      "rnot",
      "transmission" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
        immunity.barrier == tabpar[r, "immunity.barrier"] &
          modularity == tabpar[r, "modularity"] &
          infection.length == tabpar[r, "infection.length"] &
          group.size == tabpar[r, "group.size"] &
          fertility == tabpar[r, "fertility"] &
          survival == tabpar[r, "survival"] &
          rnot == tabpar[r, "rnot"] &
          transmission == tabpar[r, "transmission"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(virulence =  select[, "virulence"]), sum)
    outof <-
      aggregate(select[, "ext"], list(virulence =  select[, "virulence"]), length)
    persistence<-outof - extinction
    b<-matrix(c(extinction$x, persistence$x), nrow = 3, dimnames = list(virulence = extinction$virulence, Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
      gridim <-
        rbind(cbind(
          extinction,
          rep(Fisher$p.value, nrow(extinction))),
          gridim)
    }
  }
  
  colnames(gridim)[c(1,2)] = c(
    "virulence",
    "Extinction")
  
  
  virulent<-subset(gridim, virulence==max(gridim$virulence))
  avirulent<-subset(gridim, virulence==min(gridim$virulence))
  virulent$diff<-(virulent$Extinction - avirulent$Extinction) / 30
  
  saveRDS(virulent, path)
  
}else{
  virulent <- readRDS(path)
}


######### fertility ##################################################################################

path <- here::here("data", "fig2_fertility.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)

  tabpar <-
    expand.grid(
      unique(output[, "immunity.barrier"]),
      unique(output[, "group.size"]),
      unique(output[, "infection.length"]),
      unique(output[, "survival"]),
      unique(output[, "virulence"]),
      unique(output[, "modularity"]),
      unique(output[, "rnot"]),
      unique(output[, "transmission"])
    )
  colnames(tabpar) <-
    c(
      "immunity.barrier",
      "group.size",
      "infection.length",
      "survival",
      "virulence",
      "modularity",
      "rnot",
      "transmission" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
        immunity.barrier == tabpar[r, "immunity.barrier"] &
          modularity == tabpar[r, "modularity"] &
          infection.length == tabpar[r, "infection.length"] &
          group.size == tabpar[r, "group.size"] &
          virulence == tabpar[r, "virulence"] &
          survival == tabpar[r, "survival"] &
          rnot == tabpar[r, "rnot"] &
          transmission == tabpar[r, "transmission"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(fertility =  select[, "fertility"]), sum)
    outof <-
      aggregate(select[, "ext"], list(fertility =  select[, "fertility"]), length)
    persistence<-outof - extinction
    b<-matrix(c(extinction$x, persistence$x), nrow = 3, dimnames = list(fertility = extinction$fertility, Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
      gridim <-
        rbind(cbind(
            extinction,
            rep(Fisher$p.value, nrow(extinction))),
            gridim)
    }
  }
  
  colnames(gridim)[c(1,2)] = c(
    "fertility",
    "Extinction")
  
  fertile<-subset(gridim, fertility==max(gridim$fertility))
  sterile<-subset(gridim, fertility==min(gridim$fertility))
  fertile$diff<-(fertile$Extinction - sterile$Extinction) / 30
  
  saveRDS(fertile, path)
  
}else{
  fertile <- readRDS(path)
}


######### transmission ###############################################################################

path <- here::here("data", "fig2_transmission.Rds")

if(!file.exists(path)){
  gridim <- data.frame(NULL)

  tabpar <-
    expand.grid(
      unique(output[, "immunity.barrier"]),
      unique(output[, "group.size"]),
      unique(output[, "infection.length"]),
      unique(output[, "survival"]),
      unique(output[, "fertility"]),
      unique(output[, "modularity"]),
      unique(output[, "rnot"]),
      unique(output[, "virulence"])
    )
  colnames(tabpar) <-
    c(
      "immunity.barrier",
      "group.size",
      "infection.length",
      "survival",
      "fertility",
      "modularity",
      "rnot",
      "virulence" 
    )
  
  for (r in 1:nrow(tabpar)){
    select <-
      subset(
        output,
        immunity.barrier == tabpar[r, "immunity.barrier"] &
          modularity == tabpar[r, "modularity"] &
          infection.length == tabpar[r, "infection.length"] &
          group.size == tabpar[r, "group.size"] &
          virulence == tabpar[r, "virulence"] &
          survival == tabpar[r, "survival"] &
          rnot == tabpar[r, "rnot"] &
          fertility == tabpar[r, "fertility"]
      )
    
    extinction <-
      aggregate(select[, "ext"], list(transmission =  select[, "transmission"]), sum)
    outof <- c(30,30)
    persistence<- outof - extinction$x
    b<-matrix(c(extinction$x, persistence), nrow = 2, dimnames = list(transmission = extinction$transmission, Persistence = c("NO", "YES")))
    Fisher<-fisher.test(b)
    if(Fisher$p.value < 0.05){  
      gridim <-
        rbind(cbind(
          extinction,
          rep(Fisher$p.value, nrow(extinction))),
          gridim)
    }
  }
  colnames(gridim)[c(1,2)] = c(
    "transmission",
    "Extinction")

  frequency<-subset(gridim, transmission=="Frequency")
  ddp<-subset(gridim, transmission=="Density")
  frequency$diff<-(frequency$Extinction - ddp$Extinction) / 30
  
  saveRDS(frequency, path)
  
}else{
  frequency <- readRDS(path)
}


######################################################################################################
######### ALL PARAMETERS #############################################################################
######################################################################################################

path <- here::here("data", "fig2_allvars.Rds")

if(!file.exists(path)){
  fulldata<-rbind(cbind(modularity$diff, rep("Contact ratio", length(modularity$diff))), cbind(gpsize$diff, rep("Group size", length(gpsize$diff))), cbind(survival$diff, rep("Survival", length(survival$diff))),   cbind(fertile$diff, rep("Fertility", length(fertile$diff))), cbind(infection$diff, rep("Infection length", length(infection$diff))),
                  cbind(virulent$diff, rep("Virulence", length(virulent$diff))), cbind(frequency$diff, rep("Transmission type", length(frequency$diff))), 
                  cbind(immunity$diff, rep("Age at first between-group contact",length(immunity$diff))))
  
  colnames(fulldata)<-c("diff", "parameter")
  fulldata<-as.data.frame(fulldata)
  fulldata$diff<-as.numeric(as.character(fulldata$diff))
  fulldata$parameter<-as.character(fulldata$parameter)
  
  fulldata$type<-ifelse(fulldata$parameter %in% c("Infection length", "Virulence", "Transmission type"), "viral", "host")
  
  saveRDS(fulldata, path)
  
}else{
  fulldata <- readRDS(path)
}


######################################################################################################
## PLOT ##############################################################################################
######################################################################################################

theme_lucile <-
  theme_classic(base_size = 18, base_family = "Times New Roman") +
  theme(axis.title.x = element_text(size = 17, margin = margin(t = 15)),
        axis.title.y = element_text(size = 17, margin = margin(r = 8)),
        axis.text = element_text(size = 14),
        axis.ticks = element_line(size = .8), 
        axis.ticks.length = unit(.2, "cm"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.3),
        strip.background = element_rect(colour = "black", fill = NA, size = 1.3),
        strip.text = element_text(size = 17, face = "bold", margin = margin(10, 0, 10, 0)),
        panel.spacing.x = unit(1.3, "lines"),
        panel.spacing.y = unit(.8, "lines"),
        plot.title = element_text(face = "bold", size = 24, margin = margin(0, 0, 10, 0)),
        plot.tag = element_text(size = 24))

p_viral <- 
  fulldata %>%  
  filter(type == "viral") %>% 
  ggplot(aes(diff)) +
  geom_histogram(
    fill = "#efa253",
    color = colorspace::darken("#efa253", .4),
    binwidth = .25,
    size = 1.1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40",
    size = .8
  ) +
  scale_x_continuous(
    expand = c(.025, .025),
    limits = c(-1.5, 2),
    breaks = seq(-1.5, 2, by = 0.5)
  ) +
  scale_y_continuous(
    expand = c(.035, .035),
    limits = c(0, 1100),
    breaks = seq(0, 1000 , by = 250)
  ) +
  facet_wrap(~parameter) +
  labs(
    x = NULL,
    y = "Count",
    title = "Virus-related traits",
    tag = expression(paste("(", italic("a"), ")"))
  ) +
  theme_lucile

p_host <- 
  fulldata %>%
  filter(type == "host") %>% 
  ggplot(aes(diff)) +
  geom_histogram(
    fill = "#1a8775",
    color = colorspace::darken("#1a8775", .4),
    binwidth = .25,
    size = 1.1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40",
    size = .8
  ) +
  scale_x_continuous(
    expand = c(.025, .025),
    limits = c(-1.5, 2),
    breaks = seq(-1.5, 2, by = 0.5)
  ) +
  scale_y_continuous(
    expand = c(.035, .035),
    limits = c(0, 550),
    breaks = seq(0, 500, by = 100)
  ) +
  facet_wrap(~parameter) +
  labs(
    x = "Difference in fade-out probabilities",
    y = "Count",
    title = "Host-related traits",
    tag = expression(paste("(", italic("b"), ")"))
  ) +
  theme_lucile

p <-
  p_viral / p_host + 
  plot_layout(heights = c(.43, 1)) 

ggsave(here::here("figures", "fig_2_grey.pdf"), width = 16, height = 11, device = cairo_pdf)  


## red line
p_viral <- 
  fulldata %>%  
  filter(type == "viral") %>% 
  ggplot(aes(diff)) +
    geom_histogram(
      fill = "#efa253",
      color = colorspace::darken("#efa253", .4),
      binwidth = .25,
      size = 1.1
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#a00101", 
      size = .8
    ) +
    scale_x_continuous(
      expand = c(.025, .025),
      limits = c(-1.5, 2),
      breaks = seq(-1.5, 2, by = 0.5)
    ) +
    scale_y_continuous(
      expand = c(.035, .035),
      limits = c(0, 1100),
      breaks = seq(0, 1000 , by = 250)
    ) +
    facet_wrap(~parameter) +
    labs(
      x = NULL,
      y = "Count",
      title = "Virus-related traits",
      tag = expression(paste("(", italic("a"), ")"))
    ) +
    theme_lucile

p_host <- 
  fulldata %>%
  filter(type == "host") %>% 
  ggplot(aes(diff)) +
    geom_histogram(
      fill = "#1a8775",
      color = colorspace::darken("#1a8775", .4),
      binwidth = .25,
      size = 1.1
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#a00101",
      size = .8
    ) +
    scale_x_continuous(
      expand = c(.025, .025),
      limits = c(-1.5, 2),
      breaks = seq(-1.5, 2, by = 0.5)
    ) +
    scale_y_continuous(
      expand = c(.035, .035),
      limits = c(0, 550),
      breaks = seq(0, 500, by = 100)
    ) +
    facet_wrap(~parameter) +
    labs(
      x = "Difference in fade-out probabilities",
      y = "Count",
      title = "Host-related traits",
      tag =  expression(paste("(", italic("b"), ")"))
    ) +
    theme_lucile

p <-
  p_viral / p_host + 
  plot_layout(heights = c(.43, 1)) 

ggsave(here::here("figures", "fig_2_red.pdf"), width = 16, height = 11, device = cairo_pdf)


### alternative: all labels, smaller
p_viral <- 
  fulldata %>%  
  filter(type == "viral") %>% 
  ggplot(aes(diff)) +
    geom_histogram(
      fill = "#1a8775",
      color = colorspace::darken("#1a8775", .4),
      binwidth = .25,
      size = 1.1
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40", #efa253
      size = .8
    ) +
    scale_x_continuous(
      expand = c(.025, .025),
      limits = c(-1.5, 2),
      breaks = seq(-1.5, 2, by = 0.25),
      labels = c("-1.5", "-1.25", "-1.0", "-0.25", "-0.5", "-0.25", "0.0", "0.25", "0.5", 
                 "0.75", "1.0", "1.25", "1.5", "1.75", "2.0")
    ) +
    scale_y_continuous(
      expand = c(.025, .025),
      limits = c(0, 1100),
      breaks = seq(0, 1000 , by = 250)
    ) +
    facet_wrap(~parameter) +
    labs(
      x = NULL,
      y = "Count",
      title = "Virus-related traits",
      tag = expression(paste("(", italic("a"), ")"))
    ) +
    theme_lucile +
    theme(axis.text.x = element_text(size = 12))

p_host <- 
  fulldata %>%
  filter(type == "host") %>% 
  ggplot(aes(diff)) +
    geom_histogram(
      fill = "#efa253",
      color = colorspace::darken("#efa253", .4),
      binwidth = .25,
      size = 1.1
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40", #efa253
      size = .8
    ) +
    scale_x_continuous(
      expand = c(.025, .025),
      limits = c(-1.5, 2),
      breaks = seq(-1.5, 2, by = 0.25),
      labels = c("-1.5", "-1.25", "-1.0", "-0.25", "-0.5", "-0.25", "0.0", "0.25", "0.5", 
                 "0.75", "1.0", "1.25", "1.5", "1.75", "2.0")
    ) +
    scale_y_continuous(
      expand = c(.025, .025),
      limits = c(0, 550),
      breaks = seq(0, 500, by = 100)
    ) +
    facet_wrap(~parameter) +
    labs(
      x = "Difference in fade-out probabilities",
      y = "Count",
      title = "Host-related traits",
      tag = expression(paste("(", italic("b"), ")"))
    ) +
    theme_lucile +
    theme(axis.text.x = element_text(size = 12))

p <-
  p_viral / p_host + 
  plot_layout(heights = c(.43, 1)) 

ggsave(here::here("figures", "fig_2_alt.pdf"), width = 19, height = 12, device = cairo_pdf)

