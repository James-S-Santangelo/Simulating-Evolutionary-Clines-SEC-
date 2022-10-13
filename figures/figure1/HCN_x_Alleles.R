
library(plotly)
library(reshape2)

setwd('/Users/jamessantangelo/Documents/Academia/Doctorate_PhD/Projects/SEC_Simulating.evolutionary.clines/SEC_Visuals/HCN-x-Alleles plot')
dat <- read.csv('Allele_Freq_HCN.csv', header = T)
Cyan <- acast(dat, dat$pA ~ dat$pB, value.var="Cyan")

p <- plot_ly(x = dat$pB, y = dat$pB, z = ~Cyan, colors=c("#EBEBEB","#0A0A0A")) %>% 
  add_contour(contours = list(
    end = 1.0, 
    size = 0.05, 
    start = 0
  )) %>%
  layout(
    scene = list(
      xaxis = list(title = ""),
      yaxis = list(title = ""),
      zaxis = list(title = "")
    ))
p 

htmlwidgets::saveWidget(as_widget(p), "HCN_x_Alleles_GreyContour.html")

