# Analysis code to accompany the publication: 
# "Domestication impacts on plant-herbivore interactions: a meta-analysis"
# S.R. Whitehead, M.M. Turcotte, and K. Poveda
# Philosophical Transactions of the Royal Society B
# DOI: http://dx.doi.org/10.1098/rstb.2016.0034


#Phylogenetic Meta-Analysis of the effects of crop domestication 
#on herbivore resistance and plant defence traits

#All analyses conducted with R v. 3.2.2 -- "Fire Safety"


################################ 
#### Script Description 
#################################################

#1 Load database (in csv) and necessary libraries
#2 Load phylogeny  - prune species not in this analysis
#3 Calculate effect sizes + variances - Modify them according to direction of effect
#4- ANALYSIS--Herbivore resistance data
      #4a - Overall random effects model
      #4b - Models with moderators

#5- ANALYSIS--Plant trait data
    #5a - Overall random effects model
    #5b - Models with moderators
#6 Testing for Publication Bias 
#7 Pretty Figs for Publication 




##############################################################
#1 Load database (in csv) and libraries
##############################################################


dat <- read.csv(file="Whitehead_et_al_PTRSB_database.csv",stringsAsFactors=F)

library(ape)  #version 3.3
library(phytools) #version 0.5-20
library(metafor) #version 1.9-8
library(geiger) #version 2.0.6
library(ggplot2) #version 1.0.1
library(multcomp) #version 1.4-1
library(plyr) #version 1.8.3
library(RColorBrewer) #version 1.1-2


##############################################################
#2 Load phylogeny  - prune species not in this analysis
##############################################################

# Load in tree with all crop species
full.tree<- read.newick(file="Whitehead_et_al_PTRSB_phylogeny.new")
full.tree<- collapse.singles(full.tree) # remove single nodes
plot(full.tree, cex=0.8)

# Prune tree to spp in dataset
tree <- drop.tip(full.tree, full.tree$tip.label[ -match(unique(dat$D.Spp), full.tree$tip.label) ])


# Checking tree
    #ultrametric?
      is.ultrametric(tree)
    # do species names match between tree and data?
       treedata(tree, setNames(unique(dat$D.Spp),unique(dat$D.Spp)), warnings=T) 


######################################################  
#3 Calculate effect sizes + variances - Modify them according to direction of effect
###################################################### 

# Calculate Hedges' d and variances 

dat <-   escalc(measure= 'SMD',m1i=Mean.Dom, m2i=Mean.Wild, 
         sd1i=SD.Dom, sd2i=SD.Wild, n1i=N.Dom, n2i=N.Wild, data=dat)   
# This adds 2 columns to the dataset:  yi = Hedges d and vi=variance

# Correct direction of Hedges' d so that negative values always indicate an 
   # decrease in resistance during domestication 
   # Response.Direction variable indicates how trait change is related to resistance 
   #   0=indicates that the more resistant plant has a higher mean value (e.g. secondary metabolite concentration)
   #   1=indicates that the least resistant plant has a higher mean value (e.g. insect performance)
   # Adjust direction of change:  those with 1 values need to be switched

for (i in 1:length(dat$yi)){
  if (dat$Response.Direction[i] == 1) {
    dat$yi[i] <- -1*dat$yi[i]
  }
}



#######################################################
##   #4 Analysis--herb data
###############################################################

#  Subset the data to include only herbivore resistance measures
h.dat <- dat[dat$Response.Type %in% c('herb','herb.perform','herb.pref'),  ]

###############
# 4a  Random effects only model

temp.h.dat <- h.dat
temp.tree <- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ])
m.h.0 <- rma.mv(yi, vi, 
               random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
               R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
               data=temp.h.dat, verbose = T, method='REML')



###################
# 4b Models with moderators


### Response Type ################

temp.h.dat <- h.dat
temp.tree <- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ])

m.h.rt <- rma.mv(yi, vi, 
                mods= ~factor(Response.Type),
                random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                data=temp.h.dat, verbose = T, method='REML')

#Now fit with no intercept to get the final parameter estimates
m.h.rt.noInt <- rma.mv(yi, vi, 
                mods= ~factor(Response.Type)-1,
                random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                data=temp.h.dat, verbose = T, method='REML')

m.h.rt.noInt
summary(glht(m.h.rt.noInt, linfct=contrMat(c("herb"=1, "herb.perform"=1,"herb.pref"=1), type="Tukey")))



### Organ Measured ###############

#taking only lvs and seeds b/c these are the only organs with N (sp*study) > 3
temp.h.dat <- h.dat[h.dat$Org.Meas %in% c('l','s'), ] 
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ])

m.h.om <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Meas),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')


#Now fit with no intercept to get the final parameter estimates
m.h.om.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Meas)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')


m.h.om.noInt
summary(glht(m.h.om.noInt, linfct=contrMat(c("l"=1,"s"=1), type="Tukey")))




### Crop Use ################

temp.h.dat <- h.dat
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ])

m.h.cu <- rma.mv(yi, vi, 
                   mods= ~factor(Crop.Use),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')


#Now fit with no intercept to get the final parameter estimates
m.h.cu.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Crop.Use)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')


m.h.cu.noInt
summary(glht(m.h.cu.noInt, linfct=contrMat(c("fr"=1,"g"=1 "l"=1,"nf"=1,"o"=1, "v"=1), type="Tukey")))



### Domestication Extent (= landrace vs modern cultivar) #########################

#removing data where domestication extent was not reported
temp.h.dat<- h.dat[h.dat$Dom.Extent != "nr", ]
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ])


m.h.de<- rma.mv(yi, vi, 
                        mods= ~factor(Dom.Extent),     
                        random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                        R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                        data=temp.h.dat, verbose = T, method='REML')



#Now fit with no intercept to get the final parameter estimates
m.h.de.noInt <- rma.mv(yi, vi, 
                        mods= ~factor(Dom.Extent)-1,     
                        random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                        R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                        data=temp.h.dat, verbose = T, method='REML')


m.h.de.noInt
summary(glht(m.h.de.noInt, linfct=contrMat(c("l"=1, "m"=1), type="Tukey")))



### Life History #################

#Removing data where wild and domesticated differ in life history
temp.h.dat <- h.dat[h.dat$Life.Hist != "diff", ]
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ])

m.h.lh <- rma.mv(yi, vi, 
                   mods= ~factor(Life.Hist),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')



#Now fit with no intercept to get the final parameter estimates
m.h.lh.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Life.Hist)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')


m.h.lh.noInt
summary(glht(m.h.lh.noInt, linfct=contrMat(c("ha"=1, "hp"=1,"wp"=1), type="Tukey")))



### Organ Use (= harvested or not harvested) ##########################

temp.h.dat <- h.dat
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ])

m.h.ou <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Use),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')


#Now fit with no intercept to get the final parameter estimates
m.h.ou.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Use)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')


m.h.ou.noInt
summary(glht(m.h.ou.noInt, linfct=contrMat(c("h"=1,"nh"=1), type="Tukey")))





#######################################################
##   #5 Analysis--plant data
###############################################################

#  Subset the data to include only plant trait measures
p.dat <- dat[dat$Response.Type %in% c('chem','phys'),  ]

####################
# 5a  Random effects only model

temp.p.dat<- p.dat
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])

m.p.0 <- rma.mv(yi, vi, 
               random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
               R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
               data=temp.p.dat, verbose = T, method='REML')


###################
# 5b Models with moderators


## Response Type ###############

temp.p.dat<- p.dat
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])
m.p.rt <- rma.mv(yi, vi, 
                mods= ~factor(Response.Type),
                random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                data=temp.p.dat, verbose = T, method='REML')



#Now fit with no intercept to get the final parameter estimates
m.p.rt.noInt <- rma.mv(yi, vi, 
                mods= ~factor(Response.Type)-1,
                random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                data=temp.p.dat, verbose = T, method='REML')

m.p.rt.noInt
summary(glht(m.p.rt.noInt, linfct=contrMat(c("chem"=1, "phys"=1), type="Tukey")))




## Organ Measured #################

#taking only frts, lvs, seeds, and other vegetative b/c these are the only organs with N (sp*study) > 3
temp.p.dat <- p.dat[p.dat$Org.Meas %in% c('l','s', 'v', 'fr'), ] 
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])

m.p.om <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Meas),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


#Now fit with no intercept to get the final parameter estimates
m.p.om.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Meas)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


m.p.om.noInt
summary(glht(m.p.om.noInt, linfct=contrMat(c("fr"=1,"l"=1,"s"=1, "v"=1), type="Tukey")))



## Domestication Extent ####################

#removing data where domestication extent was not reported
temp.p.dat<- p.dat[p.dat$Dom.Extent != "nr", ]
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])


m.p.de<- rma.mv(yi, vi, 
                        mods= ~factor(Dom.Extent),     
                        random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                        R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                        data=temp.p.dat, verbose = T, method='REML')


#Now fit with no intercept to get the final parameter estimates
m.p.de.noInt <- rma.mv(yi, vi, 
                        mods= ~factor(Dom.Extent)-1,     
                        random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                        R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                        data=temp.p.dat, verbose = T, method='REML')


m.p.de.noInt
summary(glht(m.p.de.noInt, linfct=contrMat(c("l"=1, "m"=1), type="Tukey")))



## Life History ########################

#Removing data where wild and domesticated differ in life history
temp.p.dat <- p.dat[p.dat$Life.Hist != "diff", ]
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])


m.p.lh <- rma.mv(yi, vi, 
                   mods= ~factor(Life.Hist),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


#Now fit with no intercept to get the final parameter estimates
m.p.lh.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Life.Hist)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


m.p.lh.noInt
summary(glht(m.p.lh.noInt, linfct=contrMat(c("ha"=1, "hp"=1,"wp"=1), type="Tukey")))



## Crop Use #######################

#taking only frts, lvs, non-food, and vegetables b/c these are the only groups with N (sp*study) > 3
temp.p.dat <- p.dat[p.dat$Crop.Use %in% c('fr','l','nf', 'v'),  ]
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])


m.p.cu <- rma.mv(yi, vi, 
                   mods= ~factor(Crop.Use),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


#Now fit with no intercept to get the final parameter estimates
m.p.cu.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Crop.Use)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


m.p.cu.noInt
summary(glht(m.p.cu.noInt, linfct=contrMat(c("fr"=1,"g"=1, "l"=1,"nf"=1,"o"=1, "v"=1), type="Tukey")))


## Organ Use (= harvested or not harvested) ####################

temp.p.dat <- p.dat
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])


m.p.ou <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Use),     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')



#Now fit with no intercept to get the final parameter estimates
m.p.ou.noInt <- rma.mv(yi, vi, 
                   mods= ~factor(Org.Use)-1,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


m.p.ou.noInt
summary(glht(m.p.ou.noInt, linfct=contrMat(c("n"=1,"y"=1), type="Tukey")))






#######################################################################################
### #6 Testing for Publication Bias
#######################################################################


##  Funnel Plots ##################

tiff("Meta_funnel_revised.tiff", compression="lzw", width=183, height=200, units="mm", res=600)
   #use to make tiff file for pasting into supplementary materials word doc

par(mfrow=c(1,2))

funnel(m.h.0, xlab=expression("Effect Size ("~italic("d")~")"))
par(font=2)
mtext("(a)", side=3, adj=0, padj=-1)

funnel(m.p.0, xlab=expression("Effect Size ("~italic("d")~")"))
par(font=2)
mtext("(b)", side=3, adj=0, padj=-1)


dev.off()





## Study Year vs Effect Size  ##########################


temp.h.dat <- h.dat
temp.h.dat$Study.Year[temp.h.dat$Study.Year == 'in prep'] <- 2016
temp.h.dat$Study.Year <- as.numeric(substr(temp.h.dat$Study.Year, 1, 4))
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.h.dat$D.Spp), tree$tip.label) ]) 
m.h.sy <- rma.mv(yi, vi, 
                   mods= ~Study.Year,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.h.dat, verbose = T, method='REML')



temp.p.dat <- p.dat
temp.p.dat$Study.Year[temp.p.dat$Study.Year == 'in prep'] <- 2016
temp.p.dat$Study.Year <- as.numeric(substr(temp.p.dat$Study.Year, 1, 4))
temp.tree<- drop.tip(tree, tree$tip.label[ -match(unique(temp.p.dat$D.Spp), tree$tip.label) ])
m.p.sy <- rma.mv(yi, vi, 
                   mods= ~Study.Year,     
                   random = list(~ 1 | D.Spp/D.Var, ~1| D.Spp/W.Spp/W.Var, ~ 1 | factor(Study)),
                   R= list(D.Spp = vcv(temp.tree,corr=T, model=Brownian)),  
                   data=temp.p.dat, verbose = T, method='REML')


tiff("Meta_sy.tiff", compression="lzw", width=183, height=90, units="mm", res=300)
   #use to make tiff file for pasting into supplementary materials word doc

par(mfrow=c(1,2))

#For herbivore data
seq <- seq(1975,2016, 0.1)
preds.h <- predict(m.h.sy, newmods=seq)
plot(h.dat$Study.Year, h.dat$yi, pch = 1, cex = 0.5,
        xlab = "Study Year", ylab = expression("Effect Size ("~italic("d")~")"), 
        las = 1, bty = "l", cex.axis=1, cex.lab=1) 
lines(seq, preds.h$pred) 
lines(seq, preds.h$ci.lb, lty = "dashed") 
lines(seq, preds.h$ci.ub, lty = "dashed") 
#abline(h = 0, lty = "dotted")
mtext("(a)", side=3, adj=0, padj=-1, cex=1)
text(1987, 20, expression("Q"["M"]~"=0.50;"~italic("P")~"=0.48"), cex=0.8)


#For plant data
seq <- seq(1970,2016, 0.1)
preds.p <- predict(m.p.sy, newmods=seq)
plot(p.dat$Study.Year, p.dat$yi, pch = 1, cex = 0.5,
        xlab = "Study Year", ylab = expression("Effect Size ("~italic("d")~")"), 
        las = 1, bty = "l", cex.axis=1, cex.lab=1) 
lines(seq, preds.p$pred) 
lines(seq, preds.p$ci.lb, lty = "dashed") 
lines(seq, preds.p$ci.ub, lty = "dashed") 
#abline(h = 0, lty = "dotted")
mtext("(b)", side=3, adj=0, padj=-1, cex=1)
text(1982, 50, expression("Q"["M"]~"=0.20;"~italic("P")~"=0.66"), cex=0.8)

dev.off()







#######################################################################################
### #7 Pretty Figs for Publication
#######################################################################



####  Figure 1: Phylogeny 
########################


##First estimate the mean domestication effect sizes for each species individually

 #For herbivore resistance
  spp <- unique(h.dat$D.Spp)
  
  #results table
  h.res <-  data.frame(matrix(nrow=0, ncol=7, NA))
  names(h.res)<- c('D.Spp','N.Study','estimate','se','pval','ci.lb','ci.ub')
  
  # loop by species. 
  for (i in 1:length(spp)) {
    temp.dat<-  h.dat[ h.dat$D.Spp == spp[i] ,  ]
    
    if(dim(temp.dat)[1] >1){  # do this only if you have more than 1 data point
      
      meta<- rma.mv(yi, vi,   # model with no moderators or impact of domesticated spp
                    random = list(~ 1 | D.Var, ~1| W.Spp/W.Var, ~ 1 | factor(Study)),
                    data=temp.dat, verbose = F, method='REML')
      
      h.res[i, ] <- c(spp[i], meta$s.nlevels[4], meta[[1]], meta[[2]], meta[[4]],meta[[5]], meta[[6]]) # extra results and add to results table
    }else{  # species only have 1 row and thus we can't do a meta-anlysis. For these just get their yi value. 
      
      h.res[i, ]<- c(spp[i], 1, temp.dat$yi, 'NA' , 'NA' ,'NA' ,'NA' )
    }
  }
  
h.res
  

 #For plant traits
  spp <- unique(p.dat$D.Spp)
  
  #results table
  p.res <-  data.frame(matrix(nrow=0, ncol=7, NA))
  names(p.res)<- c('D.Spp','N.Study','estimate','se','pval','ci.lb','ci.ub')
  
  # loop by species. 
  for (i in 1:length(spp)) {
    temp.dat<-  p.dat[ p.dat$D.Spp == spp[i] ,  ]

   if (temp.dat$D.Spp[1] == 'Allium_cepa') { #onion is weird b/c there are no unique varieties nested within species, so we need a unique model

        meta<- rma.mv(yi, vi,   # model with no moderators or impact of domesticated spp
                    random = list(~ 1 | D.Var, ~1| W.Spp, ~ 1 | factor(Study)),
                    data=temp.dat, verbose = F, method='REML')

        p.res[i, ] <- c(spp[i], 2, meta[[1]], meta[[2]], meta[[4]],meta[[5]], meta[[6]]) # extra results and add to results table
      }else
        if(dim(temp.dat)[1] >1){  # do this only if you have more than 1 data point
      
             meta<- rma.mv(yi, vi,   # model with no moderators or impact of domesticated spp
                    random = list(~ 1 | D.Var, ~1| W.Spp/W.Var, ~ 1 | factor(Study)),
                    data=temp.dat, verbose = F, method='REML')
      
             p.res[i, ] <- c(spp[i], meta$s.nlevels[4], meta[[1]], meta[[2]], meta[[4]],meta[[5]], meta[[6]]) # extra results and add to results table
          }else{  # species only have 1 row and thus we can't do a meta-anlysis. For these just get their yi value. 
      
             p.res[i, ]<- c(spp[i], 1, temp.dat$yi, 'NA' , 'NA' ,'NA' ,'NA' )
    }
  }
  
p.res



## Set up results for each species as variables with names so they can easily be called in the plotting
h.estimate <- as.numeric(h.res$estimate)
names(h.estimate) <- h.res$D.Spp

p.estimate <- as.numeric(p.res$estimate)
names(p.estimate) <- p.res$D.Spp

h.se <- as.numeric(h.res$se)
names(h.se) <- h.res$D.Spp

p.se <- as.numeric(p.res$se)
names(p.se) <- p.res$D.Spp

h.ci.ub <- as.numeric(h.res$ci.ub)
names(h.ci.ub) <- h.res$D.Spp

h.ci.lb <- as.numeric(h.res$ci.lb)
names(h.ci.lb) <- h.res$D.Spp

p.ci.ub <- as.numeric(p.res$ci.ub)
names(p.ci.ub) <- p.res$D.Spp

p.ci.lb <- as.numeric(p.res$ci.lb)
names(p.ci.lb) <- p.res$D.Spp




## Set up colors for each family
 

#Make a color palette
display.brewer.all()
colourCount = length(unique(dat$Family))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
pal <- getPalette(colourCount)

n=seq(1:26)
barplot(n, col=pal)


#export list of species in order of phylogeny tips
write.csv(tree$tip.label, file="tips.csv")
  #I then added a column to this list where I manually assigned colors from the "pal" list to each 
  #species according to family and making sure that similar colors would not be adjacent to each other 


##Read back in modified list of colors
colors <- read.csv("spfamcol.csv")
c <- colors$Cols
names(c) <- colors$D.Spp


#Plot tree with color circles on tips for each family
plotTree(tree,fsize=0.8,ftype="i")
tiplabels(pch=21, bg=as.vector(c), cex=1)   




##  The plot  ##########


dev.new(width=183, height=235, units="mm")
## create a split plot
layout(matrix(c(1,2,3),1,3),c(0.6,0.3, 0.3))
par(oma=c(0,0.1,0,0.2)) ##c(bottom, left, top, right)
par(lwd=0.8)

#plot tree
plotTree(tree,fsize=0.82,ftype="i", offset=1, mar=c(4, 0, 2, 0), lwd=0.8)
tiplabels(pch=21, bg=as.vector(c), cex=1.6, adj=8) 
#add.scale.bar() 


## add bar plots

h.plot <- barplot(h.estimate[tree$tip.label],horiz=TRUE,width=1,space=0,
           ylim=c(1,length(tree$tip.label))-0.5,names="", xlim=c(-11, 9), 
           col=as.vector(c), xlab="", axes=FALSE, xpd=FALSE, lwd = 0.8)
par(font=1)  #1=normal, 2=bold
lines(x=c(0,0), y = c(0,73), lty=3)
axis(1,at=c(-10,-5,0,5),pos=-0.5, cex.axis=1)
mtext(expression("Effect size ("~italic("d")~")"), 
       side=1,adj=0.4, padj=1.1, cex=0.7)
mtext("Herbivore Resistance", side=3, adj=0.5, padj=1, cex=0.7)
segments(h.ci.lb[tree$tip.label], h.plot, 
         h.ci.ub[tree$tip.label],  h.plot, lwd = 0.8)

#add arrows where CIs are out of range
arrows(x0=-10.7, y0=10.5, x1=-11, length=0.05, lwd = 0.8)
arrows(x0=-10.7, y0=15.5, x1=-11, length=0.05, lwd = 0.8)
arrows(x0=-10.7, y0=65.5, x1=-11, length=0.05, lwd = 0.8)


p.plot <- barplot(p.estimate[tree$tip.label],horiz=TRUE,width=1,space=0,
          ylim=c(1,length(tree$tip.label))-0.5,names="", xlim=c(-10, 10), 
          col=as.vector(c), xlab="", axes=FALSE, xpd=TRUE, lwd = 0.8)
lines(x=c(0,0), y = c(0,73), lty=3)
axis(1,at=c(-10,-5,0,5),pos=-0.5, cex.axis=1)
par(font=1)  #1=normal, 2=bold
mtext(expression("Effect size ("~italic("d")~")"), 
       side=1,cex=0.7, adj=0.4, padj=1.1)
mtext("Plant Traits", side=3, adj=0.5, padj=1, cex=0.7)
segments(p.ci.lb[tree$tip.label] , p.plot, 
         p.ci.ub[tree$tip.label] ,  p.plot, lwd = 0.8)

#add arrows where CIs are out of range
arrows(x0=-7, y0=39, x1=-9.5, length=0.05, col="#F7864E", lwd = 0.8) #For Fragaria bar
text(-8, 38, "-32", cex=0.8) #for Fragaria
arrows(x0=-9.7, y0=40.5, x1=-10, length=0.05, lwd = 0.8) #Fragaria, LB
arrows(x0=-9.7, y0=4.5, x1=-10, length=0.05, lwd = 0.8) #Helianthus, LB
arrows(x0=-9.7, y0=42.5, x1=-10, length=0.05, lwd = 0.8) #Rubus, LB
arrows(x0=9.7, y0=17.5, x1=10, length=0.05, lwd = 0.8) #Coffea, UB
arrows(x0=9.7, y0=40.5, x1=10, length=0.05, lwd = 0.8) #Fragaria, UB
arrows(x0=9.7, y0=61.5, x1=10, length=0.05, lwd = 0.8)  #Allium, UB


dev.copy2pdf(file="Meta_phylo_revised.pdf") 








##  Figure 2  #############################
#######################

#big forest plot showing effect size estimates for each model


### First for herbivore resistance data
### collect model estimates and corresponding  confidence intervals
  
  # input the models for the figure here
  h.mod1 <- m.h.0   #Overall
  h.mod2 <- m.h.rt.noInt  #Response Type
  h.mod3 <- m.h.om.noInt  #Organ Measured
  h.mod4 <- m.h.ou.noInt  #Organ Use
  h.mod5 <- m.h.de.noInt  #Domestication Extent
  h.mod6 <- m.h.cu.noInt  #Crop use
  h.mod7 <- m.h.lh.noInt  #Life History
  
 
  # re-order coefficients if necessary
  est.cu <- c(h.mod6$b[2,], h.mod6$b[3,], h.mod6$b[5,], h.mod6$b[1,], h.mod6$b[6,], h.mod6$b[4,])
  ci.lb.cu <- c(h.mod6$ci.lb[2], h.mod6$ci.lb[3], h.mod6$ci.lb[5], h.mod6$ci.lb[1], h.mod6$ci.lb[6], h.mod6$ci.lb[4])
  ci.ub.cu <- c(h.mod6$ci.ub[2], h.mod6$ci.ub[3], h.mod6$ci.ub[5], h.mod6$ci.ub[1], h.mod6$ci.ub[6], h.mod6$ci.ub[4])

  h.estimates <- c(h.mod1[[1]], h.mod2[[1]], h.mod3[[1]], h.mod4[[1]], h.mod5[[1]], est.cu, h.mod7[[1]])
  h.ci.lb  <- c(h.mod1[[5]], h.mod2[[5]], h.mod3[[5]], h.mod4[[5]], h.mod5[[5]], ci.lb.cu, h.mod7[[5]])
  h.ci.ub <- c(h.mod1[[6]], h.mod2[[6]], h.mod3[[6]], h.mod4[[6]], h.mod5[[6]], ci.ub.cu, h.mod7[[6]])
 
  
  # create vector with labels
  h.labels <- c("", "Herbivory", "Performance", "Preference", "Leaves", "Seeds", 
                "Non-harvested","Harvested", "Landrace", "Modern", "Grains", "Legumes", "Oilseed", "Fruits", 
                "Vegetables", "Non-food", "Herb Annual", "Herb Perennial", "Woody Perennial")
  h.rows <- c(28,  25.5:23.5,  21:20,  17.5:16.5,  14:13,  10.5:5.5,  3:1)

  # create vector of sample sizes (see supplementary materials)
  h.N <- c(30, 17, 21, 10, 15, 7, 21, 13, 6, 22, 4, 6, 4, 8, 8, 3, 21, 4, 7)
  names(h.N) <- c("Overall", "Herb", "Perf", "Pref", "Lvs", "Seeds", "Non-harv", "Harv", "Landrace",
      "Modern", "Grains", "Legumes", "Oilseed", "Fruits", "Veg", "Non-food", "Herb Ann", "Herb Per", "Woody Per")



### Now for plant data
### collect model estimates and corresponding  confidence intervals

  # input the models for figure here
  p.mod1 <- m.p.0    #Overall
  p.mod2 <- m.p.rt.noInt   #Response Type
  p.mod3 <- m.p.om.noInt   #Organ Measured
  p.mod4 <- m.p.ou.noInt   #Organ Use
  p.mod5 <- m.p.de.noInt   #Domestication Extent
  p.mod6 <- m.p.cu.noInt   #Crop use
  p.mod7 <- m.p.lh.noInt   #Life History
  
 ### collect model estimates and corresponding  confidence intervals
 
  #re-order coefficients if necessary
  est.cu <- c(p.mod6$b[2,], p.mod6$b[1,], p.mod6$b[4,], p.mod6$b[3,])
  ci.lb.cu <- c(p.mod6$ci.lb[2], p.mod6$ci.lb[1], p.mod6$ci.lb[4], p.mod6$ci.lb[3])
  ci.ub.cu <- c(p.mod6$ci.ub[2], p.mod6$ci.ub[1], p.mod6$ci.ub[4], p.mod6$ci.ub[3])

  #re-order coefficients if necessary
  est.pm <- c(p.mod3$b[2,], p.mod3$b[4,], p.mod3$b[1,], p.mod3$b[3,])
  ci.lb.pm <- c(p.mod3$ci.lb[2], p.mod3$ci.lb[4], p.mod3$ci.lb[1], p.mod3$ci.lb[3])
  ci.ub.pm <- c(p.mod3$ci.ub[2], p.mod3$ci.ub[4], p.mod3$ci.ub[1], p.mod3$ci.ub[3])

  p.estimates <- c(p.mod1[[1]], p.mod2[[1]], est.pm, p.mod4[[1]], p.mod5[[1]], est.cu, p.mod7[[1]])
  p.ci.lb  <- c(p.mod1[[5]], p.mod2[[5]], ci.lb.pm, p.mod4[[5]], p.mod5[[5]], ci.lb.cu, p.mod7[[5]])
  p.ci.ub <- c(p.mod1[[6]], p.mod2[[6]], ci.ub.pm, p.mod4[[6]], p.mod5[[6]], ci.ub.cu, p.mod7[[6]])
  

  ### create vector with labels
  p.labels <- c("","Chemical", "Physical",   "Leaves", "Vegetative", "Fruits", "Seeds", 
                "Non-harvested","Harvested", "Landrace", "Modern",   "Legumes", "Fruits", 
                "Vegetables", "Non-food", "Herb Annual", "Herb Perennial", "Woody Perennial")  
  p.rows <- c(27,  24.5:23.5,  21:18,  15.5:14.5,  12:11,  8.5:5.5,  3:1)

  # create vector of sample sizes (see supplementary materials)
  p.N <- c(31, 24, 10, 18, 2, 9, 4, 22, 15, 9, 18, 5, 11, 12,6, 15, 4, 14)
  names(p.N) <- c("Overall", "Chem", "Phys", "Lvs", "Veg", "Fruits", "Seeds", "Non-harv", "Harv", "Landrace",
      "Modern", "Legumes", "Fruits", "Veg", "Non-food", "Herb Ann", "Herb Per", "Woody Per")





### Create the forest plots  ###########################

dev.new(width=183, height=183, units="mm")
par(mfrow=c(1,2))


## For herbivore resistance
  par(font=1 )  
  forest(x= h.estimates, ci.lb=h.ci.lb, ci.ub=h.ci.ub, slab=h.labels, ilab=NULL, main="",
           annotate=FALSE, rows=h.rows, ylim=c(0.5,31), xlim=c(-10, 2), alim=c(-5,1), steps=7, 
           psize=12*(h.N/sum(h.N)), xlab="", digits=0,  at=c(-5, -3, -1, 1)) 

par(font=4) #bold italic
  text (-10, c(26.5, 22, 18.5, 15, 11.5, 4), c("Response Type", "Organ Measured", "Organ Use", "Dom Extent", "Crop Use", "Life History"), pos=4) 
par(font=2)  #bold
  text (-10, 28, "Overall", pos=4)
  text (-10, 31.5, "(a) Herbivore Resistance", pos=4)
par(font=1) #normal  
text (-2, 30, substitute(paste("Effect size (",italic("d"),") with 95% CI")))
axis(1,at=c(-5:1),labels=FALSE)

# Significance lines
t <- 0.1
lines(x=c(-2.5, -2.5), y=c(24.5, 25.5) )
lines(x=c(-2.5, -2.5+t),y=c(24.5, 24.5)) 
lines(x=c(-2.5, -2.5+t),y=c(25.5, 25.5)) 
text(-2.7, y=25, "*", srt=90)

lines(x=c(-3, -3), y=c(23.5, 25.5) )
lines(x=c(-3, -3+t),y=c(23.5, 23.5)) 
lines(x=c(-3, -3+t),y=c(25.5, 25.5)) 
text(-3.2, y=24.5, "***", srt=90)

lines(x=c(-3, -3), y=c(20, 21) )
lines(x=c(-3, -3+t),y=c(20, 20)) 
lines(x=c(-3, -3+t),y=c(21, 21)) 
text(-3.2, y=20.5, "**", srt=90)


 
# For plant traits
  forest(x= p.estimates, ci.lb=p.ci.lb, ci.ub=p.ci.ub, slab=p.labels, ilab=NULL, main="",
           annotate=FALSE, rows=p.rows, ylim=c(0.5,30), xlim=c(-25, 9), alim=c(-10,8), at=c(-10, -6, -2, 2, 5, 8), 
           psize=12*(p.N/sum(p.N)), xlab="", digits=0)    
 
par(font=4) #bold italic
  text (-25, c(25.5, 22, 16.5, 13, 9.5, 4), c("Response Type", "Organ Measured", "Organ Use", "Dom Extent", "Crop Use", "Life History"), pos=4) 
par(font=2)  #bold
  text (-25, 27, "Overall", pos=4)
  text (-25, 30.5, "(b) Plant Traits", pos=4)
par(font=1) #normal  
text (-1, 29, substitute(paste("Effect size (",italic("d"),") with 95% CI")))
axis(1,at=c(-10:8),labels=FALSE)

#adding significance lines
t <- 0.3
lines(x=c(-6, -6), y=c(19, 21) )
lines(x=c(-6, -6+t),y=c(19, 19)) 
lines(x=c(-6, -6+t),y=c(21, 21)) 
text(-6.6, y=20, "***", srt=90)

lines(x=c(-7.5, -7.5), y=c(18, 21) )
lines(x=c(-7.5, -7.5+t),y=c(18, 18)) 
lines(x=c(-7.5, -7.5+t),y=c(21, 21)) 
text(-8.1, y=19.5, "**", srt=90)

lines(x=c(-6, -6), y=c(14.5, 15.5) )
lines(x=c(-6, -6+t),y=c(14.5, 14.5)) 
lines(x=c(-6, -6+t),y=c(15.5, 15.5)) 
text(-6.6, y=15, "***", srt=90)



dev.copy2pdf(file="Meta_fig2_revised.pdf") 




save.image(file='Meta_finalrun_080416.RData')



