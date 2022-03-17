rm(list=ls())

####################################
## Meta-analysis of Money Priming ##
## R-script                       ##  
## Version: Sep 13th, 2018        ##
## ------------------------------ ##
## Paul Lodder MSc                ##
## How H. Ong BSc                 ##
## Raoul P.P.P. Grasman PhD       ##
## Jelte. M. Wicherts PhD         ##
####################################


################################################################################################
### Please set the R working directory to the location of the Meta-analysis Excel database ! ###
################################################################################################

# Install and load packages
install.packages("devtools")
install.packages("rJava")
install.packages("xlsx")
install_github("robbievanaert/puniform")
install.packages("puniform")
install.packages("weightr")
install.packages("metafor")
library(devtools)
library(rJava)
library(xlsx)
library(puniform)
library(weightr)
library(metafor)
library(ggplot2)

# Read Excel data from working directory
# The script reads the 'values' sheet in the Excel file. This is a duplicate of the 'Main' sheet (select all > copy > paste special > values). 
# The R function read.xlsx could not import the 'Main' sheet, but could succesfully import the copy in the 'values' sheet. 
metdat<-read.xlsx("MoneyPrimingMetaAnalysis_Dataset_FINAL_correction2.xlsx",sheetName="values",keepformulas=FALSE,startRow=2,header=T)

# Only select studies from dataset that meet inclusion criteria
metdat<-subset(metdat,metdat$Included.==1)

# Add standard error and p-value of Hedges' g
metdat<-cbind(metdat,gse=sqrt(metdat$var.of.g))
metdat<-cbind(metdat,gp=round(1-pnorm(metdat$g/metdat$gse),4))

# Prepare colors figure 1
metdat$Preregistered.col<-metdat$Preregistered.pch<-metdat$Preregistered
metdat$Preregistered.pch[which(metdat$Preregistered==1)]<-16
metdat$Preregistered.pch[which(metdat$Preregistered==0)]<-1
metdat$Preregistered.pch[which(is.na(metdat$Preregistered))]<-1
metdat$Preregistered.col[which(metdat$Preregistered==1)]<-"darkgreen"
metdat$Preregistered.col[which(metdat$Preregistered==0)]<-"black"
metdat$Preregistered.col[which(is.na(metdat$Preregistered))]<-"black"

# Rater agreement for coding study design factors
uniquecodes<-names(table(metdat$Study.code))
metdat.code<-metdat[1:length(uniquecodes),]
for(i in 1:length(uniquecodes)){
  ro<-which(metdat$Study.code==uniquecodes[i])[1]
  metdat.code[i,]<-metdat[ro,]
}
pr.agree<-length(which(metdat.code$Prime.Type..coder.agreement.==1))/length(uniquecodes)
se.agree<-length(which(metdat.code$Setting..code.agreement.==1))/length(uniquecodes)

# Function to compute significance of p-values
sig=function(metamodel){
  output=c()
  if(metamodel$pval[1]<=0.05){output<-"*"}
  if(metamodel$pval[1]<=0.01){output<-"**"}
  if(metamodel$pval[1]<=0.001){output<-"***"}
  if(metamodel$pval[1]>0.05){output<-""}
  print(output)
}


########################################
# -------- Main meta-analysis -------- #
########################################

### -------------------------- ###
### PREREGISTERED STUDIES ONLY ###

# Subset of studies based on interaction effects
metdat.pre=subset(metdat,metdat$Preregistered==1) # Subset of preregistered studies
metdat.pre.intno<-subset(metdat.pre,metdat.pre$Interaction.==0) # Create subset of studies without interaction effects
metdat.pre.intyes<-subset(metdat.pre,metdat.pre$Interaction.==1) # Create subset of studies with interaction effects
metdat.pre.inthigh<-subset(metdat.pre.intyes,metdat.pre.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.pre.high<-rbind(metdat.pre.inthigh,metdat.pre.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
fe.pre=rma(method="FE",yi=metdat.pre$g,vi=metdat.pre$var.of.g,data=metdat.pre) # Fixed effects analysis
fe.pre.mods=rma(method="FE",yi=metdat.pre$g,vi=metdat.pre$var.of.g,mods=~metdat.pre$gse,data=metdat.pre) # Fixed effects analysis with moderators
re.pre=rma(method="REML",yi=metdat.pre$g,vi=metdat.pre$var.of.g,data=metdat.pre) # Random effects meta-analysis
re.pre.high=rma(method="REML",yi=metdat.pre.high$g,vi=metdat.pre.high$var.of.g,data=metdat.pre.high) # Random effects meta-analysis including largest predicted interaction effects 
fe.pre.high.mods=rma(method="FE",yi=metdat.pre.high$g,vi=metdat.pre.high$var.of.g,mods=~metdat.pre.high$gse,data=metdat.pre.high) #  Fixed effects analysis with standard error moderator and including largest predicted interaction effects

### ---------------------- ###
### PUBLISHED STUDIES ONLY ###

# Subset of studies based on interaction effects
metdat.pub<-subset(metdat,metdat$published==1) # Subset of published studies
metdat.pub.intno<-subset(metdat.pub,metdat.pub$Interaction.==0) # Create subset of studies without interaction effects
metdat.pub.intyes<-subset(metdat.pub,metdat.pub$Interaction.==1) # Create subset of studies with interaction effects
metdat.pub.inthigh<-subset(metdat.pub.intyes,metdat.pub.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect      
metdat.pub.high<-rbind(metdat.pub.inthigh,metdat.pub.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
re.pub=rma(method="REML",yi=metdat.pub$g,vi=metdat.pub$var.of.g,data=metdat.pub) # Random effects meta-analysis
re.pub.high=rma(method="REML",yi=metdat.pub.high$g,vi=metdat.pub.high$var.of.g,data=metdat.pub.high) # Random effects meta-analysis largest interaction effects
fe.pub.high.mods=rma(method="FE",yi=metdat.pub.high$g,vi=metdat.pub.high$var.of.g,mods=~metdat.pub.high$gse,data=metdat.pub.high) # Random effects meta-analysis largest interaction effects
fe.pub.mods=rma(method="FE",yi=metdat.pub$g,vi=metdat.pub$var.of.g,mods=~metdat.pub$gse,data=metdat.pub.high) # Fixed effects meta-analysis with standard error moderator and including largest interaction effects


### ------------------------ ###
### UNPUBLISHED STUDIES ONLY ###

# Subset of studies based on interaction effects
metdat.unp=subset(metdat,metdat$published==0) # Subset of unpublished studies
metdat.unp.intno<-subset(metdat.unp,metdat.unp$Interaction.==0) # Create subset of studies without interaction effects
metdat.unp.intyes<-subset(metdat.unp,metdat.unp$Interaction.==1) # Create subset of studies with interaction effects
metdat.unp.inthigh<-subset(metdat.unp.intyes,metdat.unp.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect             
metdat.unp.intlow<-subset(metdat.unp.intyes,metdat.unp.intyes$Interaction.Identification..1.Largest.predicted.effect.==0) # Create subset of studies with smallest predicted interaction effect
metdat.unp.high<-rbind(metdat.unp.inthigh,metdat.unp.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects
metdat.unp.low<-rbind(metdat.unp.intlow,metdat.unp.intno) # Combine studies with smallest predicted interaction effect with studies without interaction effects

# Meta-analyses
re.unp=rma(method="REML",yi=metdat.unp$g,vi=metdat.unp$var.of.g,data=metdat.unp) # Random effects meta-analysis
re.unp.high=rma(method="REML",yi=metdat.unp.high$g,vi=metdat.unp.high$var.of.g,data=metdat.unp.high) # Random effects meta-analysis largest interaction effects
fe.unp.mods=rma(method="FE",yi=metdat.unp$g,vi=metdat.unp$var.of.g,mods=~metdat.unp$gse,data=metdat.unp) # Fixed effects meta-analysis with moderators 
fe.unp.high.mods=rma(method="FE",yi=metdat.unp.high$g,vi=metdat.unp.high$var.of.g,mods=~metdat.unp.high$gse,data=metdat.unp.high) # Fixed effects meta-analysis with moderators including largest interaction effects

### ----------------------- ###
### BEHAVIORAL STUDIES ONLY ###

# Subset of studies based on interaction effects
metdat.beh=subset(metdat,metdat$Behavioral.vs..Non.behavioral==1) # Subset of behavioral studies
metdat.beh.intno<-subset(metdat.beh,metdat.beh$Interaction.==0) # Create subset of studies without interaction effects
metdat.beh.intyes<-subset(metdat.beh,metdat.beh$Interaction.==1) # Create subset of studies with interaction effects
metdat.beh.inthigh<-subset(metdat.beh.intyes,metdat.beh.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.beh.high<-rbind(metdat.beh.inthigh,metdat.beh.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
fe.beh=rma(method="FE",yi=metdat.beh$g,vi=metdat.beh$var.of.g,data=metdat.beh) # Fixed effects analysis
fe.beh.mods=rma(method="FE",yi=metdat.beh$g,vi=metdat.beh$var.of.g,mods=~metdat.beh$gse,data=metdat.beh) # Fixed effects analysis with moderators
re.beh=rma(method="REML",yi=metdat.beh$g,vi=metdat.beh$var.of.g,data=metdat.beh) # Random effects meta-analysis
re.beh.high=rma(method="REML",yi=metdat.beh.high$g,vi=metdat.beh.high$var.of.g,data=metdat.beh.high) # Random effects meta-analysis including largest predicted interaction effects 
fe.beh.high.mods=rma(method="FE",yi=metdat.beh.high$g,vi=metdat.beh.high$var.of.g,mods=~metdat.beh.high$gse,data=metdat.beh.high) #  Fixed effects analysis with standard error moderator and including largest predicted interaction effects

### ----------------------- ###
###  BEHAVIORAL PUBLISHED   ###

# Subset of studies based on interaction effects
metdat.beh.pub=subset(metdat.beh,metdat.beh$published==1) # Subset of behavioral studies
metdat.beh.pub.intno<-subset(metdat.beh.pub,metdat.beh.pub$Interaction.==0) # Create subset of studies without interaction effects
metdat.beh.pub.intyes<-subset(metdat.beh.pub,metdat.beh.pub$Interaction.==1) # Create subset of studies with interaction effects
metdat.beh.pub.inthigh<-subset(metdat.beh.pub.intyes,metdat.beh.pub.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.beh.pub.high<-rbind(metdat.beh.pub.inthigh,metdat.beh.pub.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
fe.beh.pub=rma(method="FE",yi=metdat.beh.pub$g,vi=metdat.beh.pub$var.of.g,data=metdat.beh.pub) # Fixed effects analysis
fe.beh.pub.mods=rma(method="FE",yi=metdat.beh.pub$g,vi=metdat.beh.pub$var.of.g,mods=~metdat.beh.pub$gse,data=metdat.beh.pub) # Fixed effects analysis with moderators
re.beh.pub=rma(method="REML",yi=metdat.beh.pub$g,vi=metdat.beh.pub$var.of.g,data=metdat.beh.pub) # Random effects meta-analysis
re.beh.pub.high=rma(method="REML",yi=metdat.beh.pub.high$g,vi=metdat.beh.pub.high$var.of.g,data=metdat.beh.pub.high) # Random effects meta-analysis including largest predicted interaction effects 
fe.beh.pub.high.mods=rma(method="FE",yi=metdat.beh.pub.high$g,vi=metdat.beh.pub.high$var.of.g,mods=~metdat.beh.pub.high$gse,data=metdat.beh.pub.high) #  Fixed effects analysis with standard error moderator and including largest predicted interaction effects

### ------------------------- ###
###  BEHAVIORAL UNPUBLISHED   ###

# Subset of studies based on interaction effects
metdat.beh.unp=subset(metdat.beh,metdat.beh$published==0) # Subset of behavioral studies
metdat.beh.unp.intno<-subset(metdat.beh.unp,metdat.beh.unp$Interaction.==0) # Create subset of studies without interaction effects
metdat.beh.unp.intyes<-subset(metdat.beh.unp,metdat.beh.unp$Interaction.==1) # Create subset of studies with interaction effects
metdat.beh.unp.inthigh<-subset(metdat.beh.unp.intyes,metdat.beh.unp.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.beh.unp.high<-rbind(metdat.beh.unp.inthigh,metdat.beh.unp.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
fe.beh.unp=rma(method="FE",yi=metdat.beh.unp$g,vi=metdat.beh.unp$var.of.g,data=metdat.beh.unp) # Fixed effects analysis
fe.beh.unp.mods=rma(method="FE",yi=metdat.beh.unp$g,vi=metdat.beh.unp$var.of.g,mods=~metdat.beh.unp$gse,data=metdat.beh.unp) # Fixed effects analysis with moderators
re.beh.unp=rma(method="REML",yi=metdat.beh.unp$g,vi=metdat.beh.unp$var.of.g,data=metdat.beh.unp) # Random effects meta-analysis
re.beh.unp.high=rma(method="REML",yi=metdat.beh.unp.high$g,vi=metdat.beh.unp.high$var.of.g,data=metdat.beh.unp.high) # Random effects meta-analysis including largest predicted interaction effects 
fe.beh.unp.high.mods=rma(method="FE",yi=metdat.beh.unp.high$g,vi=metdat.beh.unp.high$var.of.g,mods=~metdat.beh.unp.high$gse,data=metdat.beh.unp.high) #  Fixed effects analysis with standard error moderator and including largest predicted interaction effects

### --------------------------- ###
### NON-BEHAVIORAL STUDIES ONLY ###

# Subset of studies based on interaction effects
metdat.nob=subset(metdat,metdat$Behavioral.vs..Non.behavioral==2) # Subset of non-behavioral studies
metdat.nob.intno<-subset(metdat.nob,metdat.nob$Interaction.==0) # Create subset of studies without interaction effects
metdat.nob.intyes<-subset(metdat.nob,metdat.nob$Interaction.==1) # Create subset of studies with interaction effects
metdat.nob.inthigh<-subset(metdat.nob.intyes,metdat.nob.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.nob.high<-rbind(metdat.nob.inthigh,metdat.nob.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
fe.nob=rma(method="FE",yi=metdat.nob$g,vi=metdat.nob$var.of.g,data=metdat.nob) # Fixed effects analysis
fe.nob.mods=rma(method="FE",yi=metdat.nob$g,vi=metdat.nob$var.of.g,mods=~metdat.nob$gse,data=metdat.nob) # Fixed effects analysis with moderators
re.nob=rma(method="REML",yi=metdat.nob$g,vi=metdat.nob$var.of.g,data=metdat.nob) # Random effects meta-analysis
re.nob.high=rma(method="REML",yi=metdat.nob.high$g,vi=metdat.nob.high$var.of.g,data=metdat.nob.high) # Random effects meta-analysis including largest predicted interaction effects 
fe.nob.high.mods=rma(method="FE",yi=metdat.nob.high$g,vi=metdat.nob.high$var.of.g,mods=~metdat.nob.high$gse,data=metdat.nob.high) #  Fixed effects analysis with standard error moderator and including largest predicted interaction effects


### --------------------------- ###
###  NON-BEHAVIORAL PUBLISHED   ###

# Subset of studies based on interaction effects
metdat.nob.pub=subset(metdat.nob,metdat.nob$published==1) # Subset of behavioral studies
metdat.nob.pub.intno<-subset(metdat.nob.pub,metdat.nob.pub$Interaction.==0) # Create subset of studies without interaction effects
metdat.nob.pub.intyes<-subset(metdat.nob.pub,metdat.nob.pub$Interaction.==1) # Create subset of studies with interaction effects
metdat.nob.pub.inthigh<-subset(metdat.nob.pub.intyes,metdat.nob.pub.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.nob.pub.high<-rbind(metdat.nob.pub.inthigh,metdat.nob.pub.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
fe.nob.pub=rma(method="FE",yi=metdat.nob.pub$g,vi=metdat.nob.pub$var.of.g,data=metdat.nob.pub) # Fixed effects analysis
fe.nob.pub.mods=rma(method="FE",yi=metdat.nob.pub$g,vi=metdat.nob.pub$var.of.g,mods=~metdat.nob.pub$gse,data=metdat.nob.pub) # Fixed effects analysis with moderators
re.nob.pub=rma(method="REML",yi=metdat.nob.pub$g,vi=metdat.nob.pub$var.of.g,data=metdat.nob.pub) # Random effects meta-analysis
re.nob.pub.high=rma(method="REML",yi=metdat.nob.pub.high$g,vi=metdat.nob.pub.high$var.of.g,data=metdat.nob.pub.high) # Random effects meta-analysis including largest predicted interaction effects 
fe.nob.pub.high.mods=rma(method="FE",yi=metdat.nob.pub.high$g,vi=metdat.nob.pub.high$var.of.g,mods=~metdat.nob.pub.high$gse,data=metdat.nob.pub.high) #  Fixed effects analysis with standard error moderator and including largest predicted interaction effects


### ----------------------------- ###
###  NON-BEHAVIORAL UNPUBLISHED   ###

# Subset of studies based on interaction effects
metdat.nob.unp=subset(metdat.nob,metdat.nob$published==0) # Subset of behavioral studies
metdat.nob.unp.intno<-subset(metdat.nob.unp,metdat.nob.unp$Interaction.==0) # Create subset of studies without interaction effects
metdat.nob.unp.intyes<-subset(metdat.nob.unp,metdat.nob.unp$Interaction.==1) # Create subset of studies with interaction effects
metdat.nob.unp.inthigh<-subset(metdat.nob.unp.intyes,metdat.nob.unp.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.nob.unp.high<-rbind(metdat.nob.unp.inthigh,metdat.nob.unp.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
fe.nob.unp=rma(method="FE",yi=metdat.nob.unp$g,vi=metdat.nob.unp$var.of.g,data=metdat.nob.unp) # Fixed effects analysis
fe.nob.unp.mods=rma(method="FE",yi=metdat.nob.unp$g,vi=metdat.nob.unp$var.of.g,mods=~metdat.nob.unp$gse,data=metdat.nob.unp) # Fixed effects analysis with moderators
re.nob.unp=rma(method="REML",yi=metdat.nob.unp$g,vi=metdat.nob.unp$var.of.g,data=metdat.nob.unp) # Random effects meta-analysis
re.nob.unp.high=rma(method="REML",yi=metdat.nob.unp.high$g,vi=metdat.nob.unp.high$var.of.g,data=metdat.nob.unp.high) # Random effects meta-analysis including largest predicted interaction effects 
fe.nob.unp.high.mods=rma(method="FE",yi=metdat.nob.unp.high$g,vi=metdat.nob.unp.high$var.of.g,mods=~metdat.nob.unp.high$gse,data=metdat.nob.unp.high) #  Fixed effects analysis with standard error moderator and including largest predicted interaction effects



### -------------------- ###
### ALL STUDIES COMBINED ###

# Subset of studies based on interaction effects
metdat.intno<-subset(metdat,metdat$Interaction.==0) # Create subset of studies without interaction effects
metdat.intyes<-subset(metdat,metdat$Interaction.==1) # Create subset of studies with interaction effects
metdat.inthigh<-subset(metdat.intyes,metdat.intyes$Interaction.Identification..1.Largest.predicted.effect.==1) # Create subset of studies with largest predicted interaction effect                
metdat.high<-rbind(metdat.inthigh,metdat.intno) # Combine studies with largest predicted interaction effect with studies without interaction effects

# Meta-analyses
re.all<-rma(method="REML",yi=metdat$g,vi=metdat$var.of.g,data=metdat) # Random effects meta-analysis
re.all.mv<-rma.mv(method="REML",yi=metdat$g,V=metdat$var.of.g,random=~1|metdat$Study.code,data=metdat) # Random effects meta-analysis with shared random effects for correlated rows in the data
re.high=rma(method="REML",yi=metdat.high$g,vi=metdat.high$var.of.g,data=metdat.high) # Random effects meta-analysis largest interaction effects
re.intno=rma(method="REML",yi=metdat.intno$g,vi=metdat.intno$var.of.g,data=metdat.intno) # Random effects meta-analysis smallest interaction effects
re.inthigh=rma(method="REML",yi=metdat.inthigh$g,vi=metdat.inthigh$var.of.g,data=metdat.inthigh) # Random effects meta-analysis smallest interaction effects
re.intno=rma(method="REML",yi=metdat.intno$g,vi=metdat.intno$var.of.g,data=metdat.intno) # Random effects meta-analysis smallest interaction effects
fe.high.mods=rma(method="FE",yi=metdat.high$g,vi=metdat.high$var.of.g,mods=~metdat.high$gse,data=metdat.high) # Random effects meta-analysis largest interaction effects
fe.inthigh.mods=rma(method="FE",yi=metdat.inthigh$g,vi=metdat.inthigh$var.of.g,mods=~metdat.inthigh$gse,data=metdat.inthigh) # Add moderators standard error & publication status
fe.intno.mods=rma(method="FE",yi=metdat.intno$g,vi=metdat.intno$var.of.g,mods=~metdat.intno$gse,data=metdat.intno) # Add moderators standard error & publication status
fe.all.mod<-rma(method="FE",yi=metdat$g,vi=metdat$var.of.g,mods=~metdat$gse,data=metdat) # Fixed effects meta-analysis with standard error moderator

# Table 1: all included studies
metdat.high<-cbind(metdat.high,Prime_coded=pr.types[metdat.high$Prime.Type..final.code.],Setting_coded=se.types[metdat.high$Setting..final.code.])
table1subset<-metdat.high[which(duplicated(metdat.high$Study.short.name)==FALSE),]
table1columns<-c("Study.short.name","jrnl","Setting_coded","Prime_coded","dep")
write.xlsx(table1subset[,table1columns],"table1.xlsx")


##########################################################
## Figure NATURE feature: plot effect size across years ##

fe.high.mods=rma(method="FE",yi=metdat.high$g,vi=metdat.high$var.of.g,mods=~as.numeric(metdat.high$year.),data=metdat.high) 
summary(fe.high.mods)
yeardat<-data.frame(year=metdat.high$year.,effect=metdat.high$g,
                    published=factor(metdat.high$published,levels=c(0,1),labels=c("Unpublished","Published")),
                    Prime=factor(metdat.high$Prime.Type..final.code.,levels=c(1,2,3,4,5),labels=c("Visual","Descrambling","Handling","Thinking","Multiple types")),
                    Outcome=factor(metdat.high$Behavioral.vs..Non.behavioral,levels=c(1,2),labels=c("Behavioral","Non-behavioral")))
yeardat<-yeardat[-which(yeardat$year=="Unknown"),]
yeardat$year<-droplevels(yeardat$year, exclude = if(anyNA(levels(yeardat$year))) NULL else NA)
nyear<-as.numeric(table(yeardat$year))
years<-data.frame(year=names(table(yeardat$year)),N=nyear)
yeardat$year<-factor(yeardat$year,levels=2005:2017,labels=c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017"))
#yeardat<-rbind(yeardat,c("2008",NA,"Published","Visual","Behavioral"))
#yeardat<-rbind(yeardat,c("2011",NA,"Published","Visual","Behavioral"))

# Scatterplot published
y<-ggplot(yeardat,aes(x=as.numeric(year),y=effect,color=published,linetype=published))
#y<-y+geom_abline(intercept = fe.high.mods$beta[1], slope = fe.high.mods$beta[2],col="lightblue",lwd=1.4)
y<-y+geom_abline(intercept = 0, slope = 0,col="black",linetype="dotted")
y<-y+theme_minimal(base_size = 15)
y<-y+ylab("Effect size (Hedges' g)")
y<-y+xlab("")
y<-y+geom_jitter(width=.2)
y<-y+geom_smooth(method='lm', formula= y~x,se=F,lwd=0.5)
y<-y+theme(axis.text.x=element_blank())
#y<-y+theme(axis.ticks.x=element_line(c("2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017")))
y<-y+annotate(geom = "text", x = 1:13, y = -.69, label = levels(yeardat$year), size = 4)
y<-y+coord_cartesian(xlim=c(0,14),ylim = c(-0.5, 3), expand = FALSE, clip = "off") 
y<-y+theme(legend.title = element_blank())
y
pdf(file="MoneyPriming_Publication_regline.pdf",width=9,height=3.5)
par(mar=c(8,5,2,2))
y
dev.off()

# Scatterplot behavioral
y<-ggplot(yeardat[-which(is.na(yeardat$Outcome)),],aes(x=as.numeric(year),y=effect,color=Outcome,linetype=Outcome))
#y<-y+geom_abline(intercept = fe.high.mods$beta[1], slope = fe.high.mods$beta[2],col="lightblue",lwd=1.4)
y<-y+geom_abline(intercept = 0, slope = 0,col="black",linetype="dotted")
y<-y+theme_minimal(base_size = 15)
y<-y+ylab("Effect size (Hedges' g)")
y<-y+xlab("")
y<-y+geom_jitter(width=.2)
y<-y+geom_smooth(method='lm', formula= y~x,se=F,lwd=0.5)
y<-y+theme(axis.text.x=element_blank())
y<-y+annotate(geom = "text", x = 1:11, y = -.69, label = years$year, size = 4)
y<-y+coord_cartesian(xlim=c(0,12),ylim = c(-0.5, 3), expand = FALSE, clip = "off") 
y
pdf(file="MoneyPriming_Behavioral_regline.pdf",width=8,height=3.5)
par(mar=c(8,5,2,2))
y
dev.off()

# Scatterplot prime type
y<-ggplot(yeardat,aes(x=as.numeric(year),y=effect,color=Prime,linetype=Prime))
#y<-y+geom_abline(intercept = fe.high.mods$beta[1], slope = fe.high.mods$beta[2],col="lightblue",lwd=1.4)
y<-y+geom_abline(intercept = 0, slope = 0,col="black",linetype="dotted")
y<-y+theme_minimal(base_size = 15)
y<-y+ylab("Effect size (Hedges' g)")
y<-y+xlab("")
y<-y+geom_jitter(width=.2)
y<-y+geom_smooth(method='lm', formula= y~x,se=F,lwd=0.5)
y<-y+theme(axis.text.x=element_blank())
y<-y+annotate(geom = "text", x = 1:11, y = -.69, label = years$year, size = 4)
y<-y+coord_cartesian(xlim=c(0,12),ylim = c(-0.5, 3), expand = FALSE, clip = "off") 
y
pdf(file="MoneyPriming_Primetype_regline.pdf",width=8,height=3.5)
par(mar=c(8,5,2,2))
y
dev.off()

# Boxplot (old)
y<-ggplot(yeardat,aes(x=year,y=effect))
y<-y+geom_abline(intercept = fe.high.mods$beta[1], slope = fe.high.mods$beta[2],col="lightblue",lwd=1.4)
y<-y+geom_abline(intercept = 0, slope = 0,col="red",linetype="dotted")
y<-y+theme_minimal(base_size = 15)
y<-y+ylab("Effect size (Hedges' g)")
y<-y+xlab("")
y<-y+geom_boxplot()
y<-y+annotate(geom = "text", x = years$year, y = -1, label = years$N, size = 3)
y<-y+annotate(geom = "text", x = -.27, y = -1, label = "# studies", size = 3,fontface="bold")
y<-y+annotate(geom = "text", x = -.50, y = -.69, label = "Year", size = 4,fontface="bold")
y<-y+coord_cartesian(xlim=c(0,12),ylim = c(-0.5, 3), expand = FALSE, clip = "off") 
pdf(file="DecliningMoneyPrimingEffect.pdf",width=7,height=3.5)
par(mar=c(8,5,2,2))
y
dev.off()


################################################

# Figure 1: Funnel plots for all studies, published studies, unpublished studies and preregistered studies
setEPS()
postscript("Figure1.eps",width=8.5)
#pdf("fun_multi.pdf",width=8.5)
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T))
par(mar=c(5,5,2,.2))
funnel(re.high,shade=c("White","Grey"),level=c(95,99),cex.lab=1.8,refline=0,xlab=expression(italic("g")),ylab=expression(italic("SE")),col=metdat.high$Preregistered.col,pch=metdat.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("All studies (k=",re.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.high$b,2)," [",round(re.high$ci.lb,2),", ",round(re.high$ci.ub,2),"]",sig(re.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.high.mods$b[2],2)," (p ",ifelse(fe.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.high$b,re.high$b),c(0,1),lty=2)
lines(c(re.high$b,re.high$b-1.96),c(0,1),lty=2)
lines(c(re.high$b,re.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.pub.high,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),cex.lab=1.8,refline=0,col=metdat.pub.high$Preregistered.col,pch=metdat.pub.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Published studies (k=",re.pub.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.pub.high$b,2)," [",round(re.pub.high$ci.lb,2),", ",round(re.pub.high$ci.ub,2),"]",sig(re.pub.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.pub.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.pub.high.mods$b[2],2)," (p ",ifelse(fe.pub.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.pub.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.pub.high$b,re.pub.high$b),c(0,1),lty=2)
lines(c(re.pub.high$b,re.pub.high$b-1.96),c(0,1),lty=2)
lines(c(re.pub.high$b,re.pub.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.unp.high,shade=c("White","Grey"),level=c(95,99),cex.lab=1.8,xlab=expression(italic("g")),ylab=expression(italic("SE")),refline=0,col=metdat.unp.high$Preregistered.col,pch=metdat.unp.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Unpublished studies (k=",re.unp.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.unp.high$b,2)," [",round(re.unp.high$ci.lb,2),", ",round(re.unp.high$ci.ub,2),"]",sig(re.unp.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.unp.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.unp.high.mods$b[2],2)," (p ",ifelse(fe.unp.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.unp.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.unp.high$b,re.unp.high$b),c(0,0.55),lty=2)
lines(c(re.unp.high$b,re.unp.high$b-1.96),c(0,1),lty=2)
lines(c(re.unp.high$b,re.unp.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.pre.high,cex.lab=1.8,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),pch=16,col="black",refline=0,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Pre-registered effects (k=",re.pre.high$k,")"),cex.main=1.5)  
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.pre.high$b,2)," [",round(re.pre.high$ci.lb,2),", ",round(re.pre.high$ci.ub,2),"]",sig(re.pre.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.pre.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.pre.high.mods$b[2],2)," (p ",ifelse(fe.pre.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.pre.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.pre.high$b,re.pre.high$b),c(0,0.55),lty=2)
lines(c(re.pre.high$b,re.pre.high$b-1.96),c(0,1),lty=2)
lines(c(re.pre.high$b,re.pre.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.intno,cex.lab=1.8,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("SE")),refline=0,pch=metdat.intno$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Main effects (k=",re.intno$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.intno$b,2)," [",round(re.intno$ci.lb,2),", ",round(re.intno$ci.ub,2),"]",sig(re.intno)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.intno$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egg = ",round(fe.intno.mods$b[2],2)," (p ",ifelse(fe.intno.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.intno.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.intno$b,re.intno$b),c(0,0.55),lty=2)
lines(c(re.intno$b,re.intno$b-1.96),c(0,1),lty=2)
lines(c(re.intno$b,re.intno$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.inthigh,cex.lab=1.8,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),refline=0,col=metdat.inthigh$Preregistered.col,pch=metdat.inthigh$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Interaction effects (k=",re.inthigh$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.inthigh$b,2)," [",round(re.inthigh$ci.lb,2),", ",round(re.inthigh$ci.ub,2),"]",sig(re.inthigh)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.inthigh$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.inthigh.mods$b[2],2)," (p ",ifelse(fe.inthigh.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.inthigh.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.inthigh$b,re.inthigh$b),c(0,0.55),lty=2)
lines(c(re.inthigh$b,re.inthigh$b-1.96),c(0,1),lty=2)
lines(c(re.inthigh$b,re.inthigh$b+1.96),c(0,1),lty=2)
dev.off()


# Figure A1 (APPENDIX): Funnel plots for all studies, published studies, unpublished studies and preregistered studies

setEPS()
postscript("FigureA1.eps",width=10)
#pdf("fun_appendix.pdf",width=10)
layout(matrix(c(1,2,3,4),2,2,byrow=T))
par(mar=c(5,5,2,.2))
funnel(re.all,shade=c("White","Grey"),level=c(95,99),cex.lab=1.8,refline=0,xlab=expression(italic("g")),ylab=expression(italic("SE")),col=metdat$Preregistered.col,pch=metdat$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("All studies (k=",re.all$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.15,0.15),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.all$b,2)," [",round(re.all$ci.lb,2),", ",round(re.all$ci.ub,2),"]",sig(re.all)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.all$I2,0),"%"),cex=1,pos=4)
#text(1.5,0.15,paste0("Egger = ",round(fe.high.mods$b[2],2)," (p ",ifelse(fe.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.all$b,re.all$b),c(0,1),lty=2)
lines(c(re.all$b,re.all$b-1.96),c(0,1),lty=2)
lines(c(re.all$b,re.all$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.pub,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),cex.lab=1.8,refline=0,col=metdat.pub$Preregistered.col,pch=metdat.pub$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Published studies (k=",re.pub$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.15,0.15),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.pub$b,2)," [",round(re.pub$ci.lb,2),", ",round(re.pub$ci.ub,2),"]",sig(re.pub)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.pub$I2,0),"%"),cex=1,pos=4)
#text(1.5,0.15,paste0("Egger = ",round(fe.pub.mods$b[2],2)," (p ",ifelse(fe.pub.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.pub.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.pub$b,re.pub$b),c(0,1),lty=2)
lines(c(re.pub$b,re.pub$b-1.96),c(0,1),lty=2)
lines(c(re.pub$b,re.pub$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.unp,shade=c("White","Grey"),level=c(95,99),cex.lab=1.8,xlab=expression(italic("g")),ylab=expression(italic("SE")),refline=0,col=metdat.unp$Preregistered.col,pch=metdat.unp$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Unpublished studies (k=",re.unp$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.15,0.15),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.unp$b,2)," [",round(re.unp$ci.lb,2),", ",round(re.unp$ci.ub,2),"]",sig(re.unp)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.unp$I2,0),"%"),cex=1,pos=4)
#text(1.5,0.15,paste0("Egger = ",round(fe.unp.mods$b[2],2)," (p ",ifelse(fe.unp.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.unp.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.unp$b,re.unp$b),c(0,0.55),lty=2)
lines(c(re.unp$b,re.unp$b-1.96),c(0,1),lty=2)
lines(c(re.unp$b,re.unp$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.pre,cex.lab=1.8,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),pch=16,col="black",refline=0,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Preregistered effects (k=",re.pre$k,")"),cex.main=1.5)  
rect(c(1.5,1.5),c(0.15,0.15),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.pre$b,2)," [",round(re.pre$ci.lb,2),", ",round(re.pre$ci.ub,2),"]",sig(re.pre)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.pre$I2,0),"%"),cex=1,pos=4)
#text(1.5,0.15,paste0("Egger = ",round(fe.pre.mods$b[2],2)," (p ",ifelse(fe.pre.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.pre.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.pre$b,re.pre$b),c(0,0.55),lty=2)
lines(c(re.pre$b,re.pre$b-1.96),c(0,1),lty=2)
lines(c(re.pre$b,re.pre$b+1.96),c(0,1),lty=2)
dev.off()




# Figure 2: Funnel plots for all studies, published studies, unpublished studies and preregistered studies
setEPS()
postscript("Figure2.eps",width=8.5)
#pdf("fun_multi.pdf",width=8.5)
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T))
par(mar=c(5,5,2,.2))
funnel(re.beh.high,shade=c("White","Grey"),level=c(95,99),cex.lab=1.8,refline=0,xlab=expression(italic("g")),ylab=expression(italic("SE")),pch=metdat.beh.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Behavioral (k=",re.beh.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.beh.high$b,2)," [",round(re.beh.high$ci.lb,2),", ",round(re.beh.high$ci.ub,2),"]",sig(re.beh.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.beh.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.beh.high.mods$b[2],2)," (p ",ifelse(fe.beh.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.beh.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.beh.high$b,re.beh.high$b),c(0,1),lty=2)
lines(c(re.beh.high$b,re.beh.high$b-1.96),c(0,1),lty=2)
lines(c(re.beh.high$b,re.beh.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.nob.high,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),cex.lab=1.8,refline=0,pch=metdat.nob.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Non-behavioral (k=",re.nob.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.nob.high$b,2)," [",round(re.nob.high$ci.lb,2),", ",round(re.nob.high$ci.ub,2),"]",sig(re.nob.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.nob.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.nob.high.mods$b[2],2)," (p ",ifelse(fe.nob.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.nob.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.nob.high$b,re.nob.high$b),c(0,1),lty=2)
lines(c(re.nob.high$b,re.nob.high$b-1.96),c(0,1),lty=2)
lines(c(re.nob.high$b,re.nob.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.beh.pub.high,shade=c("White","Grey"),level=c(95,99),cex.lab=1.8,xlab=expression(italic("g")),ylab=expression(italic("SE")),refline=0,pch=metdat.beh.pub.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Published behavioral (k=",re.beh.pub.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.beh.pub.high$b,2)," [",round(re.beh.pub.high$ci.lb,2),", ",round(re.beh.pub.high$ci.ub,2),"]",sig(re.beh.pub.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.beh.pub.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.beh.pub.high.mods$b[2],2)," (p ",ifelse(fe.beh.pub.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.beh.pub.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.beh.pub.high$b,re.beh.pub.high$b),c(0,0.55),lty=2)
lines(c(re.beh.pub.high$b,re.beh.pub.high$b-1.96),c(0,1),lty=2)
lines(c(re.beh.pub.high$b,re.beh.pub.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.nob.pub.high,cex.lab=1.8,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),refline=0,pch=metdat.nob.pub.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Published non-behavioral (k=",re.nob.pub.high$k,")"),cex.main=1.5)  
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.nob.pub.high$b,2)," [",round(re.nob.pub.high$ci.lb,2),", ",round(re.nob.pub.high$ci.ub,2),"]",sig(re.nob.pub.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.nob.pub.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.nob.pub.high.mods$b[2],2)," (p ",ifelse(fe.nob.pub.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.nob.pub.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.nob.pub.high$b,re.nob.pub.high$b),c(0,0.55),lty=2)
lines(c(re.nob.pub.high$b,re.nob.pub.high$b-1.96),c(0,1),lty=2)
lines(c(re.nob.pub.high$b,re.nob.pub.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.beh.unp.high,cex.lab=1.8,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("SE")),refline=0,pch=metdat.beh.unp.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Unpublished behavioral (k=",re.beh.unp.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.beh.unp.high$b,2)," [",round(re.beh.unp.high$ci.lb,2),", ",round(re.beh.unp.high$ci.ub,2),"]",sig(re.beh.unp.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.beh.unp.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egg = ",round(fe.beh.unp.high.mods$b[2],2)," (p ",ifelse(fe.beh.unp.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.beh.unp.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.beh.unp.high$b,re.beh.unp.high$b),c(0,0.55),lty=2)
lines(c(re.beh.unp.high$b,re.beh.unp.high$b-1.96),c(0,1),lty=2)
lines(c(re.beh.unp.high$b,re.beh.unp.high$b+1.96),c(0,1),lty=2)
par(mar=c(5,5,2,.2))
funnel(re.nob.unp.high,cex.lab=1.8,shade=c("White","Grey"),level=c(95,99),xlab=expression(italic("g")),ylab=expression(italic("")),refline=0,pch=metdat.nob.unp.high$Preregistered.pch,ylim=c(0,0.5),xlim=c(-1,3),main=paste0("Unpublished non-behavioral (k=",re.nob.unp.high$k,")"),cex.main=1.5)
rect(c(1.5,1.5),c(0.2,0.2),c(3.1,3.1),c(0,0),col="white")
text(1.5,0.05,paste0("g = ",round(re.nob.unp.high$b,2)," [",round(re.nob.unp.high$ci.lb,2),", ",round(re.nob.unp.high$ci.ub,2),"]",sig(re.nob.unp.high)),cex=1,pos=4)
text(1.5,0.10,paste0("I^2 = ",round(re.nob.unp.high$I2,0),"%"),cex=1,pos=4)
text(1.5,0.15,paste0("Egger = ",round(fe.nob.unp.high.mods$b[2],2)," (p ",ifelse(fe.nob.unp.high.mods$pval[2]<.001,"<.001",paste0("= ",round(fe.nob.unp.high.mods$pval[2],3))),")"),cex=1,pos=4)
lines(c(re.nob.unp.high$b,re.nob.unp.high$b),c(0,0.55),lty=2)
lines(c(re.nob.unp.high$b,re.nob.unp.high$b-1.96),c(0,1),lty=2)
lines(c(re.nob.unp.high$b,re.nob.unp.high$b+1.96),c(0,1),lty=2)
dev.off()






######################################################
#  ---- Meta-regressions based on study design  ---- #
######################################################

# Only use this option when including besides the predicted interaction effects also all other levels of the interaction (Supplemental analyses)
# metdat.high<-metdat

# Recode aggregated DV's
dep2<-unlist(as.character(metdat.high$dep))
dep2[grep("DV1:",dep2)]<-"Aggregated"
metdat.high$dep<-dep2

# Select DV's that are used in 5 or more experiments
manydvs<-table(metdat.high$dep)[which(table(metdat.high$dep)>4)]
dv.types<-names(manydvs)

# Vector with names of money priming types
pr.types<-c("Visual",
            "Descrambling",
            "Handling",
            "Thinking",
            "Combination")

# Vector with names of study settings
se.types<-c("Lab",
            "Online",
            "Field",
            "Unknown")

# Vector with dependent measure types
beh.types<-c("Behavioral","Non-behavioral")

# Subset of studies with frequently used DV
metdat.high.dv<-metdat.high[which(metdat.high$dep %in% dv.types),]

# Create dummy categories dependent measure types 
dum1=dum2=dum3=dum4=dum5=dum6=dum7=dum8=dum9<-as.character(unlist(metdat.high.dv$dep))
dum<-cbind(dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9)
for(i in 1:dim(dum)[2]){
  for(j in 1:dim(dum)[1]){
    ifelse(as.character(dum[j,i]==dv.types[i]),dum[j,i]<-1,dum[j,i]<-0)
  }}
dum<-dum[,-1] # Remove dummy variable of reference level

# Create dummy categories money priming types (combination prime is reference)
dum.visu=dum.desc=dum.hand=dum.thnk=dum.comb<-metdat.high$Prime.Type..coder.1.
dum1<-data.frame(dum.visu,dum.desc,dum.hand,dum.thnk,dum.comb)
for(i in 1:5){
  for(j in 1:length(dum1[,1])){
    if(dum1[j,i]!=i) {dum1[j,i]<-0}
    if(dum1[j,i]==i) {dum1[j,i]<-1}
  }}
dum1<-dum1[,-5] # Remove dummy variable of reference level

# Create dummy categories study settings (online setting is reference)
dum.lab=dum.online=dum.field<-metdat.high$Setting..coder.1.
dum2<-data.frame(dum.lab,dum.online,dum.field)
for(i in 1:3){
  for(j in 1:length(dum2[,1])){
    if(dum2[j,i]!=i) {dum2[j,i]<-0}
    if(dum2[j,i]==i) {dum2[j,i]<-1}
  }}
dum2<-dum2[,-2] # Remove dummy variable of reference level

# Add dummy code columns to dataset
moddata<-cbind(metdat.high,dum1,dum2)
moddatadv<-cbind(metdat.high.dv,dum)

# Table 2: Meta-regression analyses with prime type and study setting
re.prime=rma(method="REML",yi=moddata$g,vi=moddata$var.of.g,mods=~moddata$dum.visu+moddata$dum.desc+moddata$dum.hand+moddata$dum.thnk,data=moddata)
re.set=rma(method="REML",yi=moddata$g,vi=moddata$var.of.g,mods=~moddata$dum.lab+moddata$dum.field,data=moddata)
moddata$Behavioral.vs..Non.behavioral[moddata$Behavioral.vs..Non.behavioral==2]<-0 # recode non-behavioral studies to be the reference category
moddata$Behavioral.vs..Non.behavioral[moddata$Behavioral.vs..Non.behavioral==3]<-NA # temporarily exclude studies with both DV types
re.dv=rma(method="REML",yi=moddata$g,vi=moddata$var.of.g,mods=~moddata$Behavioral.vs..Non.behavioral,data=moddata)

# Table 2: Meta-regressions other study characteristics
re.all.mod=rma(method="REML",yi=metdat$g,vi=metdat$var.of.g,
               mods=~metdat$gse+metdat$published+metdat$Preregistered+metdat$Multiple.DV.+metdat$Interaction.)
re.high.mod=rma(method="REML",yi=metdat.high$g,vi=metdat.high$var.of.g,
                mods=~metdat.high$gse+metdat.high$published+metdat.high$Preregistered+metdat.high$Multiple.DV.+metdat.high$Interaction.)



##################################################
#  ---- SEPARATE META-ANALYSES FOR DV TYPES ---- #
##################################################

# Create empty table 4
dvres<-pun<-egg<-psm<-list()
ndv<-length(dv.types)
dvtable<-data.frame(Outcome=1:ndv,k=1:ndv,g=1:ndv,gp=1:ndv,gsig=1:ndv,Q=1:ndv,Qp=1:ndv,psi=1:ndv,pp=1:ndv,tau2=1:ndv,I2=1:ndv,pu=1:ndv,puci=1:ndv,pup=1:ndv,pusi=1:ndv,egg=1:ndv,eggci=1:ndv,eggp=1:ndv,esig=1:ndv,sm=1:ndv,smp=1:ndv,smsi=1:ndv,smci=1:ndv,smpb=1:ndv,smpbsi=1:ndv)

# Fill table 4 in loop
for(i in 1:ndv){
  print(i)
  metdat.high.temp<-metdat.high.dv[which(metdat.high.dv$dep==dv.types[i]),]
  
  # Random effects meta-analysis
  dvres[[i]]<-rma(method="REML",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,data=metdat.high.temp)
  
  # P uniform
  t <- try(puniform(yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,side="right",method="P"))
  ifelse("try-error" %in% class(t),pun[[i]]<-NA,pun[[i]]<-puniform(yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,side="right",method="P"))
  
  # Fixed effects with SE moderator
  t <- try(rma(method="FE",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,mods=~metdat.high.temp$gse,data=metdat.high.temp))
  ifelse("try-error" %in% class(t),egg[[i]]<-NA,egg[[i]]<-rma(method="FE",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,mods=~metdat.high.temp$gse,data=metdat.high.temp))
  
  # Selection model
  t <- try(psm[[i]]<-weightfunct(effect=metdat.high.temp$g,v=metdat.high.temp$var.of.g))
  ifelse("try-error" %in% class(t),psm[[i]]<-NA,psm[[i]]<-weightfunct(effect=metdat.high.temp$g,v=metdat.high.temp$var.of.g))
  rm(tmp)
  tmp <- suppressWarnings(capture.output(psm[[i]]))
  tmp.d <- suppressWarnings(as.numeric(strsplit(tmp[20], " ")[[1]]))
  tmp.pb <- suppressWarnings(as.numeric(strsplit(tmp[24], " ")[[1]]))
  tmp.d <- tmp.d[which(is.na(tmp.d) == FALSE)]
  tmp.pb <- tmp.pb[which(is.na(tmp.pb) == FALSE)]
  dvtable$sm[i]<-tmp.d[1]
  dvtable$smp[i]<-tmp.d[4]
  dvtable$smci[i]<-paste0(round(tmp.d[1],2)," [",round(tmp.d[5],2),", ",round(tmp.d[6],2),"]")
  dvtable$smpb[i]<-tmp.pb
  dvtable[i,]$gsig[which(dvres[[i]]$pval<=.05)]<-"*"
  dvtable[i,]$gsig[which(dvres[[i]]$pval<.01)]<-"**"
  dvtable[i,]$gsig[which(dvres[[i]]$pval<.001)]<-"***"
  dvtable[i,]$gsig[which(dvres[[i]]$pval>.05)]<-""
  dvtable[i,]$Outcome<-dv.types[i]
  dvtable[i,]$gp<-round(dvres[[i]]$pval,3)
  dvtable[i,]$k<-dvres[[i]]$k
  dvtable[i,]$g<-paste0(round(dvres[[i]]$b,2)," [",round(dvres[[i]]$ci.lb,2),", ",round(dvres[[i]]$ci.ub,2),"]")
  dvtable[i,]$Q<-round(dvres[[i]]$QE,2)
  dvtable[i,]$Qp=round(dvres[[i]]$QEp,3)
  dvtable[i,]$tau2<-round(dvres[[i]]$tau2,2)
  dvtable[i,]$I2<-round(dvres[[i]]$I2,1)
  dvtable[i,]$egg<-round(egg[[i]]$b[2],2)
  dvtable[i,]$eggci<-paste0(round(egg[[i]]$b[2],2)," [",round(egg[[i]]$ci.lb[2],2),", ",round(egg[[i]]$ci.ub[2],2),"]")
  dvtable[i,]$eggp<-round(egg[[i]]$pval[2],3)
  dvtable[i,]$pu<-round(pun[[i]]$est,2)
  dvtable[i,]$puci<-paste0(round(pun[[i]]$est,2)," [",round(pun[[i]]$ci.lb,2),", ",round(pun[[i]]$ci.ub,2),"]")
  dvtable[i,]$pup<-round(pun[[i]]$pval.pb,3)
  dvtable[i,]$pp<-round(pun[[i]]$pval.0,3)
  dvtable[i,]$esig[which(dvtable$eggp[i]<=.05)]<-"*"
  dvtable[i,]$esig[which(dvtable$eggp[i]<.01)]<-"**"
  dvtable[i,]$esig[which(dvtable$eggp[i]<.001)]<-"***"
  dvtable[i,]$esig[which(dvtable$eggp[i]>.05)]<-""
  dvtable[i,]$pusi[which(dvtable$pup[i]<=.05)]<-"*"
  dvtable[i,]$pusi[which(dvtable$pup[i]<.01)]<-"**"
  dvtable[i,]$pusi[which(dvtable$pup[i]<.001)]<-"***"
  dvtable[i,]$pusi[which(dvtable$pup[i]>.05)]<-""
  dvtable[i,]$psi[which(dvtable$pp[i]<=.05)]<-"*"
  dvtable[i,]$psi[which(dvtable$pp[i]<.01)]<-"**"
  dvtable[i,]$psi[which(dvtable$pp[i]<.001)]<-"***"
  dvtable[i,]$psi[which(dvtable$pp[i]>.05)]<-""
  dvtable[i,]$smpbsi[which(dvtable$smpb[i]<=.05)]<-"*"
  dvtable[i,]$smpbsi[which(dvtable$smpb[i]<.01)]<-"**"
  dvtable[i,]$smpbsi[which(dvtable$smpb[i]<.001)]<-"***"
  dvtable[i,]$smpbsi[which(dvtable$smpb[i]>.05)]<-""
  dvtable[i,]$smsi[which(dvtable$smp[i]<=.05)]<-"*"
  dvtable[i,]$smsi[which(dvtable$smp[i]<.01)]<-"**"
  dvtable[i,]$smsi[which(dvtable$smp[i]<.001)]<-"***"
  dvtable[i,]$smsi[which(dvtable$smp[i]>.05)]<-""
}

# Summarise results
dvresult<-data.frame(Outcome=dvtable$Outcome,
                     k=dvtable$k,
                     g=paste0(dvtable$g,dvtable$gsig),
                     I2=paste0(round(dvtable$I2,1),"%"),
                     Egger=paste0(dvtable$eggci,dvtable$esig),
                     Puniform=paste0(dvtable$puci,dvtable$psi),
                     ThreePSM=paste(dvtable$smci,dvtable$smsi))
View(dvresult)


# Figure 4: Funnel plots of DV subsets with 5 or more studies
setEPS()
postscript("Figure4.eps",width=8,height=6)
#pdf(file="dvfun.pdf",width=8,height=6)
par(mar=c(4,4,4,1))
layout(matrix(seq(1,9,1),3,3,byrow=T))
for(i in 1:ndv){
  funnel(dvres[[i]],refline=0,xlim=c(-1,3),xlab=expression(italic("g")),ylab=expression(italic("SE")),cex.lab=1.3,cex=1.3,ylim=c(.5,0),shade=c("White","Grey"),level=c(95,99))
  title(main=dv.types[i])
  rect(c(1.3,1.3),c(0.34,0.34),c(3.1,3.1),c(0,0),col="white")
  text(1.3,0.05,paste0("g = ",round(dvres[[i]]$b,2),dvtable$gsig[i]),cex=1,pos=4)
  text(1.3,0.11,paste0("I^2 = ",round(dvres[[i]]$I2,0),"%"),cex=1,pos=4)
  text(1.3,0.17,ifelse(dvtable$esig[i]=="***","Egger < .001***", paste0("Egger = ",round(dvtable$eggp[i],3),dvtable$esig[i])),cex=1,pos=4)
  text(1.3,0.23,ifelse(dvtable$smpbsi[i]=="***","3PSM < .001***", paste0("3PSM = ",round(dvtable$smpb[i],3),dvtable$smpbsi[i])),cex=1,pos=4)
  text(1.3,0.29,ifelse(dvtable$pusi[i]=="***","Punif < .001***", paste0("Punif = ",round(dvtable$pup[i],3),dvtable$pusi[i])),cex=1,pos=4)
  lines(c(dvres[[i]]$b,dvres[[i]]$b),c(0,1),lty=2)
  lines(c(dvres[[i]]$b,dvres[[i]]$b-1.96),c(0,1),lty=2)
  lines(c(dvres[[i]]$b,dvres[[i]]$b+1.96),c(0,1),lty=2)
  par(font=2)
  rect(c(-1.15,-1.15),c(0.110,0.110),c(-.75,-.75),c(0,0),col="white")
  text(-.95,0.05,i,cex=1.3)
  par(font=1)
}
dev.off()



#####################################################
#  ---- SEPARATE META-ANALYSES FOR PRIME TYPES ---- #
#####################################################

# Subset of studies with frequently used DV
pr.types2<-1:5
metdat.high.pr<-metdat.high[which(metdat.high$Prime.Type..final.code. %in% pr.types2),]
prres<-pun<-egg<-psm<-list()
npr<-length(pr.types2)
prtable<-data.frame(Outcome=1:npr,k=1:npr,g=1:npr,gp=1:npr,gsig=1:npr,Q=1:npr,Qp=1:npr,tau2=1:npr,I2=1:npr,pu=1:npr,puci=1:npr,pup=1:npr,pusi=1:npr,egg=1:npr,eggci=1:npr,eggp=1:npr,esig=1:npr)
for(i in 1:npr){
  metdat.high.temp<-metdat.high.pr[which(metdat.high.pr$Prime.Type..final.code.==pr.types2[i]),]
  
  # Random effects meta-analysis
  prres[[i]]<-rma(method="REML",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,data=metdat.high.temp)
  
  # P uniform
  t <- try(puniform(yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,side="right",method="P"))
  ifelse("try-error" %in% class(t),pun[[i]]<-NA,pun[[i]]<-puniform(yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,side="right",method="P"))
  
  # Fixed effects with SE moderator
  t <- try(rma(method="FE",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,mods=~metdat.high.temp$gse,data=metdat.high.temp))
  ifelse("try-error" %in% class(t),egg[[i]]<-NA,egg[[i]]<-rma(method="FE",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,mods=~metdat.high.temp$gse,data=metdat.high.temp))
  
  # Selection model
  psm[[i]]<-weightfunct(effect=metdat.high.temp$g,v=metdat.high.temp$var.of.g)
  #psm.out[[i]]<-capture.output(psm[[i]])
  
  prtable[i,]$gsig[which(prres[[i]]$pval<=.05)]<-"*"
  prtable[i,]$gsig[which(prres[[i]]$pval<.01)]<-"**"
  prtable[i,]$gsig[which(prres[[i]]$pval<.001)]<-"***"
  prtable[i,]$gsig[which(prres[[i]]$pval>.05)]<-""
  prtable[i,]$Outcome<-pr.types2[i]
  prtable[i,]$gp<-round(prres[[i]]$pval,3)
  prtable[i,]$k<-prres[[i]]$k
  prtable[i,]$g<-paste0(round(prres[[i]]$b,2)," [",round(prres[[i]]$ci.lb,2),", ",round(prres[[i]]$ci.ub,2),"]")
  prtable[i,]$Q<-round(prres[[i]]$QE,2)
  prtable[i,]$Qp=round(prres[[i]]$QEp,3)
  prtable[i,]$tau2<-round(prres[[i]]$tau2,2)
  prtable[i,]$I2<-round(prres[[i]]$I2,0)
  prtable[i,]$egg<-round(egg[[i]]$b[2],2)
  prtable[i,]$eggci<-paste0(round(egg[[i]]$b[2],2)," [",round(egg[[i]]$ci.lb[2],2),", ",round(egg[[i]]$ci.ub[2],2),"]")
  prtable[i,]$eggp<-round(egg[[i]]$pval[2],3)
  prtable[i,]$pu<-round(pun[[i]]$est,2)
  prtable[i,]$puci<-paste0(round(pun[[i]]$est,2)," [",round(pun[[i]]$ci.lb,2),", ",round(pun[[i]]$ci.ub,2),"]")
  prtable[i,]$pup<-round(pun[[i]]$pval.pb,3)
  prtable[i,]$esig[which(prtable$eggp[i]<=.05)]<-"*"
  prtable[i,]$esig[which(prtable$eggp[i]<.01)]<-"**"
  prtable[i,]$esig[which(prtable$eggp[i]<.001)]<-"***"
  prtable[i,]$esig[which(prtable$eggp[i]>.05)]<-""
  prtable[i,]$pusi[which(prtable$pup[i]<=.05)]<-"*"
  prtable[i,]$pusi[which(prtable$pup[i]<.01)]<-"**"
  prtable[i,]$pusi[which(prtable$pup[i]<.001)]<-"***"
  prtable[i,]$pusi[which(prtable$pup[i]>.05)]<-""
}


# Extra Figure (not in paper): Funnel plots of prime types

pdf(file="prfun.pdf",width=14,height=2.5)
par(mar=c(3,2,4,2))
layout(matrix(seq(1,5,1),1,5,byrow=T))
for(i in 1:npr){
  funnel(prres[[i]],refline=0,xlim=c(-1,3),ylim=c(.5,0))
  title(main=paste(pr.types[i]))
  rect(c(1,1),c(0.2,0.2),c(3,3),c(0,0),col="white")
  text(1.2,0.04,paste0("g = ",round(prres[[i]]$b,2),prtable$gsig[i]),cex=1,pos=4)
  text(1.2,0.08,paste0("I^2 = ",round(prres[[i]]$I2,0),"%"),cex=1,pos=4)
  text(1.2,0.12,paste("Egg=",prtable$egg[i],prtable$esig[i]),cex=1,pos=4)
  text(1.2,0.16,paste("Pu=",prtable$pu[i],prtable$pusi[i]),cex=1,pos=4)
  lines(c(prres[[i]]$b,prres[[i]]$b),c(0,1),lty=2)
  lines(c(prres[[i]]$b,prres[[i]]$b-1.96),c(0,1),lty=2)
  lines(c(prres[[i]]$b,prres[[i]]$b+1.96),c(0,1),lty=2)
}
dev.off()





#######################################################
#  ---- SEPARATE META-ANALYSES FOR SETTING TYPES ---- #
#######################################################

# Subset of studies with frequently used DV
set.types2<-1:3
metdat.high.set<-metdat.high[which(metdat.high$Setting..final.code. %in% set.types2),]
setres<-pun<-egg<-psm<-list()
nset<-length(set.types2)
settable<-data.frame(Outcome=1:nset,k=1:nset,g=1:nset,gp=1:nset,gsig=1:nset,Q=1:nset,Qp=1:nset,tau2=1:nset,I2=1:nset,pu=1:nset,puci=1:nset,pup=1:nset,pusi=1:nset,egg=1:nset,eggci=1:nset,eggp=1:nset,esig=1:nset)
for(i in 1:nset){
  metdat.high.temp<-metdat.high.set[which(metdat.high.set$Setting..final.code.==set.types2[i]),]
  
  # Random effects meta-analysis
  setres[[i]]<-rma(method="REML",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,data=metdat.high.temp)
  
  # P uniform
  t <- try(puniform(yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,side="right",method="P"))
  ifelse("try-error" %in% class(t),pun[[i]]<-NA,pun[[i]]<-puniform(yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,side="right",method="P"))
  
  # Fixed effects with SE moderator
  t <- try(rma(method="FE",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,mods=~metdat.high.temp$gse,data=metdat.high.temp))
  ifelse("try-error" %in% class(t),egg[[i]]<-NA,egg[[i]]<-rma(method="FE",yi=metdat.high.temp$g,vi=metdat.high.temp$var.of.g,mods=~metdat.high.temp$gse,data=metdat.high.temp))
  
  # Selection model
  psm[[i]]<-weightfunct(effect=metdat.high.temp$g,v=metdat.high.temp$var.of.g)
  
  settable[i,]$gsig[which(setres[[i]]$pval<=.05)]<-"*"
  settable[i,]$gsig[which(setres[[i]]$pval<.01)]<-"**"
  settable[i,]$gsig[which(setres[[i]]$pval<.001)]<-"***"
  settable[i,]$gsig[which(setres[[i]]$pval>.05)]<-""
  settable[i,]$Outcome<-set.types2[i]
  settable[i,]$gp<-round(setres[[i]]$pval,3)
  settable[i,]$k<-setres[[i]]$k
  settable[i,]$g<-paste0(round(setres[[i]]$b,2)," [",round(setres[[i]]$ci.lb,2),", ",round(setres[[i]]$ci.ub,2),"]")
  settable[i,]$Q<-round(setres[[i]]$QE,2)
  settable[i,]$Qp=round(setres[[i]]$QEp,3)
  settable[i,]$tau2<-round(setres[[i]]$tau2,2)
  settable[i,]$I2<-round(setres[[i]]$I2,0)
  settable[i,]$egg<-round(egg[[i]]$b[2],2)
  settable[i,]$eggci<-paste0(round(egg[[i]]$b[2],2)," [",round(egg[[i]]$ci.lb[2],2),", ",round(egg[[i]]$ci.ub[2],2),"]")
  settable[i,]$eggp<-round(egg[[i]]$pval[2],3)
  settable[i,]$pu<-round(pun[[i]]$est,2)
  settable[i,]$puci<-paste0(round(pun[[i]]$est,2)," [",round(pun[[i]]$ci.lb,2),", ",round(pun[[i]]$ci.ub,2),"]")
  settable[i,]$pup<-round(pun[[i]]$pval.pb,3)
  settable[i,]$esig[which(settable$eggp[i]<=.05)]<-"*"
  settable[i,]$esig[which(settable$eggp[i]<.01)]<-"**"
  settable[i,]$esig[which(settable$eggp[i]<.001)]<-"***"
  settable[i,]$esig[which(settable$eggp[i]>.05)]<-""
  settable[i,]$pusi[which(settable$pup[i]<=.05)]<-"*"
  settable[i,]$pusi[which(settable$pup[i]<.01)]<-"**"
  settable[i,]$pusi[which(settable$pup[i]<.001)]<-"***"
  settable[i,]$pusi[which(settable$pup[i]>.05)]<-""
}


# Extra figure (not in paper): Funnel plots of setime types

pdf(file="setfun.pdf",width=8.4,height=2.5)
par(mar=c(3,2,4,2))
layout(matrix(seq(1,3,1),1,3,byrow=T))
for(i in 1:nset){
  funnel(setres[[i]],refline=0,xlim=c(-1,3),ylim=c(.5,0))
  title(main=paste(se.types[i]))
  rect(c(1,1),c(0.2,0.2),c(3,3),c(0,0),col="white")
  text(1.2,0.04,paste0("g = ",round(setres[[i]]$b,2),settable$gsig[i]),cex=1,pos=4)
  text(1.2,0.08,paste0("I^2 = ",round(setres[[i]]$I2,0),"%"),cex=1,pos=4)
  text(1.2,0.12,paste("Egg=",settable$egg[i],settable$esig[i]),cex=1,pos=4)
  text(1.2,0.16,paste("Pu=",settable$pu[i],settable$pusi[i]),cex=1,pos=4)
  lines(c(setres[[i]]$b,setres[[i]]$b),c(0,1),lty=2)
  lines(c(setres[[i]]$b,setres[[i]]$b-1.96),c(0,1),lty=2)
  lines(c(setres[[i]]$b,setres[[i]]$b+1.96),c(0,1),lty=2)
}
dev.off()




############################################################################
#  ---- CREATE HOMOGENEOUS SUBSETS BASED ON SETTING, PRIME & DV TYPES ---- #
############################################################################

# Before creating subsets, exclude studies with unkonwn design types
metdat.high<-metdat.high[-which(metdat.high$Setting..final.code.==99),]
metdat.high<-metdat.high[-which(metdat.high$Behavioral.vs..Non.behavioral==3),]

# Split dataset based on setting, prime and dv (behavioral or non-behavioral)
moddatas<-split(metdat.high,list(metdat.high$Prime.Type..final.code.,metdat.high$Setting..final.code.,metdat.high$Behavioral.vs..Non.behavioral))

# Create empty lists to store results of subset meta-analyses
moddatas.fe=moddatas.fm=moddatas.re=moddatas.pu<-list() 

# Create empty dataframe to store important estimates of subset meta-analyses
moddatas.base=data.frame(k=rep(0,length(moddatas)),
                         Setting=rep(0,length(moddatas)),Prime=rep(0,length(moddatas)),DV=rep(0,length(moddatas)),
                         I2=rep(0,length(moddatas)),Qm=rep(0,length(moddatas)),Qm.p=rep(0,length(moddatas)),
                         Pu=rep(0,length(moddatas)),Pu.p=rep(0,length(moddatas)),Pu.ci=rep(0,length(moddatas)),Pb.p=rep(0,length(moddatas)),
                         Egg=rep(0,length(moddatas)),Eggp=rep(0,length(moddatas)),Eggci=rep(0,length(moddatas)),
                         sm=rep(0,length(moddatas)),smp=rep(0,length(moddatas)),smci=rep(0,length(moddatas)),smpb=rep(0,length(moddatas)),
                         g=rep(0,length(moddatas)),gp=rep(0,length(moddatas)),gci=rep(0,length(moddatas)),mean.p=rep(0,length(moddatas)))

# Run meta-analyses on each split dataset (fixed effects, fixed effects with se mod, random effects, puniform)
for(i in 1:length(moddatas)){
  print(i)
  moddatas.base$k[i]<-dim(moddatas[[i]])[1]
  moddatas.base$Setting[i]<-se.types[moddatas[[i]]$Setting..final.code.[1]]
  moddatas.base$Prime[i]<-pr.types[moddatas[[i]]$Prime.Type..final.code.[1]]  
  moddatas.base$DV[i]<-moddatas[[i]]$Behavioral.vs..Non.behavioral[1]
  
  
  # Create NA if no studies in subset  
  if(dim(moddatas[[i]])[1]==0){
    moddatas.fe[[i]]<-moddatas.fm[[i]]<-moddatas.re[[i]]<-moddatas.pu[[i]]<-NA
  }
  
  # Create NA if error, if not then insert statistics in dataframe
  if(dim(moddatas[[i]])[1]>0){
    
    # Fixed effects model
    t <- try(rma(method="FE",yi=moddatas[[i]]$g,vi=moddatas[[i]]$var.of.g,data=moddatas[[i]]))
    ifelse("try-error" %in% class(t),moddatas.fe[[i]]<-NA,moddatas.fe[[i]]<-rma(method="FE",yi=moddatas[[i]]$g,vi=moddatas[[i]]$var.of.g,data=moddatas[[i]]))
    
    # Random effects model
    t <- try(rma(method="REML",yi=moddatas[[i]]$g,vi=moddatas[[i]]$var.of.g,data=moddatas[[i]]))
    ifelse("try-error" %in% class(t),moddatas.re[[i]]<-NA,moddatas.re[[i]]<-t)
    
    # P uniform
    t <- try(puniform(yi=moddatas[[i]]$g,vi=moddatas[[i]]$var.of.g,side="right",method="P"),silent=T)
    ifelse("try-error" %in% class(t),moddatas.pu[[i]]<-NA,moddatas.pu[[i]]<-puniform(yi=moddatas[[i]]$g,vi=moddatas[[i]]$var.of.g,side="right",method="P"))
    
    # Fixed effects with SE moderator
    t <- try(rma(method="FE",yi=moddatas[[i]]$g,vi=moddatas[[i]]$var.of.g,mods=~moddatas[[i]]$gse,data=moddatas[[i]]))
    ifelse("try-error" %in% class(t),moddatas.fm[[i]]<-NA,moddatas.fm[[i]]<-rma(method="FE",yi=moddatas[[i]]$g,vi=moddatas[[i]]$var.of.g,mods=~moddatas[[i]]$gse,data=moddatas[[i]]))
    
    ifelse(length(which(moddatas[[i]]$gp<.05))>0,moddatas.base$mean.p[i]<-mean(moddatas[[i]]$gp[which(moddatas[[i]]$gp<.05)]),moddatas.base$mean.p[i]<-0)
    
    # Selection model
    t <- try(weightfunct(effect=moddatas[[i]]$g,v=moddatas[[i]]$var.of.g),silent=T)
    ifelse("try-error" %in% class(t),psm[[i]]<-NA,psm[[i]]<-weightfunct(effect=moddatas[[i]]$g,v=moddatas[[i]]$var.of.g))
    
  }
}

# Select only subsets with 5 or more studies
nr<-which(moddatas.base$k>4)

# Combine results in presentable dataframe
moddatas.def<-moddatas.base[nr,]
for(i in 1:dim(moddatas.def)[1]){
  print(i)
  moddatas.def$g[i]<-moddatas.re[[nr[i]]]$b[1]
  moddatas.def$gp[i]<-moddatas.re[[nr[i]]]$pval[1]
  moddatas.def$gci[i]<-paste0(round(moddatas.re[[nr[i]]]$b[1],2)," [",round(moddatas.re[[nr[i]]]$ci.lb[1],2),", ",round(moddatas.re[[nr[i]]]$ci.ub[1],2),"]")
  moddatas.def$I2[i]<-moddatas.re[[nr[i]]]$I2
  moddatas.def$Qm[i]<-moddatas.fm[[nr[i]]]$QM  
  moddatas.def$Qm.p[i]<-moddatas.fm[[nr[i]]]$QMp  
  moddatas.def$Egg[i]<-moddatas.fm[[nr[i]]]$b[2] 
  moddatas.def$Eggp[i]<-moddatas.fm[[nr[i]]]$pval[2]
  moddatas.def$Eggci[i]<-paste0(round(moddatas.fm[[nr[i]]]$b[2],2)," [",round(moddatas.fm[[nr[i]]]$ci.lb[2],2),", ",round(moddatas.fm[[nr[i]]]$ci.ub[2],2),"]")
  
  # P-uniform
  t <- try(moddatas.pu[[nr[i]]]$est)
  ifelse("try-error" %in% class(t),moddatas.def$Pu[i]<-NA,moddatas.def$Pu[i]<-moddatas.pu[[nr[i]]]$est)
  ifelse("try-error" %in% class(t),moddatas.def$Pu.ci[i]<-NA,moddatas.def$Pu.ci[i]<-paste0(round(moddatas.pu[[nr[i]]]$est,2)," [",round(moddatas.pu[[nr[i]]]$ci.lb,2),", ",round(moddatas.pu[[nr[i]]]$ci.ub,2),"]"))
  t <- try(moddatas.pu[[nr[i]]]$pval.0)
  ifelse("try-error" %in% class(t),moddatas.def$Pu.p[i]<-NA,moddatas.def$Pu.p[i]<-moddatas.pu[[nr[i]]]$pval.0)
  t <- try(moddatas.pu[[nr[i]]]$pval.pb)
  ifelse("try-error" %in% class(t),moddatas.def$Pb.p[i]<-NA,moddatas.def$Pb.p[i]<-moddatas.pu[[nr[i]]]$pval.pb)
  
  # Selection model
  rm(tmp)
  tmp <- suppressWarnings(capture.output(psm[[nr[i]]]))
  tmp.d <- suppressWarnings(as.numeric(strsplit(tmp[20], " ")[[1]]))
  tmp.pb <- suppressWarnings(as.numeric(strsplit(tmp[24], " ")[[1]]))
  tmp.d <- tmp.d[which(is.na(tmp.d) == FALSE)]
  tmp.pb <- tmp.pb[which(is.na(tmp.pb) == FALSE)]
  moddatas.def$sm[i]<-tmp.d[1]
  moddatas.def$smp[i]<-tmp.d[4]
  moddatas.def$smci[i]<-paste0(round(tmp.d[1],2)," [",round(tmp.d[5],2),", ",round(tmp.d[6],2),"]")
  moddatas.def$smpb[i]<-tmp.pb
  
}



# Prepare funnel plot dataset
mdata<-moddatas.def
Qsig=gsig=psig=smsig=smpbsig=pbsig=esig<-as.vector(rep("",dim(mdata)[1]))
mdata<-data.frame(mdata,Qsig,gsig,esig,psig,pbsig,smsig,smpbsig,stringsAsFactors = FALSE)

# Create significance notation in funnel plot
mdata$Qsig[which(mdata$Qm.p<=.05)]<-"*"
mdata$Qsig[which(mdata$Qm.p<.01)]<-"**"
mdata$Qsig[which(mdata$Qm.p<.001)]<-"***"
mdata$gsig[which(mdata$gp<=.05)]<-"*"
mdata$gsig[which(mdata$gp<.01)]<-"**"
mdata$gsig[which(mdata$gp<.001)]<-"***"
mdata$psig[which(mdata$Pu.p<=.05)]<-"*"
mdata$psig[which(mdata$Pu.p<.01)]<-"**"
mdata$psig[which(mdata$Pu.p<.001)]<-"***"
mdata$smsig[which(mdata$smp<=.05)]<-"*"
mdata$smsig[which(mdata$smp<.01)]<-"**"
mdata$smsig[which(mdata$smp<.001)]<-"***"
mdata$smpbsig[which(mdata$smpb<=.05)]<-"*"
mdata$smpbsig[which(mdata$smpb<.01)]<-"**"
mdata$smpbsig[which(mdata$smpb<.001)]<-"***"
mdata$pbsig[which(mdata$Pb.p<=.05)]<-"*"
mdata$pbsig[which(mdata$Pb.p<.01)]<-"**"
mdata$pbsig[which(mdata$Pb.p<.001)]<-"***"
mdata$esig[which(mdata$Eggp<=.05)]<-"*"
mdata$esig[which(mdata$Eggp<.01)]<-"**"
mdata$esig[which(mdata$Eggp<.001)]<-"***"

# Table 3: Meta-analytic statistics of money priming subsets with 4 or more studies
print(mdata)
mdata$DV<-ifelse(mdata$DV==1,mdata$DV<-"Behavioral","Non-behavioral")
subres<-data.frame(Prime=mdata$Prime,
                   Setting=mdata$Setting,
                   DV=mdata$DV,
                   k=mdata$k,
                   g=paste0(mdata$gci,mdata$gsig),
                   I2=paste0(round(mdata$I2,1),"%"),
                   Egger=paste0(mdata$Eggci,mdata$esig),
                   Puniform=paste0(mdata$Pu.ci,mdata$psig),
                   ThreePSM=paste(mdata$smci,mdata$smsig))
View(subres)


# Figure 3: Funnel plots of subsets with 5 or more studies
flist<-moddatas.re[nr]
setEPS()
postscript("Figure3.eps",width=9.5,height=12)
#pdf(file="modfun.pdf",width=9.5,height=9)
par(mar=c(4,4,4,1))
layout(matrix(seq(1,15,1),5,3,byrow=T))
for(i in 1:dim(mdata)[1]){
  funnel(flist[[i]],refline=0,xlim=c(-1,3),xlab=expression(italic("g")),ylab=expression(italic("SE")),cex.lab=1.3,cex=1.3,ylim=c(.5,0),shade=c("White","Grey"),level=c(95,99))
  title(main=paste( mdata$DV[i],mdata$Setting[i],"Studies",
                    "\nPrime = ",mdata$Prime[i]), cex.main=1.3)
  lines(c(mdata$g[i],mdata$g[i]),c(0,1),lty=2)
  lines(c(mdata$g[i],mdata$g[i]-1.96),c(0,1),lty=2)
  lines(c(mdata$g[i],mdata$g[i]+1.96),c(0,1),lty=2)
  rect(c(1.4,1.4),c(0.28,0.28),c(3.05,3.05),c(0,0),col="white")
  text(1.4,0.04,paste0("g = ",round(mdata$g[i],2),mdata$gsig[i]),cex=1,pos=4)
  text(1.4,0.09,paste0(expression(I^2)," = ",round(mdata$I2[i],0),"%"),cex=1,pos=4)
  text(1.4,0.14,ifelse(mdata$Qsig[i]=="***","Egger < .001***", paste0("Egger = ",round(mdata$Eggp[i],3),mdata$esig[i])),cex=1,pos=4)
  text(1.4,0.24,ifelse(mdata$pbsig[i]=="***","Punif < .001***", paste0("Punif = ",round(mdata$Pb.p[i],3),mdata$pbsig[i])),cex=1,pos=4)
  text(1.4,0.19,ifelse(mdata$smpbsig[i]=="***","3PSM < .001***", paste0("3PSM = ",round(mdata$smpb[i],3),mdata$smpbsig[i])),cex=1,pos=4)
  par(font=2)
  rect(c(-1.15,-1.15),c(0.110,0.110),c(-.75,-.75),c(0,0),col="white")
  text(-.95,0.05,i,cex=1.3)
  par(font=1)
}
dev.off()




