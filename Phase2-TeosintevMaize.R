---
title: "Teosinte and Maize Comparison"
author: "Alonso Favela"
date: "5/20/2019"
----
  
##Running Needed packages
library(phyloseq)
library(ggplot2)
library(vegan)

##Set working directory 
setwd('~/Documents/Projects/Research/Maize Microbiome/MM2015/Data Analysis /CSV folder')


###16S rRNA Statisitcs ======

##Load Data 
MM.16S <- read.csv("16S_STATs_table_rare34000.csv")

##Subset the genotypes of interest 
TvMcomparison=subset(MM.16S,Genotype == "A632" |Genotype == "B73"
| Genotype == "LH82"
| Genotype == "PH207"
| Genotype == "Mo17"
| Genotype == "Pa91"
| Genotype == "Ames21786"
| Genotype == "Ames21789"
| Genotype == "Ames21809"
| Genotype == "PI566674"
| Genotype == "PI566677"
| Genotype == "PI566680"
#| Genotype == "Bulk" 
)
##This last step modulates the inclusion of the bulk soils sample. Not sure I need to include it in the genetic and domesticated analysis 


##Subsetting the values of interest 
MM.16S.treatment <- TvMcomparison[,1:12]
MM.16S.spec <- TvMcomparison[,13:15131]

#Genotype comparison on th microibome || I need to remove Bulk
adonis(MM.16S.spec ~ Genotype, MM.16S.treatment)
#Types effect on the whole microbiome unaveraged. R2 0.1067
#This model includes Genotype,subspecies (Type), and domestication (Family)
##Unaverged includes all sample models in analysis 
adonis(MM.16S.spec ~ Genotype+Type+Family, MM.16S.treatment)
#Subspecies model
adonis(MM.16S.spec ~ Type, MM.16S.treatment)
#Domestication
adonis(MM.16S.spec ~ Family, MM.16S.treatment)

#Here will be the mean values version of this figure 
#Make sure to change all of the names this can come to be annoying
TvMcomparison.GenotypicMeans<-aggregate(TvMcomparison[, (14:14926)], list(Genotype=TvMcomparison$Genotype, 
                                                                          Family=TvMcomparison$Family,
                                                                          Type=TvMcomparison$Type), mean)
#This is to pull out the data 
Genotype.16S.treatment <- TvMcomparison.GenotypicMeans[,1:3]
Genotype.16S.spec <- TvMcomparison.GenotypicMeans[,4:14916]

##This type factor explained 66% of the microbial community, 
adonis(Genotype.16S.spec ~ Type, Genotype.16S.treatment)
adonis(Genotype.16S.spec ~ Family, Genotype.16S.treatment)

##Supplmental Table X: ANOVA microbiome =====
##These models below includes the bluk soil which give an over estimate. 
# Type : sub species 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Type       3   0.37355 0.124515  6.0025 0.66676  0.001 ***
#   Residuals  9   0.18669 0.020744         0.33324           
# Total     12   0.56024                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##The residuals don't completely make sense here. 
# ## I should use the family test===== 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Family     2   0.33884 0.16942  7.5365 0.60117  0.001 ***
#   Residuals 10   0.22480 0.02248         0.39883           
# Total     12   0.56364                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###I wanted to add some NMDS Analysis 

MM.16S.treatment <- MM.16S[,1:12]
MM.16S.spec <- MM.16S[,13:15131]

ord<-metaMDS(MM.16S.spec)
plot(ord, type = "t")

NMDSAxis<-ord$points
NMDSAxis<-as.data.frame(NMDSAxis)
res.aov<-aov(NMDSAxis$MDS1 ~ MM.16S.treatment$Genotype)

summary(res.aov)

##NMDS1 0.7898/(0.7898+1.9912)==0.2839986
# Df Sum Sq  Mean Sq F value   Pr(>F)    
# MM.16S.treatment$Genotype  36 0.7898 0.021939    3.57 5.48e-10 ***
#   Residuals                 324 1.9912 0.006146                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##All MC NMDS1 1.123/(1.123+1.323)= 0.4591169
# Df Sum Sq  Mean Sq F value Pr(>F)    
# MM.16S.treatment$Genotype  36  1.123 0.031206   7.643 <2e-16 ***
#   Residuals                 324  1.323 0.004083                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##This is just the 16S communtiy 
##Broad Sense H2 SSmean/SStotal = H2  NMDS NMDS 0.5190/(0.6013+0.5190) =0.4632688
# Df Sum Sq Mean Sq F value   Pr(>F)    
# MM.16S.treatment$Genotype  11 0.5190 0.04719   8.475 1.31e-10 ***
#   Residuals                 108 0.6013 0.00557                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##NMDS2 SSmean/SStotal = H2  NMDS NMDS 0.2376/(0.2376+0.7063) =0.2517216
# Df Sum Sq Mean Sq F value   Pr(>F)    
# MM.16S.treatment$Genotype  11 0.2376 0.02160   3.304 0.000614 ***
#   Residuals                 108 0.7063 0.00654                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##This suggest that the communtiy is under some degree of selecion., This is our communtiy level hereitbiltiy 

##Fungal ITS Statisitcs =====
MM.AITS <- read.csv("otu_table_ITS_MaizeGDB.csv")
MM.AITS.treatment <- MM.AITS[,1:18]
MM.AITS.spec <- MM.AITS[,19:1045]
##Caluclating communtiy heretiability 
ord<-metaMDS(MM.AITS.spec)
#plot(ord, type = "t")

NMDSAxis<-ord$points
NMDSAxis<-as.data.frame(NMDSAxis)
res.aov<-aov(NMDSAxis$MDS1 ~ MM.AITS.treatment$Genotype)
summary(res.aov)
## 2.438/(2.438+8.303)=0.2269807
# Df Sum Sq Mean Sq F value   Pr(>F)    
# MM.AITS.treatment$Genotype  35  2.438 0.06967    2.71 2.34e-06 ***
#   Residuals                  323  8.303 0.02571                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
res.aov<-aov(NMDSAxis$MDS2 ~ MM.AITS.treatment$Genotype)
summary(res.aov)
# 3.047/(3.047+7.016)=0.3027924
# MM.AITS.treatment$Genotype  35  3.047 0.08704   4.007 1.37e-11 ***
#   Residuals                  323  7.016 0.02172                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##Subset the genotypes of interest 
TvM.ITS=subset(MM.AITS,Genotype == "A632" |Genotype == "B73"
                     | Genotype == "LH82"
                     | Genotype == "PH207"
                     | Genotype == "Mo17"
                     | Genotype == "Pa91"
                     | Genotype == "Ames21786"
                     | Genotype == "Ames21789"
                     | Genotype == "Ames21809"
                     | Genotype == "PI566674"
                     | Genotype == "PI566677"
                     | Genotype == "PI566680"
                     #| Genotype == "Bulk" )


##Subsetting the values of interest 
MM.ITS.treatment <- TvM.ITS[,1:18]
MM.ITS.spec <- TvM.ITS[,19:1045]

#Genotype comparison on th microibome
adonis(MM.ITS.spec ~ Genotype, MM.ITS.treatment)
adonis(MM.ITS.spec ~ Type, MM.ITS.treatment)
adonis(MM.ITS.spec~Family,MM.ITS.treatment)
#Types effect on the whole ITS unaveraged. R2 0.5

#Here will be the mean values version of this figure 
#Make sure to change all of the names this can come to be annoying
TvM.ITS.Means<-aggregate(TvM.ITS[, (19:1045)], list(Genotype=TvM.ITS$Genotype, 
                                                                          Family=TvM.ITS$Family,
                                                                          Type=TvM.ITS$Type), mean)

MM.ITS.treatment <- TvM.ITS.Means[,1:3]
MM.ITS.spec <- TvM.ITS.Means[,4:1030]

adonis(MM.ITS.spec ~ Type, MM.ITS.treatment)
# erms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# Type       2   0.44138 0.22069  2.0635 0.31439  0.002 **
#   Residuals  9   0.96254 0.10695         0.68561          
# Total     11   1.40393                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.ITS.spec ~ Family, MM.ITS.treatment)
##Use this for supplmental tables 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# Family     1    0.2981 0.29810  2.6957 0.21233  0.004 **
#   Residuals 10    1.1058 0.11058         0.78767          
# Total     11    1.4039                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Begining of Functional Genes Amplicon Data======
###Bacterial amoA ====
MM.BamoA <- read.csv("otu_table_AmoA-Bact.csv")


#Subet genotypes
TvM.BamoA=subset(MM.BamoA,Genotype == "A632" |Genotype == "B73"
               | Genotype == "LH82"
               | Genotype == "PH207"
               | Genotype == "Mo17"
               | Genotype == "Pa91"
               | Genotype == "Ames21786"
               | Genotype == "Ames21789"
               | Genotype == "Ames21809"
               | Genotype == "PI566674"
               | Genotype == "PI566677"
               | Genotype == "PI566680"
               | Genotype == "Bulk" )
##Subsetting the values of interest 
MM.BamoA.treatment <- TvM.BamoA[,1:13]
MM.BamoA.spec <- TvM.BamoA[,14:204]

#Genotype comparison on th microibome
adonis(MM.BamoA.spec ~ Genotype, MM.BamoA.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Genotype   11    2.4427 0.22206  1.1769 0.10884  0.098 .
# Residuals 106   20.0011 0.18869         0.89116         
# Total     117   22.4438                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis(MM.BamoA.spec ~ Type, MM.BamoA.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Type        2    0.6558 0.32788  1.7306 0.02922  0.019 *
#   Residuals 115   21.7880 0.18946         0.97078         
# Total     117   22.4438                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TvM.BamoA.Means<-aggregate(TvM.BamoA[, (14:204)], list(Genotype=TvM.BamoA$Genotype, 
                                                   Family=TvM.BamoA$Family,
                                                   Type=TvM.BamoA$Type), mean)

MM.bamoa.treatment <- TvM.BamoA.Means[,1:3]
MM.bamoa.spec <- TvM.BamoA.Means[,4:194]

##BamoA Functional Gene statisitc==== 
adonis(MM.bamoa.spec ~ Type, MM.bamoa.treatment)
# Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)  
# Type       2   0.13369 0.066844  1.7596 0.2811  0.032 *
#   Residuals  9   0.34189 0.037988         0.7189         
# Total     11   0.47558                  1.0000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.bamoa.spec ~ Family, MM.bamoa.treatment)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# Family     1   0.10311 0.103111  2.7683 0.21681  0.009 **
#   Residuals 10   0.37247 0.037247         0.78319          
# Total     11   0.47558                  1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


MM.BamoA.treatment <- MM.BamoA[,1:13]
MM.BamoA.spec <- MM.BamoA[,14:204]

ord<-metaMDS(MM.BamoA.spec)
#plot(ord, type = "t")

NMDSAxis<-ord$points
NMDSAxis<-as.data.frame(NMDSAxis)
res.aov<-aov(NMDSAxis$MDS1 ~ MM.BamoA.treatment$Genotype)
summary(res.aov)

res.aov<-aov(NMDSAxis$MDS2 ~ MM.BamoA.treatment$Genotype)
summary(res.aov)



###Arch amoA ====
MM.AamoA <- read.csv("otu_table_amoA_Arch.csv")
MM.AamoA.treatment <- MM.AamoA[,1:13]
MM.AamoA.spec <- MM.AamoA[,14:105]

#Subet genotypes
TvM.AamoA=subset(MM.AamoA,Genotype == "A632" |Genotype == "B73"
                 | Genotype == "LH82"
                 | Genotype == "PH207"
                 | Genotype == "Mo17"
                 | Genotype == "Pa91"
                 | Genotype == "Ames21786"
                 | Genotype == "Ames21789"
                 | Genotype == "Ames21809"
                 | Genotype == "PI566674"
                 | Genotype == "PI566677"
                 | Genotype == "PI566680"
                 | Genotype == "Bulk" )
##Subsetting the values of interest 
MM.AamoA.treatment <- TvM.AamoA[,1:13]
MM.AamoA.spec <- TvM.AamoA[,14:105]

#Genotype comparison on th microibome
adonis(MM.AamoA.spec ~ Genotype, MM.AamoA.treatment)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Genotype   11    1.4324 0.130214  1.7788 0.15338  0.001 ***
#   Residuals 108    7.9061 0.073205         0.84662           
# Total     119    9.3385                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis(MM.AamoA.spec ~ Type, MM.AamoA.treatment)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Type        2    0.4433 0.221665  2.9156 0.04747  0.001 ***
#   Residuals 117    8.8951 0.076027         0.95253           
# Total     119    9.3385                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 

TvM.AamoA.Means<-aggregate(TvM.AamoA[, (14:105)], list(Genotype=TvM.AamoA$Genotype, 
                                                       Family=TvM.AamoA$Family,
                                                       Type=TvM.AamoA$Type), mean)

MM.AamoA.treatment <- TvM.AamoA.Means[,1:3]
MM.AamoA.spec <- TvM.AamoA.Means[,4:94]

##AamoA Functional Gene statisitc==== 
adonis(MM.AamoA.spec ~ Type, MM.AamoA.treatment)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# Type       2  0.047654 0.023827  1.5871 0.26073  0.087 .
# Residuals  9  0.135122 0.015014         0.73927         
# Total     11  0.182776                  1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.AamoA.spec ~ Family, MM.AamoA.treatment)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Family     1  0.020754 0.020754  1.2809 0.11355  0.224
# Residuals 10  0.162023 0.016202         0.88645       
# Total     11  0.182776                  1.00000       



ord<-metaMDS(MM.AamoA.spec)
#plot(ord, type = "t")

NMDSAxis<-ord$points
NMDSAxis<-as.data.frame(NMDSAxis)
res.aov<-aov(NMDSAxis$MDS1 ~ MM.AamoA.treatment$Genotype)
summary(res.aov)
#3.408/(3.408+13.520)=0.20
res.aov<-aov(NMDSAxis$MDS2 ~ MM.AamoA.treatment$Genotype)
summary(res.aov)
#2.323/(2.323+12.671)=0.154928
##NirS ====
MM.nirS <- read.csv("otu_table_nirS.csv")

MM.nir.treatment <- MM.nir[,1:15]
MM.nir.spec <- MM.nir[,16:358]

#Subet genotypes
TvM.nirS=subset(MM.nirS, Genotype == "A632" |Genotype == "B73"
                 | Genotype == "LH82"
                 | Genotype == "PH207"
                 | Genotype == "Mo17"
                 | Genotype == "Pa91"
                 | Genotype == "Ames21786"
                 | Genotype == "Ames21789"
                 | Genotype == "Ames21809"
                 | Genotype == "PI566674"
                 | Genotype == "PI566677"
                 | Genotype == "PI566680"
                 | Genotype == "Bulk" )
##Subsetting the values of interest 
MM.nirS.treatment <- TvM.nirS[,1:15]
MM.nirS.spec <- TvM.nirS[,16:358]

#Genotype comparison on th microibome
adonis(MM.nirS.spec ~ Genotype, MM.nirS.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   11     5.595 0.50864  1.3864 0.12373  0.001 ***
#   Residuals 108    39.624 0.36689         0.87627           
# Total     119    45.219                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis(MM.nirS.spec ~ Type, MM.nirS.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type        2     1.761 0.88064  2.3709 0.03895  0.001 ***
#   Residuals 117    43.458 0.37143         0.96105           
# Total     119    45.219                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TvM.nirS.Means<-aggregate(TvM.nirS[, (16:358)], list(Genotype=TvM.nirS$Genotype, 
                                                       Family=TvM.nirS$Family,
                                                       Type=TvM.nirS$Type), mean)

MM.nirS.treatment <- TvM.nirS.Means[,1:3]
MM.nirS.spec <- TvM.nirS.Means[,4:346]

##nirS Functional Gene statisitc==== 
adonis(MM.nirS.spec ~ Type, MM.nirS.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type       2   0.54183 0.27092  1.5355 0.25441  0.001 ***
#   Residuals  9   1.58794 0.17644         0.74559           
# Total     11   2.12978                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.nirS.spec ~ Family, MM.nirS.treatment)
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)   
# Family     1   0.33949 0.33949  1.8963 0.1594  0.004 **
#   Residuals 10   1.79029 0.17903         0.8406          
# Total     11   2.12978                 1.0000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


###NirK =======
MM.nirK <- read.csv("otu_table_nirK.csv")
MM.nir.treatment <- MM.nir[,1:16]
MM.nir.spec <- MM.nir[,17:7892]


#Subet genotypes
TvM.nirK=subset(MM.nirK, Genotype == "A632" |Genotype == "B73"
                | Genotype == "LH82"
                | Genotype == "PH207"
                | Genotype == "Mo17"
                | Genotype == "Pa91"
                | Genotype == "Ames21786"
                | Genotype == "Ames21789"
                | Genotype == "Ames21809"
                | Genotype == "PI566674"
                | Genotype == "PI566677"
                | Genotype == "PI566680"
                | Genotype == "Bulk" )
##Subsetting the values of interest 
MM.nirK.treatment <- TvM.nirK[,1:16]
MM.nirK.spec <- TvM.nirK[,17:7892]

#Genotype comparison on th microibome
adonis(MM.nirK.spec ~ Genotype, MM.nirK.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   11     5.499 0.49992  1.1987 0.10881  0.001 ***
#   Residuals 108    45.040 0.41704         0.89119           
# Total     119    50.539                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis(MM.nirK.spec ~ Type, MM.nirK.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type        2     1.296 0.64821  1.5401 0.02565  0.001 ***
#   Residuals 117    49.243 0.42088         0.97435           
# Total     119    50.539                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TvM.nirK.Means<-aggregate(TvM.nirK[, (17:7892)], list(Genotype=TvM.nirK$Genotype, 
                                                     Family=TvM.nirK$Family,
                                                     Type=TvM.nirK$Type), mean)

MM.nirK.treatment <- TvM.nirK.Means[,1:3]
MM.nirK.spec <- TvM.nirK.Means[,4:7879]

##nirK Functional Gene statisitc==== 
 adonis(MM.nirK.spec ~ Type, MM.nirK.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type       2   0.58647 0.29324  1.2644 0.21934  0.001 ***
#   Residuals  9   2.08729 0.23192         0.78066           
# Total     11   2.67376                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.nirK.spec ~ Family, MM.nirK.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Family     1   0.30708 0.30708  1.2975 0.11485  0.001 ***
#   Residuals 10   2.36668 0.23667         0.88515           
# Total     11   2.67376                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


###NOSZ Statistic ====
MM.nosZ <- read.csv("otu_table_nosZ.csv")
MM.nos.treatment <- MM.nosZ[,1:16]
MM.nos.spec <- MM.nosZ[,17:1674]

#Subet genotypes
TvM.nosZ=subset(MM.nosZ, Genotype == "A632" |Genotype == "B73"
                | Genotype == "LH82"
                | Genotype == "PH207"
                | Genotype == "Mo17"
                | Genotype == "Pa91"
                | Genotype == "Ames21786"
                | Genotype == "Ames21789"
                | Genotype == "Ames21809"
                | Genotype == "PI566674"
                | Genotype == "PI566677"
                | Genotype == "PI566680"
                | Genotype == "Bulk" )
##Subsetting the values of interest 
MM.nosZ.treatment <- TvM.nosZ[,1:16]
MM.nosZ.spec <- TvM.nosZ[,17:1674]

#Genotype comparison on th microibome
adonis(MM.nosZ.spec ~ Genotype, MM.nosZ.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   11     6.441 0.58551  1.2899 0.11708  0.001 ***
#   Residuals 107    48.570 0.45392         0.88292           
# Total     118    55.011                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis(MM.nosZ.spec ~ Type, MM.nosZ.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type        2     1.688 0.84384  1.8357 0.03068  0.001 ***
#   Residuals 116    53.323 0.45968         0.96932           
# Total     118    55.011                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TvM.nosZ.Means<-aggregate(TvM.nosZ[, (17:1674)], list(Genotype=TvM.nosZ$Genotype, 
                                                      Family=TvM.nosZ$Family,
                                                      Type=TvM.nosZ$Type), mean)

MM.nosZ.treatment <- TvM.nosZ.Means[,1:3]
MM.nosZ.spec <- TvM.nosZ.Means[,4:1661]

##nosZ Functional Gene statisitc==== 
adonis(MM.nosZ.spec ~ Type, MM.nosZ.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# Type       2    0.8683 0.43417  1.3244 0.22738  0.003 **
#   Residuals  9    2.9505 0.32783         0.77262          
# Total     11    3.8188                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.nosZ.spec ~ Family, MM.nosZ.treatment)

# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)   
# Family     1    0.4995 0.49949  1.5048 0.1308  0.006 **
#   Residuals 10    3.3193 0.33193         0.8692          
# Total     11    3.8188                 1.0000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




###nifH Statistic ====
MM.nifH <- read.csv("otu_table_nifH.csv")
# MM.nif.treatment <- MM.nifH[,1:16]
# MM.nif.spec <- MM.nifH[,17:1674]

#Subet genotypes
TvM.nifH=subset(MM.nifH, Genotype == "A632" |Genotype == "B73"
                | Genotype == "LH82"
                | Genotype == "PH207"
                | Genotype == "Mo17"
                | Genotype == "Pa91"
                | Genotype == "Ames21786"
                | Genotype == "Ames21789"
                | Genotype == "Ames21809"
                | Genotype == "PI566674"
                | Genotype == "PI566677"
                | Genotype == "PI566680"
                #| Genotype == "Bulk" )
##Subsetting the values of interest 
MM.nifH.treatment <- TvM.nifH[,1:13]
MM.nifH.spec <- TvM.nifH[,14:14926]

#Genotype comparison on th microibome
adonis(MM.nifH.spec ~ Genotype, MM.nifH.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   11    2.4439 0.22217  2.2023 0.18321  0.001 ***
# Residuals 108   10.8951 0.10088         0.81679           
# Total     119   13.3390                 1.00000           
# ---

adonis(MM.nifH.spec ~ Type, MM.nifH.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type        2     1.080 0.53998  5.1536 0.08096  0.001 ***
#   Residuals 117    12.259 0.10478         0.91904           
# Total     119    13.339                 1.00000           

adonis(MM.nifH.spec ~ Family, MM.nifH.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Family      1    0.7717 0.77175  7.2463 0.05786  0.001 ***
#   Residuals 118   12.5673 0.10650         0.94214           
# Total     119   13.3390                 1.00000           

TvM.nifH.Means<-aggregate(TvM.nifH[, (14:14926)], list(Genotype=TvM.nifH$Genotype, 
                                                      Family=TvM.nifH$Family,
                                                      Type=TvM.nifH$Type), mean)

MM.nifH.treatment <- TvM.nifH.Means[,1:3]
MM.nifH.spec <- TvM.nifH.Means[,4:14916]

##nifH Functional Gene statisitc==== 
adonis(MM.nifH.spec ~ Type, MM.nifH.treatment)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Type       2   0.13365 0.066823  2.8791 0.39017  0.001 ***
#   Residuals  9   0.20888 0.023209         0.60983           
# Total     11   0.34253                  1.00000           
# ---
adonis(MM.nifH.spec ~ Family, MM.nifH.treatment)

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# Family     1   0.09092 0.090918  3.6134 0.26543  0.004 **
#   Residuals 10   0.25161 0.025161         0.73457          
# Total     11   0.34253                  1.00000          
# ---

###I should probably do the analysis for norB and nosZ





## Quantative PCR Genes Analyis=====
setwd('~/Documents/Projects/Research/Maize Microbiome/MM2015/Data Analysis /CSV folder')
QPCR.All <- read.csv("MM2015.Functional.Gene.Samples.csv")
MM.qPCR.treatment <- QPCR.All[,1:14]
MM.qPCR.spec <- QPCR.All[,15:24]


TvM.qPCR=subset(QPCR.All, Genotype == "A632" |Genotype == "B73"
                | Genotype == "LH82"
                | Genotype == "PH207"
                | Genotype == "Mo17"
                | Genotype == "Pa91"
                | Genotype == "Ames21786"
                | Genotype == "Ames21789"
                | Genotype == "Ames21809"
                | Genotype == "PI566674"
                | Genotype == "PI566677"
                | Genotype == "PI566680"
                #| Genotype == "Bulk"
                )

MM.qPCR.treatment <- TvM.qPCR[,1:14]
MM.qPCR.spec <- TvM.qPCR[,17:24]

adonis(MM.qPCR.spec ~ Genotype, MM.qPCR.treatment)
adonis(MM.qPCR.spec ~ Type, MM.qPCR.treatment)
adonis(MM.qPCR.spec ~ Family, MM.qPCR.treatment)

TvM.aPCR.Means<-aggregate(TvM.qPCR[, (17:24)], list(Genotype=TvM.qPCR$Genotype, 
                                                      Family=TvM.qPCR$Family,
                                                      Type=TvM.qPCR$Type), mean)
MM.qPCR.treatment <- TvM.aPCR.Means[,1:3]
MM.qPCR.spec <- TvM.aPCR.Means[,4:11]

adonis(MM.qPCR.spec ~ Type, MM.qPCR.treatment)
adonis(MM.qPCR.spec ~ Family, MM.qPCR.treatment)
####All funcitonal genes 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type       3   1.07450 0.35817  22.233 0.88111  0.001 ***
#   Residuals  9   0.14498 0.01611         0.11889           
# Total     12   1.21948                 1.00000           

##There are clear differences in composition in functional gene and qPCR abundance in the greenhouse setting! 
##QPCR Gene abundance 
#This ma
ArchAmoAEffectAll<-TvM.qPCR%>%
  group_by(Genotype,Type)%>%
  summarise(mean.Nit=mean(ARCH.AMOA),
            var.Nit=var(ARCH.AMOA),
            Genotypic.effect.N=mean(ARCH.AMOA)-mean(TvM.qPCR$ARCH.AMOA),
            PercentChange=(mean(ARCH.AMOA)-mean(TvM.qPCR$ARCH.AMOA))/mean(TvM.qPCR$ARCH.AMOA)
  )%>%ungroup()

BACTERIAL.AMOAEffectAll<-TvM.qPCR%>%
  group_by(Genotype,Type)%>%
  summarise(mean.Nit=mean(BACTERIAL.AMOA),
            var.Nit=var(BACTERIAL.AMOA),
            Genotypic.effect.N=mean(BACTERIAL.AMOA)-mean(TvM.qPCR$BACTERIAL.AMOA),
            PercentChange=(mean(BACTERIAL.AMOA)-mean(TvM.qPCR$BACTERIAL.AMOA))/mean(TvM.qPCR$BACTERIAL.AMOA)
  )%>%ungroup()

NIFHEffectAll<-TvM.qPCR%>%
  group_by(Genotype,Type)%>%
  summarise(mean.Nit=mean(NIFH),
            var.Nit=var(NIFH),
            Genotypic.effect.N=mean(NIFH)-mean(TvM.qPCR$NIFH),
            PercentChange=(mean(NIFH)-mean(TvM.qPCR$NIFH))/mean(TvM.qPCR$NIFH)
  )%>%ungroup()

NIRKEffectAll<-TvM.qPCR%>%
  group_by(Genotype,Type)%>%
  summarise(mean.Nit=mean(NIRK),
            var.Nit=var(NIRK),
            Genotypic.effect.N=mean(NIRK)-mean(TvM.qPCR$NIRK),
            PercentChange=(mean(NIRK)-mean(TvM.qPCR$NIRK))/mean(TvM.qPCR$NIRK)
  )%>%ungroup()

NIRSEffectAll<-TvM.qPCR%>%
  group_by(Genotype,Type)%>%
  summarise(mean.Nit=mean(NIRS),
            var.Nit=var(NIRS),
            Genotypic.effect.N=mean(NIRS)-mean(TvM.qPCR$NIRS),
            PercentChange=(mean(NIRS)-mean(TvM.qPCR$NIRS))/mean(TvM.qPCR$NIRS)
  )%>%ungroup()

NORBEffectAll<-TvM.qPCR%>%
  group_by(Genotype,Type)%>%
  summarise(mean.Nit=mean(NORB),
            var.Nit=var(NORB),
            Genotypic.effect.N=mean(NORB)-mean(TvM.qPCR$NORB),
            PercentChange=(mean(NORB)-mean(TvM.qPCR$NORB))/mean(TvM.qPCR$NORB)
  )%>%ungroup()

NOSZEffectAll<-TvM.qPCR%>%
  group_by(Genotype,Type)%>%
  summarise(mean.Nit=mean(NOSZ),
            var.Nit=var(NOSZ),
            Genotypic.effect.N=mean(NOSZ)-mean(TvM.qPCR$NOSZ),
            PercentChange=(mean(NOSZ)-mean(TvM.qPCR$NOSZ))/mean(TvM.qPCR$NOSZ)
  )%>%ungroup()

###ASREML-R analysis comparing the QPCR Data=====
library(asreml)

##Load Data=
QPCR.All <- read.csv("MM2015.Functional.Gene.Samples.csv")
MM.qPCR.treatment <- QPCR.All[,1:14]
MM.qPCR.spec <- QPCR.All[,15:24]


TvM.qPCR=subset(QPCR.All, Genotype == "A632" |Genotype == "B73"
                | Genotype == "LH82"
                | Genotype == "PH207"
                | Genotype == "Mo17"
                | Genotype == "Pa91"
                | Genotype == "Ames21786"
                | Genotype == "Ames21789"
                | Genotype == "Ames21809"
                | Genotype == "PI566674"
                | Genotype == "PI566677"
                | Genotype == "PI566680"
                #| Genotype == "Bulk"
)

str(TvM.qPCR)

##nifH statisitcal model
nifH.Model.Random<-asreml(data=TvM.qPCR,
                         fixed=NIFH~Family,
                         random=,
                         na.action = na.method(x = 'fail', y = 'omit'))

wald(nifH.Model.Random)
summary(nifH.Model.Random)
plot(nifH.Model.Random)
plot(varioGram(nifH.Model.Random))
summary(nifH.Model.Random)$bic


#Here we can predict the LS means of the factors
predict(nifH.Model.Random, classify="Family",sed=T)
#These Values can be used in plotting
nifH.LSMeans<-predict(nifH.Model.Random, classify="Family")
nifH.LSMeans<-nifH.LSMeans$pvals

##Nice! This is what I need to plot as a factor in the model
ggplot(nifH.LSMeans, aes(x=reorder(Family,predicted.value), y=predicted.value)) + 
  geom_errorbar(aes(ymin=predicted.value-std.error*1.96, ymax=predicted.value+std.error*1.96), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()


###Looking at the average nifH community -ASREML r====
TvM.aPCR.Means<-aggregate(TvM.qPCR[, (17:24)], list(Genotype=TvM.qPCR$Genotype, 
                                                    Family=TvM.qPCR$Family,
                                                    Type=TvM.qPCR$Type), mean)

nifH.Model.Random<-asreml(data=TvM.aPCR.Means,
                          fixed=NIFH~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(nifH.Model.Random)
summary(nifH.Model.Random)
plot(nifH.Model.Random)
plot(varioGram(nifH.Model.Random))
summary(nifH.Model.Random)$bic
#There is no difference in the requitment of nifH across these time scales- they are almost identiical 

#Bact amoA
qPCR.Model.Random<-asreml(data=TvM.qPCR,
                          fixed=BACTERIAL.AMOA~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(qPCR.Model.Random)
summary(qPCR.Model.Random)

#amoA
qPCR.Model.Random<-asreml(data=TvM.qPCR,
                          fixed=ARCH.AMOA~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(qPCR.Model.Random)
summary(qPCR.Model.Random)


#nirK
qPCR.Model.Random<-asreml(data=TvM.qPCR,
                          fixed=NIRK~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(qPCR.Model.Random)
summary(qPCR.Model.Random)

#nirS
qPCR.Model.Random<-asreml(data=TvM.qPCR,
                          fixed=NIRS~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(qPCR.Model.Random)
summary(qPCR.Model.Random)


#nosZ
qPCR.Model.Random<-asreml(data=TvM.qPCR,
                          fixed=NOSZ~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(qPCR.Model.Random)
summary(qPCR.Model.Random)

#nrfA
qPCR.Model.Random<-asreml(data=TvM.qPCR,
                          fixed=NRFA~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(qPCR.Model.Random)
summary(qPCR.Model.Random)

#norB
qPCR.Model.Random<-asreml(data=TvM.qPCR,
                          fixed=NORB~Family,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))

wald(qPCR.Model.Random)
summary(qPCR.Model.Random)

###Figure 3B. Here Im going to figure out how to replot the Bar-graphs from the paper to exclude norB and nrfA 
library(tidyverse)

TvM.aPCR.Means<-aggregate(TvM.qPCR[, (17:24)], list(Genotype=TvM.qPCR$Genotype, 
                                                    Family=TvM.qPCR$Family,
                                                    Type=TvM.qPCR$Type), mean)
TvM.aPCR.sd<-aggregate(TvM.qPCR[, (17:24)], list(Genotype=TvM.qPCR$Genotype, 
                                                    Family=TvM.qPCR$Family,
                                                    Type=TvM.qPCR$Type), sd)

TvM.aPCR.Means.for<-TvM.aPCR.Means %>% gather(FunctionalGenes, Abundance, ARCH.AMOA:NRFA)
TvM.aPCR.sd.for<-TvM.aPCR.sd %>% gather(FunctionalGenes, Abundance, ARCH.AMOA:NRFA)

TvM.Mean.SD.qPCR<-TvM.qPCR%>%
  group_by(Genotype,Family,Domesticated)%>%
    summarise(#nifH
              NIFH.mean=mean(NIFH),
              NIFH.sd=sd(NIFH),
              #Archamoa
              ARCH.AMOA.mean=mean(ARCH.AMOA),
              ARCH.AMOA.sd=sd(ARCH.AMOA),
              #Bacteriaamoa
              BACTERIAL.AMOA.mean=mean(BACTERIAL.AMOA),
              BACTERIAL.AMOA.sd=sd(BACTERIAL.AMOA),
              #NIRK
              NIRK.mean=mean(NIRK),
              NIRK.sd=sd(NIRK),
              #NIRS
              NIRS.mean=mean(NIRS),
              NIRS.sd=sd(NIRS),
              #NIRS
              NOSZ.mean=mean(NOSZ),
              NOSZ.sd=sd(NOSZ)        
    )%>%ungroup()


TvM.Mean.SD.qPCR<-TvM.qPCR%>%
  group_by(Domesticated)%>%
  summarise(#nifH
    NIFH.mean=mean(NIFH),
    NIFH.sd=sd(NIFH),
    #Archamoa
    ARCH.AMOA.mean=mean(ARCH.AMOA),
    ARCH.AMOA.sd=sd(ARCH.AMOA),
    #Bacteriaamoa
    BACTERIAL.AMOA.mean=mean(BACTERIAL.AMOA),
    BACTERIAL.AMOA.sd=sd(BACTERIAL.AMOA),
    #NIRK
    NIRK.mean=mean(NIRK),
    NIRK.sd=sd(NIRK),
    #NIRS
    NIRS.mean=mean(NIRS),
    NIRS.sd=sd(NIRS),
    #NIRS
    NOSZ.mean=mean(NOSZ),
    NOSZ.sd=sd(NOSZ)        
  )%>%ungroup()

FunctionalGenes.FigureFormat<-TvM.Mean.SD.qPCR %>% gather(FunctionalGenes, Abundance, NIFH.mean:NOSZ.sd)

FunctionalGenes.FigureFormat<-TvM.Mean.SD.qPCR %>% gather(FunctionalGenes, Abundance, -Domesticated, -'**.sd')
##Figure 3A====
library(ggplot2)
library(reshape2)
alpha <- read.csv("~/Documents/Projects/Research/Maize Microbiome/MM2015/Data Analysis /qRT-PCR/Averages/Barplotpqcr.csv")
alpha<-alpha[-c(1,4,7,10,13,16,19,22,17,18,23,24),]
dodge <- position_dodge(width=0.9)

ggplot(data=alpha, aes(gene, mean, fill = type)) +
  geom_bar(stat = "identity", position = "dodge",colour="black") +
  scale_fill_manual(values = c( "grey30","white")) +
  geom_errorbar(aes(ymin = mean - std_err,
                    ymax = mean + std_err), 
                width = 0.25, position = dodge) + xlab("Functional genes ") + ylab("Gene abundance (copies/ng)")+
  theme_bw()#+ theme(legend.position = "none")
  #facet_zoom(ylim = c(0, 250), switch=x)
  #ggtitle("Fuctional Gene Abundance")
  
ggsave(file = "Figure.3A.Final",width = 6,height =4)


##3B Alpha Diversity Bar plots ====

alpha <- read.csv("~/Documents/Projects/Research/Maize Microbiome/MM2015/Data Analysis /CSV folder/Figure.3B.data.csv")
alpha<-subset(alpha,type!="Hybrid")
alpha<-subset(alpha,gene!="16S")
alpha<-subset(alpha,gene!="nrfA")
alpha2<-subset(alpha,gene!="norB")


dodge <- position_dodge(width=0.9)

ggplot(data=alpha2, aes(gene, chao1, fill = type)) +
  geom_bar(stat = "identity", position = "dodge",colour="black") +
  scale_fill_manual(values = c( "grey30","white")) +
  geom_errorbar(aes(ymin = chao1 - chao1_err,
                    ymax = chao1 + chao1_err), 
                width = 0.25, position = dodge) + xlab("Functional genes ") + ylab("Gene richness (chao1)")+
  theme_bw()+facet_zoom(ylim = c(0, 100))#+ theme(legend.position = "none")

# for (i in x){print(as.character(colnames(TvM.qPCR[i])))}
# 
# x<-17:24
# y<-(colnames(TvM.qPCR[i]))
# for (i in x){print(y)}
#                    
# for (i in x){print(c(y,".mean"))}
# 
# colnames(TvM.qPCR[17])
# x<-17:24
# 
# TvM.qPCR[1,1]
# 
# TvM.Mean.SD.qPCR<-TvM.qPCR%>%
#   group_by(Genotype,Family,Domesticated)%>%
#   for (i in x){summarise(TvM.qPCR[1,i].mean=mean(TvM.qPCR[1,i]),
#             TvM.qPCR[1,i].sd=sd(TvM.qPCR[1,i])},
#   )%>%ungroup()
# 
#     
#     
#     for (i in x){summarise(colnames(TvM.qPCR[i]).mean=mean(colnames(TvM.qPCR[i])),
#                            colnames(TvM.qPCR[i]).sd=sd(colnames(TvM.qPCR[i]))}
# 
#     as.character(colnames(TvM.qPCR[1]))
#####Phyloseq TvM Comparison -----
#####Phyloseq Figure 1-2 After this point ====
##Written 1.15.21
library(phyloseq)
library(ggplot2)

##16S Analysis ====
setwd("~/Documents/Projects/Research/Maize Microbiome/MM2015/Phyloseq_MM2015")
#setwd('~/Documents/Projects/Research/Maize Microbiome/MM2015/Data Analysis /CSV folder')

#Laoding Files 
#biom_file = "Data/16S_OTU_NEW.biom"

biom_file = "Data/16S_NEW_table_34000.biom"
map_file = "Data/MappingFile_MaizeGDB.txt"



##Making Phyloseq Object 
biomOTU = import_biom(biom_file,treefilename = 'Data/rep_seq_aligned_pfiltered.tre', 
                      parseFunction = parse_taxonomy_greengenes)
Map = import_qiime_sample_data("Data/MappingFile_MaizeGDB.txt")
OTU16S = merge_phyloseq(biomOTU, Map)

##Susetting Genotypes used in comparison
TvM.OTU<- subset_samples(OTU16S,Genotype == "A632" |Genotype == "B73"
               | Genotype == "LH82"
               | Genotype == "PH207"
               | Genotype == "Mo17"
               | Genotype == "Pa91"
               | Genotype == "Ames21786"
               | Genotype == "Ames21789"
               | Genotype == "Ames21809"
               | Genotype == "PI566674"
               | Genotype == "PI566677"
               | Genotype == "PI566680"
               #| Genotype == "Bulk" 
               )

TvM.OTU<-transform_sample_counts(TvM.OTU, log1p)
##Set theme===
theme_set(theme_classic())
#Ploting Alpha Diversity====
p2<-plot_richness(TvM.OTU, x = "Family", measures = "Chao1") +
     theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()+
  theme(strip.text.x = element_blank())+   xlab("Treatment") + ylab("Alpha diversity (chao1)")+
  ggtitle("Prokaryotic 16S rRNA")

##NMDS Figure 1=====
GP.ord <- ordinate(TvM.OTU, "NMDS", "bray",)
p1 = plot_ordination(TvM.OTU, GP.ord, type="Samples", color="Family",  title="Bacterial 16S Ordination")+
  stat_ellipse(aes(color=Family, group=Family),type = "norm")+
  geom_point(size=3)+scale_fill_brewer(palette="BuPu")+
  geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)
  
print(p1)

#Here we re removing some outliers #0337, 0338
to_remove <- c("MM2015-0337", "MM2015-0338")

z <- prune_samples(!(sample_names(TvM.OTU) %in% to_remove), TvM.OTU)

GP.ord <- ordinate(z, "NMDS", "bray",)
p1 = plot_ordination(z, GP.ord, type="Samples", color="Family",  title="Prokaryotic 16S rRNA")+
  geom_point(size=3)+scale_fill_brewer(palette="BuPu")+
  labs(col="Treatment")+  annotate("text", x=-.2, y=-.15, label= "R2=0.27, p<0.001")

print(p1)


##Plotting the two figure together
par(mfrow=c(1,2))
plot(p1)
plot(p2)

plot(p2,p1)

##Large differences in the microbiome compostion between maize and teosinte in the greenhouse setting. This treatment factor explains about 80% of the varaiton 

##NMDS Average Figure 1====
TvM.OTU <- subset_samples(TvM.OTU, Family != "Bulk")

OTU_merge = merge_samples(TvM.OTU, "Genotype", fun = mean) 
sample_data(OTU_merge)$Genotype <- levels(sample_data(TvM.OTU)$Genotype)
sample_data(OTU_merge)$Type <- levels(sample_data(TvM.OTU)$Type)

GenotypeAve = transform_sample_counts(OTU_merge, function(x)  x/10 )

GP.ord <- ordinate(GenotypeAve, "NMDS", "bray",)
p1 = plot_ordination(GenotypeAve, GP.ord, type="Samples", color="Family",  title="Bacterial 16S Ordination")+
  geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)
print(p1)

##It seems like the Typing lable is off for some reason- but the Family typing is not off. 
##Im going to remove the NA from data
TvM.OTU = subset_taxa(TvM.OTU, Phylum != "NA")


##DESEQ2 Analysis determining which taxa are different. Prokaryotes====
library("DESeq2")
diagdds = phyloseq_to_deseq2(TvM.OTU, ~ Family)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TvM.OTU)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#[1] 346  19
##The first number is the number of taxa that are present in the 
#These are the taxa that were changing in response to treatment 

##Grab code and ordering form other data set to make sure Im doing it right
#write.csv(sigtab,'DESEQ2.0.01.Results.csv')

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

#Order
x= tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

##Plot ##Need to remove the NA from files 
##Manuscript Figure 2====
ggplot(sigtab, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

##Here we recaptilualed the first plot. I think Angela wants me to use this model for the stack plot
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#######
##Here Is the code that I should use to 
MM15.Average = merge_samples(TvM.OTU, "Family", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(MM15.Average)$Family <- levels(sample_data(TvM.OTU)$Family)
#relative
MM15.Average = transform_sample_counts(MM15.Average, function(x) 100 * x/sum(x) )

##Here we will be filterign out these enriched taxa for stack plot
TopNOTUs<-c(rownames(sigtab))

#Here we are getting the enriched list of taxa to subset them- them 
MM20156x610 = prune_taxa(TopNOTUs, MM15.Average)
##This commented out code will be used if I want to plot the Raw abudannce.
#MM20156x610 = prune_taxa(TopNOTUs, TvM.OTU)

MM20156x610 = subset_taxa(MM20156x610, Order != "NA")


##I orginally plotted them in abudnace now I will plot them in Relvativee abudnance
##Here I am pulling the signifcant OTUS and using them in a stack plot
plot_bar(MM20156x610, "sample_Family", fill = "Class", facet_grid = ~Phylum)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(legend.position = "none",strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_grey()


plot_bar(MM20156x610, "sample_Family", fill = "Class", facet_grid = ~Phylum)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_grey()


##Figure alterations Angela requested 
#MM20156x6_Man = subset_taxa(MM20156x610, Phylum != "Armatimonadetes"|Phylum != "Chloroflexi"|Phylum != "Crenarchaeota"|Phylum != "Fibrobacteres" |Phylum != "Gemmatimonadetes"|Phylum != "Planctomycetes")
MM20156x6_Man = subset_taxa(MM20156x610, Phylum == "Acidobacteria"|Phylum == "Actinobacteria"|Phylum == "Bacteroidetes"|Phylum == "Firmicutes" |Phylum == "Proteobacteria"|Phylum == "Verrucomicrobia")

plot_bar(MM20156x6_Man, "sample_Family", fill = "Class", facet_grid = ~Phylum)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(legend.position = "none",strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+scale_fill_grey()

plot_bar(MM20156x6_Man, "sample_Family", fill = "Class", facet_grid = ~Phylum)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+scale_fill_grey()



##Major changes in the proteobacteria- these may disagree with the field findings but Im not super sure.
##Im going to pull out the fermicutes as those seem like an interesting grouping
##Teosinte showed increases in Firmicutes, Verrucomicrobia, and Actinobacteria
##Inbred showed increases in Proteobacteria, Bacteroidetes 
##The rest were really minor in changes in the commuminty 

get_taxa_unique(MM20156x610, "Phylum")
get_taxa_unique(MM20156x610, "Order")
get_taxa_unique(MM20156x610, "Class")
get_taxa_unique(MM20156x610, "Genus")

###Supplemental Composition: Enriched in Teosinte
Firmicutes = subset_taxa(MM20156x610, Phylum == "Firmicutes")
#King, Phyl, Order, Family, Genus, Sp
plot_bar(Firmicutes, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment",title="Firmicutes")+
  scale_y_continuous(labels = scales::percent_format(scale=1))
##Seems like this enriches Paenibacillus which has been shown to be a plant promoting bacteria

Verrucomicrobia = subset_taxa(MM20156x610, Phylum == "Verrucomicrobia")
#King, Phyl, Order, Family, Genus, Sp
plot_bar(Verrucomicrobia, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment",title="Verrucomicrobia")+
  scale_y_continuous(labels = scales::percent_format(scale=1))
##Here we see increases in Chthoniobacterales

Actinobacteria = subset_taxa(MM20156x610, Phylum == "Actinobacteria")
#King, Phyl, Order, Family, Genus, Sp
plot_bar(Actinobacteria, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment",title="Actinobacteria")+
  scale_y_continuous(labels = scales::percent_format(scale=1))
##Actinomycetales is larger in the teosinte Rhizosphere this was a pervious finging
##Gaiellales is also present ieosinte rhizosphere. these are interesting fingd.

# Actinobacteria = subset_taxa(MM20156x610, Phylum == "Actinobacteria")
# #King, Phyl, Order, Family, Genus, Sp
# plot_bar(Actinobacteria, "sample_Family", fill = "Family", facet_grid = ~Order)+
#   geom_bar(stat="identity")+
#   scale_fill_hue()+
#   labs(y="Relative abundance", x="Treatment")+
#   scale_y_continuous(labels = scales::percent_format(scale=1))
# ##Now I will be plotting the taxa that were increased by Inbred maize rhizosphere

#####Enriched in maize Rhizosphere compared to teosinte


Proteobacteria = subset_taxa(MM20156x610, Phylum == "Proteobacteria")
#King, Phyl, Order, Family, Genus, Sp
plot_bar(Proteobacteria, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment",title="Proteobacteria")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(legend.position = "none",strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_grey()


plot_bar(Proteobacteria, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment",title="Proteobacteria")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_grey()

## In most casses inbred maize showed to have larger differences in the micrbiome 
##Interestingly Alteromonadaceae was enriched in the teosinte rhizosphere. 
Bacteroidetes = subset_taxa(MM20156x610, Phylum == "Bacteroidetes")
#King, Phyl, Order, Class, Family, Genus, Sp
plot_bar(Bacteroidetes, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment",title="Bacteroidetes")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  facet_wrap(vars(Order),nrow=1)+theme(strip.text.x = element_text(size = 6))


##Just for Fun Im going to look at nitrosomonas babe|| This didnt work because the data didnt show signifcant differences between the two 
Nitrfiers = subset_taxa(TvM.OTU, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")
Nitrfiers = subset_taxa(MM15.Average, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")
Nitrfiers = subset_taxa(TvM.OTU, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")

#King, Phyl, Order, Class, Family, Genus, Sp
plot_bar(Nitrfiers, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))

##The Nitrosophaerales were enriched in the teosinte||

plot_bar(Nitrfiers, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Abundance", x="Treatment")

#This is kinda of interesting but these taxa make up a really small percentage of the MC

#Ploting Alpha Diversity====
plot_richness(Nitrfiers, x = "Family", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

##NMDS Figure 1=====
GP.ord <- ordinate(Nitrfiers, "NMDS", "bray",)
p1 = plot_ordination(TvM.OTU, GP.ord, type="Samples", color="Family",  title="Nitrifer 16S Ordination")+
  stat_ellipse(aes(color=Family, group=Family),type = "norm")+
  geom_point(size=3)#+
geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)

print(p1)
##Differences in the microbial communites are within nitrosomonadales and nitrosophaerales
MM15.Average.Type = merge_samples(TvM.OTU, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(MM15.Average.Type)$Type <- levels(sample_data(TvM.OTU)$Type)
#relative
MM15.Average.Type = transform_sample_counts(MM15.Average.Type, function(x) 100 * x/sum(x) )

Nitrfiers = subset_taxa(MM15.Average.Type, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")

plot_bar(Nitrfiers, "Type", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))


##Here we are going to look at the preditory bacteria 16S=====

MM15.Average
get_taxa_unique(MM15.Average, "Order")
get_taxa_unique(MM15.Average, "Class")
get_taxa_unique(MM15.Average, "Phylum")
get_taxa_unique(MM15.Average, "Genus")
get_taxa_unique(MM15.Average, "Family")

#Myxococcales
Deltaproteobacteria = subset_taxa(MM15.Average, Class == "Deltaproteobacteria")

get_taxa_unique(Deltaproteobacteria, "Family")

plot_bar(Deltaproteobacteria, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))


Gammaproteobacteria = subset_taxa(MM15.Average, Class == "Gammaproteobacteria")

plot_bar(Gammaproteobacteria, "sample_Family", fill = "Family", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))


get_taxa_unique(Gammaproteobacteria, "Family")


Preditors = subset_taxa(OTU16S, Genus == "Lysobacter"|Family == "Myxococcaceae" |Family == "Bdellovibrionaceae")

plot_bar(Preditors, "Genotype", fill = "Genus", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))

MM15.Average = merge_samples(OTU16S, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(MM15.Average)$Type <- levels(sample_data(OTU16S)$Type)
#relative
MM15.Average = transform_sample_counts(MM15.Average, function(x) 100 * x/sum(x) )



Preditors = subset_taxa(MM15.Average, Genus == "Lysobacter"|Family == "Myxococcaceae" |Family == "Bdellovibrionaceae")

plot_bar(Preditors, "Type", fill = "Genus", facet_grid = ~Order)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))

Preditors = subset_taxa(TvM.OTU, Genus == "Lysobacter"|Family == "Myxococcaceae" |Family == "Bdellovibrionaceae")




#Fungal ITS Anlaysis ====
# import data into R studio
#biom_file = import_biom("Data/otu_table_ITS.biom")
biom_file = import_biom("Data/otu_table_AllOTUN.biom")
##May have renamed to all Fungal OTU BS
Map = import_qiime_sample_data("Data/MappingFile_MaizeGDB.txt")

OTUITS = merge_phyloseq(biom_file, Map)

##Susetting Genotypes used in comparison
TvM.ITS.OTU<- subset_samples(OTUITS,Genotype == "A632" |Genotype == "B73"
                         | Genotype == "LH82"
                         | Genotype == "PH207"
                         | Genotype == "Mo17"
                         | Genotype == "Pa91"
                         | Genotype == "Ames21786"
                         | Genotype == "Ames21789"
                         | Genotype == "Ames21809"
                         | Genotype == "PI566674"
                         | Genotype == "PI566677"
                         | Genotype == "PI566680"
                         | Genotype == "Bulk" 
)

TvM.ITS.OTU<-transform_sample_counts(TvM.ITS.OTU, log1p)

##Set theme===
theme_set(theme_classic())
#Ploting Alpha Diversity====
plot_richness(TvM.ITS.OTU, x = "Family", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()+theme(strip.text.x = element_blank())+   xlab("Treatment") + 
  ylab("Alpha diversity (chao1)")+ggtitle("Fungal ITS")

##NMDS Figure 1 Fungi=====
GP.ord <- ordinate(TvM.ITS.OTU, "NMDS", "bray",)
p1 = plot_ordination(TvM.ITS.OTU, GP.ord, type="Samples", color="Family",  title="Fungal ITS")+
  stat_ellipse(aes(color=Family, group=Family),type = "norm")+
  geom_point(size=3)+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1)
plot(p1)


theme_set(theme_bw())

##Manuscript Figure 1B=====
GP.ord <- ordinate(TvM.ITS.OTU, "NMDS", "bray",)
p1 = plot_ordination(TvM.ITS.OTU, GP.ord, type="Samples", color="Family",  title="Fungal ITS")+
  geom_point(size=3)+scale_fill_brewer(palette="BuPu")+
  labs(col="Treatment")+  annotate("text", x=-1, y=-1, label= "R2=0.21, p<0.004")

print(p1)

TvM.ITS.OTU = subset_taxa(TvM.ITS.OTU, Rank2 != "p__unidentified")

TvM.ITS.OTU = subset_samples(TvM.ITS.OTU, Genotype != "Bulk")

library("DESeq2")

TvM.ITS.OTU.1 = transform_sample_counts(TvM.ITS.OTU, function(x) x+1 )
diagdds = phyloseq_to_deseq2(TvM.ITS.OTU.1, ~ Family)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TvM.ITS.OTU.1)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)
#1] 48 13
##The first number is the number of taxa that are present in the 

#write.csv(sigtab,'Fungi.DESEQ2.0.01.Results.csv')


#Order=Rank6
x= tapply(sigtab$log2FoldChange, sigtab$Rank3, function(x) max(x))
x = sort(x, TRUE)
sigtab$Rank3 = factor(as.character(sigtab$Rank3), levels=names(x))

##Figures|| I can make these figures Need to reorder
ggplot(sigtab, aes(x=Rank3, y=log2FoldChange, color=Rank2)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  labs(x="Class", color="Phylum",title="Fungal ITS")


####Fungal Stack plots

MM15.ITS.Average = merge_samples(TvM.ITS.OTU, "Family", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(MM15.ITS.Average)$Family <- levels(sample_data(TvM.ITS.OTU)$Family)
#relative
MM15.ITS.Average = transform_sample_counts(MM15.ITS.Average, function(x) 100 * x/sum(x) )

TopNOTUs<-c(rownames(sigtab))

#Here we are getting the enriched list of taxa to subset them- them 
MM2015.ITS = prune_taxa(TopNOTUs, MM15.ITS.Average)
##This commented out code will be used if I want to plot the Raw abudannce.
#MM20156x610 = prune_taxa(TopNOTUs, TvM.OTU)

MM20156x610 = subset_taxa(MM20156x610, Order != "NA")


##I orginally plotted them in abudnace now I will plot them in Relvativee abudnance
##Here I am pulling the signifcant OTUS and using them in a stack plot
##Fungal Stack Plots====
plot_bar(MM2015.ITS, "Family", fill = "Rank3", facet_grid = ~Rank2)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(legend.position = "none",strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_grey()
##The scale doesn look right here but im not sure.
plot_bar(MM2015.ITS, "Family", fill = "Rank3", facet_grid = ~Rank2)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_grey()

###Here Imgoing to put the reduced ones Angela Requested
MM2015.ITS_Man = subset_taxa(MM2015.ITS, Rank2 == "p__Ascomycota"|Rank2 == "p__Basidiomycota"|Rank2 == "p__Rozellomycota"|Rank2 == "p__Zygomycota")

MM2015.ITS@tax_table


plot_bar(MM2015.ITS_Man, "Family", fill = "Rank3", facet_grid = ~Rank2)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()+theme(legend.position = "none",strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+scale_fill_grey()


plot_bar(MM2015.ITS_Man, "Family", fill = "Rank3", facet_grid = ~Rank2)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  scale_y_continuous(labels = scales::percent_format(scale=1))+
  theme_bw()#+theme(legend.position = "none",strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+scale_fill_grey()

##Phyloseq of Bacterial amoA Functional Genes====
biom_file = import_biom("Data/otu_table_BamoA.biom")
Map = import_qiime_sample_data("Data/MappingFile_MaizeGDB.txt")

OTU_AMOB = merge_phyloseq(biom_file, Map)


TvM.AMOA.OTU<- subset_samples(OTU_AMOB,Genotype == "A632" |Genotype == "B73"
                             | Genotype == "LH82"
                             | Genotype == "PH207"
                             | Genotype == "Mo17"
                             | Genotype == "Pa91"
                             | Genotype == "Ames21786"
                             | Genotype == "Ames21789"
                             | Genotype == "Ames21809"
                             | Genotype == "PI566674"
                             | Genotype == "PI566677"
                             | Genotype == "PI566680")
                             #| Genotype == "Bulk" )
  

plot_richness(TvM.AMOA.OTU, x = "Family", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

                             
plot_bar(TvM.AMOA.OTU, "Family", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_bar(TvM.AMOA.OTU, "Type", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



MM15.Average.Type = merge_samples(TvM.AMOA.OTU, "Family", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(MM15.Average.Type)$Family <- levels(sample_data(TvM.AMOA.OTU)$Family)
#relative
TvM.AMOA.OTU.Ave = transform_sample_counts(MM15.Average.Type, function(x) 100 * x/sum(x) )


##Thing Might be a figure I would want to present==== 
plot_bar(TvM.AMOA.OTU.Ave, "Type", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(title="Genus"))+title(main=)


##NMDS amoA=====
GP.ord <- ordinate(TvM.AMOA.OTU, "NMDS", "bray",)
p1 = plot_ordination(TvM.AMOA.OTU, GP.ord, type="Samples", color="Family",  title="BactamoA Gene")+
  stat_ellipse(aes(color=Family, group=Family),type = "norm")+
  geom_point(size=3)#+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1
plot(p1)

library("DESeq2")

TvM.AMOA.OTU.1 = transform_sample_counts(TvM.AMOA.OTU, function(x) x+1 )
diagdds = phyloseq_to_deseq2(TvM.AMOA.OTU.1, ~ Family)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05 #### Results look a little different with an alpha of 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TvM.AMOA.OTU.1)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)


TopNOTUs<-c(rownames(sigtab))

#Here we are getting the enriched list of taxa to subset them
TvM.AMOA.OTU.2 = prune_taxa(TopNOTUs, TvM.AMOA.OTU)
##This commented out code will be used if I want to plot the Raw abudannce.
#MM20156x610 = prune_taxa(TopNOTUs, TvM.OTU)

##Supplemental Figure results Bac amoA
plot_bar(TvM.AMOA.OTU.2, "Family", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative Abundance", x="Treatment",title="Bacterial amoA DESEQ OTUs")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  guides(fill=guide_legend(title="Species"))

plot_bar(TvM.AMOA.OTU, "Family", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative Abundance", x="Treatment",title="Bacterial amoA")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  guides(fill=guide_legend(title="Species"))


plot_bar(TvM.AMOA.OTU.2, "Type", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


plot_bar(TvM.AMOA.OTU.2.Ave, "Type", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


TvM.AMOA.OTU.1 = transform_sample_counts(TvM.AMOA.OTU.2, function(x) x+1 )

GP.ord <- ordinate(TvM.AMOA.OTU.1, "NMDS", "bray",)
p1 = plot_ordination(TvM.AMOA.OTU.2, GP.ord, type="Samples", color="Type",  title="BactamoA Gene")+
  stat_ellipse(aes(color=Family, group=Family),type = "norm")+
  geom_point(size=3)#+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1
plot(p1)

##Archeal amoA =====

biom_file = import_biom("Data/otu_table_AamoA.biom")
Map = import_qiime_sample_data("Data/MappingFile_MaizeGDB.txt")

OTU_AMOA = merge_phyloseq(biom_file, Map)


TvM.AMOA.OTU<- subset_samples(OTU_AMOA,Genotype == "A632" |Genotype == "B73"
                              | Genotype == "LH82"
                              | Genotype == "PH207"
                              | Genotype == "Mo17"
                              | Genotype == "Pa91"
                              | Genotype == "Ames21786"
                              | Genotype == "Ames21789"
                              | Genotype == "Ames21809"
                              | Genotype == "PI566674"
                              | Genotype == "PI566677"
                              | Genotype == "PI566680")
#| Genotype == "Bulk" )


plot_richness(TvM.AMOA.OTU, x = "Family", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()


plot_bar(TvM.AMOA.OTU, "Family", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

GP.ord <- ordinate(TvM.AMOA.OTU, "NMDS", "bray",)
p1 = plot_ordination(TvM.AMOA.OTU, GP.ord, type="Samples", color="Family",  title="BactamoA Gene")+
  stat_ellipse(aes(color=Family, group=Family),type = "norm")+
  geom_point(size=3)#+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1
plot(p1)

library("DESeq2")

TvM.AMOA.OTU.1 = transform_sample_counts(TvM.AMOA.OTU, function(x) x+1 )
diagdds = phyloseq_to_deseq2(TvM.AMOA.OTU.1, ~ Family)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TvM.AMOA.OTU.1)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)


TopNOTUs<-c(rownames(sigtab))

#Here we are getting the enriched list of taxa to subset them- them 
TvM.AMOA.OTU.2 = prune_taxa(TopNOTUs, TvM.AMOA.OTU)
##This commented out code will be used if I want to plot the Raw abudannce.
#MM20156x610 = prune_taxa(TopNOTUs, TvM.OTU)

plot_bar(TvM.AMOA.OTU.2, "Family", fill = "Rank7", facet_grid=~Rank4)+
  geom_bar(stat="identity")+
  scale_fill_hue()+
  labs(y="Relative abundance", x="Treatment")+
  theme_bw()+theme(strip.text.x = element_text(angle=90),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





###Below Here I will be listing all of the information for the Simple Regressions I perfromed between
##The functional gene data and diversity
library(vegan)
library(rlm)
library(Hmisc)
setwd('~/Documents/Projects/Research/Maize Microbiome/MM2015/Data Analysis /CSV folder')
#Loading Data 
#These Data consist of sampkes 
FunctionGeneAverages<- read.csv("MM2015.6x6.ManuscriptRegression.csv")
#FunctionGeneAverages<- read.csv("PhenotypesCorrelation2.csv")


CorrTable<-cor(FunctionGeneAverages[,2:94])
#Calcualted
#After adding the time factor to this model it seemsclear the year has a stong effect on root traits
#We see a stong negative relationship on AI, UMLW, what root traits are these

cor.test(FunctionGeneAverages$FGD1, FunctionGeneAverages$Average.Functional.Gene.Abundance)
# t = -2.1811, df = 9, p-value = 0.05707
# cor 
# -0.5880401 
#Missing the NorB values
cor.test(FunctionGeneAverages$FGD2, FunctionGeneAverages$Average.Functional.Gene.Abundance)
# t = -2.2956, df = 9, p-value = 0.04734
# cor 
# -0.6076948 
# This was the averages that were perviously calculated
cor.test(FunctionGeneAverages$Functional.Gene.Diversity, FunctionGeneAverages$Average.Functional.Gene.Abundance)
# t = -3.1121, df = 9, p-value = 0.01248
#   cor 
# -0.7199521 

#All of these correlations are significant and strong. I'm regressing the Genotype's avarage gene richenss against average functional gene abundance 

ggplot((FunctionGeneAverages$Functional.Gene.Diversity~FunctionGeneAverages$Average.Functional.Gene.Abundance), main="Scat")


plot(FunctionGeneAverages, aes(x=Average.Functional.Gene.Abundance, y=FGD1)) +
  geom_smooth(method=lm) +geom_point(shape=1)


ggplot(FunctionGeneAverages, aes(x=Average.Functional.Gene.Abundance, y=FGD2)) +
  geom_smooth(method=lm) +geom_point(shape=1)

###Need to look at the PERMAOVA for the Diversity Data====
AbuDiv<-read.csv("MM2015.6x6.Diversity.QPCRAbundance.csv")

FGD.treatment <- AbuDiv[,1:2]
FGD.spec <- AbuDiv[,5:10]

FGD.spec.omit<-na.omit(FGD.spec)
#Genotype comparison on th microibome
adonis(FGD.spec ~ Type, FGD.treatment)

##I MAY NEED TO CORRECT THIS 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# Type       1  0.091678 0.091678  13.703 0.60357  0.003 **
#   Residuals  9  0.060215 0.006691         0.39643          
# Total     10  0.151892                  1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Manuscript Figure S5####
ggplot(FunctionGeneAverages, aes(x=Average.Functional.Gene.Abundance, y=Functional.Gene.Diversity)) +
  geom_smooth(method=lm, color="black")+
  geom_point(aes(colour = factor(Type)))+theme_classic()+
  labs( y="Average Functional Gene Diversity", title="Genotypes Diversity-Abundance Relationship", x="Average Functional Gene Abundance")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(col="Treatment")+annotate("text", x=300, y=50, label= "r=0.731, p<0.01")


ggplot(FunctionGeneAverages, aes(x=nosZ.qpcr, y=nosZ)) +
  geom_smooth(method=lm, color="black")+
  geom_point(aes(colour = factor(Type)))+theme_classic()+
  labs( y="Average Functional Gene Diversity (chao1)", title="nosZ Diversity-Abundance Relationship", x="Average Functional Gene Abundance (copies/ng)")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(col="Treatment")+annotate("text", x=150, y=10, label= "r=0.3423839, p=0.3027")
cor.test(FunctionGeneAverages$nosZ.qpcr, FunctionGeneAverages$nosZ)
# data:  FunctionGeneAverages$nosZ.qpcr and FunctionGeneAverages$nosZ
# t = 1.0932, df = 9, p-value = 0.3027
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3240461  0.7817062
# sample estimates:
#   cor 
# 0.3423839 


ggplot(FunctionGeneAverages, aes(x=nirS.qpcr, y=nirS)) +
  geom_smooth(method=lm, color="black")+
  geom_point(aes(colour = factor(Type)))+theme_classic()+
  labs( y="Average Functional Gene Diversity (chao1)", title="nirS Diversity-Abundance Relationship", x="Average Functional Gene Abundance (copies/ng)")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(col="Treatment")+annotate("text", x=60, y=10, label= "r=0.4531642, p=0.139")
cor.test(FunctionGeneAverages$nirS.qpcr, FunctionGeneAverages$nirS)
# data:  FunctionGeneAverages$nirS.qpcr and FunctionGeneAverages$nirS
# t = 1.6076, df = 10, p-value = 0.139
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1631745  0.8150852
# sample estimates:
#   cor 
# 0.4531642 

ggplot(FunctionGeneAverages, aes(x=nirK.qpcr, y=nirK)) +
  geom_smooth(method=lm, color="black")+
  geom_point(aes(colour = factor(Type)))+theme_classic()+
  labs( y="Average Functional Gene Diversity (chao1)", title="nirK Diversity-Abundance Relationship", x="Average Functional Gene Abundance (copies/ng)")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(col="Treatment")+annotate("text", x=500, y=270, label= "r=-0.2412213, p=0.4501")
cor.test(FunctionGeneAverages$nirK.qpcr, FunctionGeneAverages$nirK)
# data:  FunctionGeneAverages$nirK.qpcr and FunctionGeneAverages$nirK
# t = -0.78602, df = 10, p-value = 0.4501
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.7160016  0.3861359
# sample estimates:
#   cor 
# -0.2412213


ggplot(FunctionGeneAverages, aes(x=AamoA.qpcr, y=amoA)) +
  geom_smooth(method=lm, color="black")+
  geom_point(aes(colour = factor(Type)))+theme_classic()+
  labs( y="Average Functional Gene Diversity (chao1)", title="Archeal amoA Diversity-Abundance Relationship", x="Average Functional Gene Abundance (copies/ng)")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(col="Treatment")+annotate("text", x=800, y=19, label= "r=-0.3357908, p=0.2859")
cor.test(FunctionGeneAverages$AamoA.qpcr, FunctionGeneAverages$amoA)
# t = -1.1273, df = 10, p-value = 0.2859
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.7627099  0.2949511
# sample estimates:
#   cor 
# -0.3357908 

ggplot(FunctionGeneAverages, aes(x=BamoA.qpcr, y=amoB)) +
  geom_smooth(method=lm, color="black")+
  geom_point(aes(colour = factor(Type)))+theme_classic()+
  labs( y="Average Functional Gene Diversity (chao1)", title="Bacterial amoA Diversity-Abundance Relationship", x="Average Functional Gene Abundance (copies/ng)")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(col="Treatment")+annotate("text", x=400, y=8, label= "r=0.4895889, p=0.1062")
cor.test(FunctionGeneAverages$BamoA.qpcr, FunctionGeneAverages$amoB)
# data:  FunctionGeneAverages$BamoA.qpcr and FunctionGeneAverages$amoB
# t = -1.7756, df = 10, p-value = 0.1062
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.830219  0.117260
# sample estimates:
#   cor 
# -0.4895889

ggplot(FunctionGeneAverages, aes(x=nifH.qpcr, y=nifH)) +
  geom_smooth(method=lm, color="black")+
  geom_point(aes(colour = factor(Type)))+theme_classic()+
  labs( y="Average Functional Gene Diversity (chao1)", title="nifH Diversity-Abundance Relationship", x="Average Functional Gene Abundance (copies/ng)")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(col="Treatment")+annotate("text", x=500, y=10, label= "r=0.4477685, p=0.1444")
cor.test(FunctionGeneAverages$nifH.qpcr, FunctionGeneAverages$nifH)
# t = 1.5836, df = 10, p-value = 0.1444
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1697562  0.8128006
# sample estimates:
#   cor 
# 0.4477685 

cor()

##ASREML -r === Checking the Regressions
library(asreml)


Function.Model<-asreml(data=FunctionGeneAverages,
                          fixed=Average.Functional.Gene.Abundance~Functional.Gene.Diversity,
                          random=,
                          na.action = na.method(x = 'fail', y = 'omit'))


wald(Function.Model)
summary(Function.Model)
plot(Function.Model)
plot(varioGram(Function.Model))
summary(Function.Model)$bic
