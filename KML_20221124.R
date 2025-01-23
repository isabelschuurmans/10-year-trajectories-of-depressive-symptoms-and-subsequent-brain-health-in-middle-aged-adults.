######################################
# DIFFERENTIATING THE MATERNAL BRAIN #
######################################

# Final MAIN ANALYSIS

# Authors: I.K. Schuurmans, ..., A. Luik. 
# 2020-03-18

#----------------------------------------------------------
# LIBRARY
#----------------------------------------------------------

# Open library
library(foreign)
library(magrittr)
library(dplyr)
library(kml)
library(ggplot2)
library(xlsx)
library(mice)
library(haven)

#----------------------------------------------------------
# FUNCTIONS
#----------------------------------------------------------

source('O:/medewerkers/042647 Schuurmans, I/Project_3_DEP_MRI/R Code/Function_20220322.r')

#----------------------------------------------------------
# READ IN DATA
#----------------------------------------------------------

setwd("O:/Medewerkers/042647 Schuurmans, I/Project_3_DEP_MRI")
df <- readRDS('data_depr_project.rds')

#----------------------------------------------------------
# recode variables
#----------------------------------------------------------

# wml recoden
df$wml_log <- log(df$wm_hypointensities_vol)
df$wml_log <- as.numeric(scale(ifelse(df$wml_log == -Inf, -6, df$wml_log)))

#----------------------------------------------------------
# KML
#----------------------------------------------------------

# GET CLUSTERS

# define dataset
ldl <- clusterLongData(traj = df[,c('mp', 'm3', 'm9')], #trajectory depression
                       idAll = as.character(df[,'id_dicom'])) #identifier
plot(ldl)

# define trajectories
set.seed(12345)
#kml(ldl)
#choice(ldl)

# plot the different trajectories and criterion
#plotMeans(ldl, 4)
#plotAllCriterion(ldl)

# get the clusters
#df$clus <- getClusters(ldl, 4)

# save
#saveRDS(df, 'df.rds')
#df <- readRDS("df.rds") #--> does not include all ages
df <- readRDS("df_1.rds")

#----------------------------------------------------------
# recode variables
#----------------------------------------------------------

# dichotomize depression # women cutoff > .80; men cutoff > .71

## pregnancy
df$mp_d <- as.factor(ifelse(df$type == 2,
                            ifelse(df$mp > .8, 1, 0),
                            ifelse(df$mp > .71, 1, 0)))

## age 3
df$m3_d <- as.factor(ifelse(df$type == 2,
                            ifelse(df$m3 > .8, 1, 0),
                            ifelse(df$m3 > .71, 1, 0)))

## age 9
df$m9_d <- as.factor(ifelse(df$type == 2,
                            ifelse(df$m9 > .8, 1, 0),
                            ifelse(df$m9 > .71, 1, 0)))

## age 13
df$cesd_d <- as.factor(ifelse(df$depression > 16, 1, 0))

#----------------------------------------------------------
# GET INFORMATION TABLE 1
#----------------------------------------------------------

# recode
df$edu <- ifelse(df$edu_high == 0 & df$edu_intermediate == 0, 1, 0)
df$mardich <- df$mardich -1

# prepare df
df_table1 <- df[,c('agepre_1', 'type', 'ethni', 'edu', 'edu_intermediate', 'edu_high','mardich', 'income_dich', 'mp', 'm3', 'm9', 'depression', 'mp_d', 'm3_d', 'm9_d', 'cesd_d')]
df_table1$type <- as.factor(ifelse(df_table1$type==2, 1, 0))
df_table1[,c('ethni', 'edu',"edu_intermediate","edu_high", 'mardich', 'income_dich')] <- apply(df_table1[,c('ethni', 'edu',"edu_intermediate","edu_high", 'mardich', 'income_dich')], 2, as.factor)

# get info
table1 <- data.frame(get_table1(df_table1[df$clus=="A",]),
                     get_table1(df_table1[df$clus=="B",]),
                     get_table1(df_table1[df$clus=="C",]),
                     get_table1(df_table1[df$clus=="D",]))
rownames(table1) <- c('N',colnames(df_table1))
colnames(table1) <- rep(c('low symptoms', 'increasing symptoms', 'decreasing symptoms', 'always high'))

# write out
write.xlsx(table1, 'Tables/Table1.xlsx')  

# set edu_low to NULL
df$edu <- NULL

#----------------------------------------------------------
# PREPARE AND EXECUTE MULTIPLE IMPUTATION
#----------------------------------------------------------

# select only variables you need
df_imp <- df[,c("mother","partner","id_dicom","type","ageparent","totalgrayvol",
                "cerebralwhitemattervol","amygdala","hippocampus","thalamus", "caudate",                                           
                "putamen","pallidum","accumbens", "etiv","wml_log", "li","cmb","mp",
                "m3","m9" ,"smoking","ethni","edu_intermediate","edu_high", 'depression', 'clus')]

# check missingness
summary(df_imp)

# set up run
imp0 <- mice(df_imp, maxit = 0, defaultMethod = c('pmm', 'pmm', 'pmm', 'pmm'))

# consider the methods for the variables with missings
meth <- imp0$method
meth
# variables with a imp method -> mp, m3, m9, smoking, ethni, edu

# change predictor matrix
predictormatrix <- imp0$predictorMatrix
# COL - we don't want to use in imputation models
predictormatrix[,c('mother', 'partner','id_dicom', 'depression')] <- 0
# ROW - we don't want to impute: depression
predictormatrix['depression',] <- 0

# visit the secquence
VisSeq <- imp0$visitSequence

# imputations
#implist <- mice(df_imp, m = 30, maxit = 60, seed = 08121996, method = meth, visitSequence = VisSeq, predictorMatrix = predictormatrix)

# save it
#saveRDS(implist, "imp_data_depr_20220322.rds")

# read in the implist - no need to run again
implist <- readRDS("imp_data_depr_20220322.rds")

#----------------------------------------------------------
# FIGURE 1: MEAN TRAJECTORIES DEPRESSION
#----------------------------------------------------------

# means per timepoint per group

## empty dataframe
depr_plot <- as.data.frame(matrix(NA, 12, 2))
colnames(depr_plot) <- c('Time', 'BSI')

## add time and trajectory variable
depr_plot$Time <- rep(c(0, 3, 10), each = 4)
depr_plot$Trajectory <- rep(c('A', 'B', 'C', 'D'))

## scale
df$mp_s <- scale(df$mp)
df$m3_s <- scale(df$m3)
df$m9_s <- scale(df$m9)

## get mean BSI
depr_plot[depr_plot$Time == 0, 'BSI'] <- 
  aggregate(df$mp_s, list(df$clus), FUN=mean, na.rm=T)[,2] # pregnancy
depr_plot[depr_plot$Time == 3, 'BSI'] <- 
  aggregate(df$m3_s, list(df$clus), FUN=mean, na.rm=T)[,2] # at age 3
depr_plot[depr_plot$Time == 10, 'BSI'] <- 
  aggregate(df$m9_s, list(df$clus), FUN=mean, na.rm=T)[,2] # at age 9

## get SD BSI
depr_plot$SD <- NA
depr_plot[depr_plot$Time == 0, 'SD'] <- 
  aggregate(df$mp_s, list(df$clus), FUN=sd, na.rm=T)[,2] # pregnancy
depr_plot[depr_plot$Time == 3, 'SD'] <- 
  aggregate(df$m3_s, list(df$clus), FUN=sd, na.rm=T)[,2] # at age 3
depr_plot[depr_plot$Time == 10, 'SD'] <- 
  aggregate(df$m9_s, list(df$clus), FUN=sd, na.rm=T)[,2] # at age 9

## get N per group
depr_plot$N <- NA
depr_plot[depr_plot$Time == 0, 'N'] <- 
  aggregate(df$mp[complete.cases(df$mp_s)], list(df$clus[complete.cases(df$mp)]), FUN=length)[,2] # pregnancy
depr_plot[depr_plot$Time == 3, 'N'] <- 
  aggregate(df$m3[complete.cases(df$m3_s)], list(df$clus[complete.cases(df$m3)]), FUN=length)[,2] # at age 3
depr_plot[depr_plot$Time == 10, 'N'] <- 
  aggregate(df$m9[complete.cases(df$m9_s)], list(df$clus[complete.cases(df$m9)]), FUN=length)[,2] # at age 9

## get upper and lower limits
depr_plot$LL <- depr_plot$BSI - 1.96 * depr_plot$SD / sqrt(depr_plot$N)
depr_plot$UL <- depr_plot$BSI + 1.96 * depr_plot$SD / sqrt(depr_plot$N)

## write out
write.xlsx(depr_plot, 'Tables/Figure1.xlsx')

## rename
depr_plot$`Depressive symptoms trajectory` <- depr_plot$Trajectory
depr_plot$`Depressive symptoms trajectory`[depr_plot$`Depressive symptoms trajectory` =='A'] <- 'Low symptoms (N = 1,190)'
depr_plot$`Depressive symptoms trajectory`[depr_plot$`Depressive symptoms trajectory` =='B'] <- 'Low increasing symptoms (N = 318)'
depr_plot$`Depressive symptoms trajectory`[depr_plot$`Depressive symptoms trajectory` =='C'] <- 'Decreasing symptoms (N = 128)'
depr_plot$`Depressive symptoms trajectory`[depr_plot$`Depressive symptoms trajectory` =='D'] <- 'High increasing symptoms (N = 40)'

png("Figure1.png", width = 4000, height = 2500,res = 600,type = "cairo-png")

## plot
ggplot(data = depr_plot, aes(x=Time, y=BSI, group=`Depressive symptoms trajectory`)) +
  geom_line(aes(color=`Depressive symptoms trajectory`), size = 0.75)+
  geom_errorbar(aes(ymin=LL, ymax=UL, color=`Depressive symptoms trajectory`), size = 0.75, width = 1, position = position_dodge(0.1)) +
  geom_point(aes(shape=`Depressive symptoms trajectory`, color = `Depressive symptoms trajectory`), size = 2.5,position = position_dodge(0.1))+
  scale_x_continuous(breaks=c(0,3,10))+
  scale_shape_manual(values = c(15, 17, 19, 18))+
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73"))+
  xlab('Time (years)')+
  ylab('Depressive symptoms (z-score)')+
  theme_classic()

dev.off()

# get all data
plotdata_0 <- df[,c("id_dicom",'mp_s', "clus")]
plotdata_3 <- df[,c("id_dicom",'m3_s', "clus")]
plotdata_10 <- df[,c("id_dicom",'m9_s', "clus")]

# create time
plotdata_0$time <- 0
plotdata_3$time <- 3
plotdata_10$time <- 10

# fix names
names(plotdata_0) <- names(plotdata_3) <- names(plotdata_10) <- 
  c("id_dicom", "depression", "clus","time")

# bind
plotdata <- rbind(plotdata_0, plotdata_3, plotdata_10)

## rename
plotdata$`Depressive symptoms trajectory` <- as.character(plotdata$clus)
plotdata$`Depressive symptoms trajectory`[plotdata$`Depressive symptoms trajectory` =='A'] <- 'Low symptoms (N = 1,190)'
plotdata$`Depressive symptoms trajectory`[plotdata$`Depressive symptoms trajectory` =='B'] <- 'Low increasing symptoms (N = 318)'
plotdata$`Depressive symptoms trajectory`[plotdata$`Depressive symptoms trajectory` =='C'] <- 'Decreasing symptoms (N = 128)'
plotdata$`Depressive symptoms trajectory`[plotdata$`Depressive symptoms trajectory` =='D'] <- 'High increasing symptoms (N = 40)'

#------FIGURE S1

library(ggnewscale)

png("FigureS1.png", width = 4500, height = 3000,res = 600,type = "cairo-png")

ggplot(data = depr_plot, aes(x=Time, y=BSI, group=`Depressive symptoms trajectory`))+
  geom_line(data = plotdata, aes(x=time, y=depression, group=id_dicom, color=`Depressive symptoms trajectory`))+
  scale_color_manual(values = c("gray79", "lightgoldenrod1", "lightblue2", "palegreen2"))+
  new_scale_color()+
  geom_line(aes(x=Time, y=BSI, group=`Depressive symptoms trajectory`, color=`Depressive symptoms trajectory`), size = 1)+
  geom_errorbar(aes(ymin=LL, ymax=UL, color=`Depressive symptoms trajectory`), size = 1, width = 1, position = position_dodge(0.1)) +
  geom_point(aes(shape=`Depressive symptoms trajectory`, color = `Depressive symptoms trajectory`), size = 2, position = position_dodge(0.1))+
  scale_x_continuous(breaks=c(0,3,10))+
  scale_shape_manual(values = c(15, 17, 19, 18))+
  scale_color_manual(values = c("black", "goldenrod3", "deepskyblue3", "springgreen4"))+
  xlab('Time (years)')+
  ylab('Depressive symptoms (z-score)')+
  theme_classic()

dev.off()

#----------------------------------------------------------
# FIGURE 2: THE ASSOCIATION BETWEEN TRAJECTORIES AND BRAIN VOLUME
#----------------------------------------------------------

# if * emerge in bse table - check how significant + which predictor

# ADJUSTED FOR AGE, SEX (AND ETIV) 

## global outcomes
glob_1 <- coefficienttable(outcomes = df[,c("totalgrayvol","cerebralwhitemattervol")],
                           pred = df$clus,
                           confounders = df[,c('ageparent', 'type')],
                           neff = 11)

glob_adj_1 <- bse_table(outcomes = df[,c("totalgrayvol","cerebralwhitemattervol")],
                            pred = df$clus,
                            confounders = df[,c('ageparent', 'type')])

## subcortical outcomes
subc_1 <- coefficienttable(outcomes = df[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                           pred = df$clus,
                           confounders = df[,c('ageparent', 'type', 'etiv')],
                           neff = 11)

subc_adj_1 <- bse_table(outcomes = df[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                            pred = df$clus,
                            confounders = df[,c('ageparent', 'type', 'etiv')])

## WML
wml_1 <- coefficienttable(outcomes = df[,"wml_log"],
                          pred = df$clus,
                          confounders = df[,c('ageparent', 'type')],
                          neff = 11)

wml_adj_1 <- bse_table(outcomes = df[,"wml_log"],
                          pred = df$clus,
                          confounders = df[,c('ageparent', 'type')])
names(wml_adj_1) <- 'wml'

## CMB
cmb_1 <- coefficienttable(outcomes = df[,c("cmb")],
                          pred = df$clus,
                          confounders = df[,c('ageparent', 'type')],
                          analysis = 'glm',
                          nef = 11)


cmb_adj_1 <- bse_table(outcomes = df[,c("cmb")],
                          pred = df$clus,
                          confounders = df[,c('ageparent', 'type')],
                          analysis = 'glm')
names(cmb_adj_1) <- 'cmb'


## combine and write out
structural_1 <- rbind(glob_1, subc_1, wml_1, cmb_1)
structural_adj_1a <- cbind(glob_adj_1, wml_adj_1, cmb_adj_1)
structural_adj_1b <- subc_adj_1
write.xlsx(structural_1, 'Tables/Table2_Table3_fullestimates_model1.xlsx') # the regression models

# ADJUSTED FOR AGE, SEX, ETHNICITY, EDUCATION, SMOKING (AND ETIV)

## global outcomes
glob_2 <- coefficienttable(outcomes = df[,c("totalgrayvol","cerebralwhitemattervol")],
                           pred = df$clus,
                           confounders = 'global',
                           neff = 11,
                           MI = T, implist = implist)

glob_adj_2 <- bse_table(outcomes = df[,c("totalgrayvol","cerebralwhitemattervol")],
                            pred = df$clus,
                            confounders = 'global',
                            MI = T, implist = implist)
 
## subcortical outcomes
subc_2 <- coefficienttable(outcomes = df[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                           pred = df$clus,
                           confounders = 'subcortical',
                           neff = 11,
                           MI = T, implist = implist)

subc_adj_2 <- bse_table(outcomes = df[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                            pred = df$clus,
                            confounders = 'subcortical',
                            MI = T, implist = implist)

# WHITE MATTER HYPERINTENSITIES
wml_2 <- coefficienttable(outcomes = df[,"wml_log"],
                          pred = df$clus,
                          confounders = 'global',
                          neff = 11,
                          MI = T, implist = implist)

wml_adj_2 <- bse_table(outcomes = df[,"wml_log"],
                          pred = df$clus,
                          confounders = 'global',
                          MI = T, implist = implist)
names(wml_adj_2) <- 'wml'

# MICROBLEEDS
cmb_2 <- coefficienttable(outcomes = df[,c("cmb")],
                          pred = df$clus,
                          confounders = 'global',
                          analysis = 'glm',
                          neff = 11,
                          MI = T, implist = implist)

cmb_adj_2 <- bse_table(outcomes = df[,c("cmb")],
                          pred = df$clus,
                          confounders = 'global',
                          analysis = 'glm',
                          MI = T, implist = implist)
names(cmb_adj_2) <- 'cmb'

## combine and write out
structural_2 <- rbind(glob_2, subc_2, wml_2, cmb_2)
structural_adj_2a <- cbind(glob_adj_2, wml_adj_2, cmb_adj_2)
structural_adj_2b <- subc_adj_2
write.xlsx(structural_2, 'Tables/Table2_Table3_fullestimates_model2.xlsx') # the regression models

# combine models
structural_adja <- rbind(structural_adj_1a, structural_adj_2a)
structural_adjb <- rbind(structural_adj_1b, structural_adj_2b)
write.xlsx(structural_adja, 'Tables/Table2.xlsx') 
write.xlsx(structural_adjb, 'Tables/Table3.xlsx') 

#----------------------------------------------------------
# FIGURE 3: QDECR (also for sensitivity analyses)
#----------------------------------------------------------

# write out data for qdecr, perform analyses on server

# model 1
df_towrite <- df[,c('id_dicom', 'clus', 'depression', 'ageparent','type')]
saveRDS(df_towrite, "df_qdecr.RData")

# model 1 depression
df_towrite_depr <- df[!is.na(df$depression),c('id_dicom', 'clus', 'depression', 'ageparent','type')]
df_towrite_depr$depression <- as.numeric(scale(df_towrite_depr$depression))
saveRDS(df_towrite_depr, "df_qdecr_depr.RData")

# model 1 women
df_towrite_women <- df[df$type == 2,c('id_dicom', 'clus', 'depression', 'ageparent','type')]
saveRDS(df_towrite_women, "df_qdecr_women.RData")


# model 2
imp2list.mids <- function(x) lapply(seq_len(x$m), function(y) mice::complete(x, y))
maindataimp_clus.complete <- imp2list.mids(implist)
saveRDS(maindataimp_clus.complete, 'maindataimp_clus.complete.rds')

# model 2 depression sens analysis
mids_compl <- complete(implist, action = 'long', include = T)
ids_depr <- unique(mids_compl[mids_compl$.imp==0 & !is.na(mids_compl$depression),'.id'])
mids_compl_depr <- mids_compl[mids_compl$.id %in% ids_depr, ]
mids_compl_depr$depression <- as.numeric(scale(mids_compl_depr$depression))
implist_depr <- as.mids(mids_compl_depr)
maindataimp_depr.complete <- imp2list.mids(implist_depr)
saveRDS(maindataimp_depr.complete, 'maindataimp_depr.complete.rds')
df_depr <- df[!is.na(df$depression),]

# model 2 sensitivity women
mids_compl_w <- mids_compl[which(mids_compl$type==2),]
implist_women <- as.mids(mids_compl_w)
maindataimp_clus_w.complete <- imp2list.mids(implist_women)
saveRDS(maindataimp_clus_w.complete, 'maindataimp_clus.complete_women.rds')

#----------------------------------------------------------
# TABLE S3 AND S4: SEX INTERACTION STRUCTURAL
#----------------------------------------------------------

# ADJUSTED FOR AGE, SEX (AND ETIV) 

# Get implist/df for only women
df_women <- df[df$type == 2, ]
implist_c <- complete(implist, action = 'long', include = T)
implist_women <- as.mids(implist_c[implist_c$type == 2,])

## global outcomes
glob_1_sexint_full <- coefficienttable(outcomes = df_women[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                                       pred = df_women$clus,
                                       confounders = df_women[,c('ageparent')],
                                       neff = 11)

glob_1_sexint <- bse_table(outcomes = df_women[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                           pred = df_women$clus,
                           confounders = df_women[,c('ageparent')],
                           neff = 11)

## subcortical outcomes
subc_1_sexint_full <- coefficienttable(outcomes = df_women[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                                       pred = df_women$clus,
                                       confounders = df_women[,c('ageparent', 'etiv')],
                                       neff = 11)

subc_1_sexint <- bse_table(outcomes = df_women[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                           pred = df_women$clus,
                           confounders = df_women[,c('ageparent', 'etiv')],
                           neff = 11)

# ADJUSTED FOR AGE, SEX, ETHNICITY, EDUCATION, SMOKING (AND ETIV)

## global outcomes
glob_2_sexint_full <- coefficienttable(outcomes = df_women[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                                         pred = df_women$clus,
                                         confounders = 'global',
                                         MI = T, implist = implist_women)

glob_2_sexint <- bse_table(outcomes = df_women[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                           pred = df_women$clus,
                           confounders = 'global',
                           neff = 11,
                           MI = T, implist = implist_women)

## subcortical outcomes
subc_2_sexint_full <- coefficienttable(outcomes = df_women[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                                       pred = df_women$clus,
                                       confounders = 'subcortical',
                                       neff = 11,
                                       MI = T, implist = implist_women)

subc_2_sexint <- bse_table(outcomes = df_women[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                           pred = df_women$clus,
                           confounders = 'subcortical',
                           neff = 11,
                           MI = T, implist = implist_women)

tables1 <- rbind(glob_1_sexint, glob_2_sexint)
tables2 <- rbind(subc_1_sexint, subc_2_sexint)
sexint_full_1 <- rbind(glob_1_sexint_full, subc_1_sexint_full)
sexint_full_2 <- rbind(glob_2_sexint_full, subc_2_sexint_full)
write.xlsx(tables1, 'Tables/TableS5.xlsx')
write.xlsx(tables2, 'Tables/TableS6.xlsx')
write.xlsx(sexint_full_1, 'Tables/TableS5_TableS6_fullestimates_model1.xlsx')
write.xlsx(sexint_full_2, 'Tables/TableS5_TableS6_fullestimates_model2.xlsx')

#----------------------------------------------------------
# TABLE S5: ASSOCIATION CROSS SECTIONAL DEPRESSION WITH VOLUME
#----------------------------------------------------------

## global outcomes
glob_depr_1_full <- coefficienttable_depr(outcomes = df_depr[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                              pred = scale(df_depr$depression),
                              confounders = df_depr[,c('ageparent', 'type')],
                              neff = 11)

glob_depr_1 <- bse_table_depr(outcomes = df_depr[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                                     pred = scale(df_depr$depression),
                                     confounders = df_depr[,c('ageparent', 'type')],
                                     neff = 11)

## cmb
cmb_depr_1_full <- coefficienttable_depr(outcomes = df_depr[,c("cmb")],
                                    pred = scale(df_depr$depression),
                                    confounders = df_depr[,c('ageparent', 'type')],
                                    analysis = 'glm',
                                    neff = 11)

cmb_depr_1 <- bse_table_depr(outcomes = df_depr[,c("cmb")],
                             pred = scale(df_depr$depression),
                             confounders = df_depr[,c('ageparent', 'type')],
                             analysis = 'glm',
                             neff = 11)

## subcortical outcomes
subc_depr_1_full <- coefficienttable_depr(outcomes = df_depr[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                                     pred = scale(df_depr$depression),
                                     confounders = df_depr[,c('ageparent', 'type', 'etiv')],
                                     neff = 11)

subc_depr_1 <- bse_table_depr(outcomes = df_depr[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                              pred = scale(df_depr$depression),
                              confounders = df_depr[,c('ageparent', 'type', 'etiv')],
                              neff = 11)


# ADJUSTED FOR AGE, SEX, ETHNICITY, EDUCATION, SMOKING (AND ETIV)

## global outcomes
glob_depr_2_full <- coefficienttable_depr(outcomes = df_depr[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                                     pred = scale(df_depr$depression),
                                     confounders = 'global',
                                     neff = 11,
                                     MI = T, implist = implist_depr)

glob_depr_2 <- bse_table_depr(outcomes = df_depr[,c("totalgrayvol","cerebralwhitemattervol", "wml_log")],
                              pred = scale(df_depr$depression),
                              confounders = 'global',
                              neff = 11,
                              MI = T, implist = implist_depr)


# cmb
cmb_depr_2_full <- coefficienttable_depr(outcomes = df_depr$cmb,
                                    pred = scale(df_depr$depression),
                                    confounders = 'global',
                                    analysis = 'glm', MI = T, implist = implist_depr,
                                    neff = 11)

cmb_depr_2 <- bse_table_depr(outcomes = df_depr$cmb,
                             pred = scale(df_depr$depression),
                             confounders = 'global',
                             analysis = 'glm', MI = T, implist = implist_depr,
                             neff = 11)

## subcortical outcomes
subc_depr_2_full <- coefficienttable_depr(outcomes = df_depr[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                                     pred = scale(df_depr$depression),
                                     confounders = 'subcortical',
                                     neff = 11,
                                     MI = T, implist = implist_depr)

subc_depr_2 <- bse_table_depr(outcomes = df_depr[,c("amygdala","hippocampus", "thalamus", "caudate", "putamen","pallidum","accumbens")],
                              pred = scale(df_depr$depression),
                              confounders = 'subcortical',
                              neff = 11,
                              MI = T, implist = implist_depr)

## combine and write out
tables4 <- rbind(cbind(glob_depr_1, cmb_depr_1), cbind(glob_depr_2, cmb_depr_2))
tables5 <- rbind(subc_depr_1, subc_depr_2)

depr_full_1 <- rbind(glob_depr_1_full, cmb_depr_1_full, subc_depr_1_full)
depr_full_2 <- rbind(glob_depr_2_full, cmb_depr_2_full, subc_depr_2_full)

write.xlsx(tables4, 'Tables/tableS2.xlsx')
write.xlsx(tables5, 'Tables/tableS3.xlsx')
write.xlsx(depr_full_1, 'Tables/TableS2_TableS3_fullestimates_model1.xlsx')
write.xlsx(depr_full_2, 'Tables/TableS2_TableS3_fullestimates_model2.xlsx')

#----------------------------------------------------------
# AGE DURING GR1003
#----------------------------------------------------------

# read in data
momage <- read_sav("O:/medewerkers/042647 Schuurmans, I/DATA/20211126_invulleeftijd_GR1003.sav")
momagepreg <- read_sav("O:/medewerkers/042647 Schuurmans, I/DATA/MOTHERAGE_BIRTHCHILD_11062017.sav")
general <- read_sav("O:/medewerkers/042647 Schuurmans, I/DATA/CHILD-ALLGENERALDATA_12112020.sav")

# rename
momagepreg$age_birth <- momagepreg$agemother_birthchild
momage$age_quest <- momage$Leeftijdmoeder_invullenGR1003
general$totalpreg <- general$GESTBIR/56

# merge
df1 <- merge(general, momagepreg, by = 'IDC')
df2 <- merge(df1, momage, by = 'IDM')

# retrieve our participants
ids <- df$idparent[df$type==2]
df3 <- df2[df2$MOTHER %in% ids,]

conception <- df3$age_birth-df3$totalpreg
weeks <- (df3$age_quest - conception)*56
weeks <- weeks[which(weeks>0)]

# merge with normal df for post mri
df4 <- merge(df_women, df3, by.x = 'mother', by.y = 'MOTHER')
years <- df4$ageparent - df4$age_birth

#----------------------------------------------------------
# TABLES
#----------------------------------------------------------

# table 1
table1

# figure 1
depr_plot

# table 2
structural_adja

# table 3
structural_adjb

#----------------------------------------------------------
# RESULTS IN TEXT
#----------------------------------------------------------

# ABSTRACT
round(prop.table(table(df$type))*100,1)
round(mean(df$agepre, na.rm=T), 1)
round(sd(df$agepre, na.rm = T), 1)
round(c(mean(years, na.rm = T), sd(years, na.rm = T)), 1)

# METHODS
# flowchart in in merge code

# age 
round(c(mean(weeks, na.rm = T), sd(weeks, na.rm = T)), 1)
round(c(mean(df$age3y/12, na.rm = T), sd(df$age3y/12, na.rm = T)), 1)
round(c(mean(df$age9y, na.rm = T), sd(df$age9y, na.rm = T)), 1)
round(c(mean(years, na.rm = T), sd(years, na.rm = T)), 1)

# age intake and outcome
mean(df$agepre_1, na.rm = T); sd(df$agepre_1, na.rm = T)
mean(df$age3y_1, na.rm = T); sd(df$age3y_1, na.rm = T)
mean(df$age9y_1, na.rm = T); sd(df$age9y_1, na.rm = T)
mean(df$ageparent, na.rm = T); sd(df$ageparent, na.rm = T)

# RESULTS

# study population

## sex
table(df$type)
round(prop.table(table(df$type)),3)*100

# age intake and outcome
mean(df$agepre_1, na.rm = T); sd(df$agepre_1, na.rm = T)
mean(df$age3y_1, na.rm = T); sd(df$age3y_1, na.rm = T)
mean(df$age9y_1, na.rm = T); sd(df$age9y_1, na.rm = T)
mean(df$ageparent, na.rm = T); sd(df$ageparent, na.rm = T)

# education
edu_low <- ifelse(df$edu_intermediate == 0 & df$edu_high == 0, 1, 0)
table(edu_low)
table(df$edu_intermediate)
table(df$edu_high)
round(prop.table(table(df$edu_intermediate)),3)*100

# ethnicity
table(df$ethni)
round(prop.table(table(df$ethni)),3)*100

# marital status
table(df$mardich)
round(prop.table(table(df$mardich)),3)*100

# income
table(df$income_dich)
round(prop.table(table(df$income_dich)),3)*100

# differences between trajectories regarding demographics

## age
summary(aov(df$ageparent ~ df$clus))
aggregate(df$ageparent, list(df$clus), mean)

## sex (2 = women)
chisq.test(df$type, df$clus)
prop.table(table(df$type, df$clus), 2)

## ehtnicity (1 = non-western)
chisq.test(df$ethni, df$clus)
prop.table(table(df$ethni, df$clus), 2)

## education (1 = low education)
chisq.test(df$edu, df$clus)
prop.table(table(df$edu, df$clus), 2)

## marital status
chisq.test(df$mardich, df$clus)
prop.table(table(df$mardich, df$clus), 2)

## income
chisq.test(df$income_dich, df$clus)
prop.table(table(df$income_dich, df$clus), 2)


# 95% CI cortical thickness
ci <- function(M, SE) print(c(round(M-1.96*SE,2), round(M+1.96*SE,2)))
ci(0.0611448094918294,	0.0159585467316554)
ci(0.0624273087950009, 0.0159471239215932)

# GET p adjusted
structural_1$p_corrected <- p.adjust(structural_1$p, method = 'BH')
structural_2$p_corrected <- p.adjust(structural_2$p, method = 'BH')
sexint_full_1$p_corrected <- round(p.adjust(sexint_full_1$p, method = 'BH'),3)
sexint_full_1$p <- round(sexint_full_1$p, 3)
sexint_full_2$p_corrected <- round(p.adjust(sexint_full_2$p, method = 'BH'),3)
sexint_full_2$p <- round(sexint_full_2$p, 3)
depr_full_1$p_corrected <- round(p.adjust(depr_full_1$p, method = 'BH'),3)
depr_full_1$p <- round(depr_full_1$p, 3)
depr_full_2$p_corrected <- round(p.adjust(depr_full_2$p, method = 'BH'),3)
depr_full_2$p <- round(depr_full_2$p, 3)

# DISCUSSION

# infarcts
round(prop.table(table(df$li)),3)

# SUPPLEMENTARY

# 95% CI cortical thickness/area depression
ci(-0.0214771218835927,	0.00568039369689057)
ci(-0.00739146033023943, 0.0021484110625148)
ci(-0.0122337214910914,	0.00327698495293642)


# 95% CI cortical thickness/area women
ci(0.0564659502822906,	0.0149064753825466)
ci(0.0740257910321761,	0.0188167581065583)
ci(0.0745965205194032,	0.0187750111166082)

## add medication use
meds <-  read_sav("O:/medewerkers/042647 Schuurmans, I/DATA/MEDICATIONSELFREPORTPREGNANCY_30112017.sav")
meds <- merge(general[,c('IDM', 'MOTHER')], meds, by = 'IDM')
meds$medication <- 0
meds$medication[meds$SSRITOT < 5] <- 1
meds$medication[which(meds$Q01OTHR8 == 1)] <- 1
meds$medication[which(meds$Q02OTHR8 == 1)] <- 1
meds$medication[which(meds$Q03OTHR8 == 1)] <- 1
meds$medication[which(meds$Q01TCA == 1)] <- 1
meds$medication[which(meds$Q02TCA == 1)] <- 1
meds$medication[which(meds$Q03TCA == 1)] <- 1
meds_rm <- meds[meds$medication == 1,]
meds_rm <- meds_rm[!duplicated(meds_rm$MOTHER),]
df_meds <- merge(meds_rm[,c('MOTHER', 'medication')], df_women, by.x = 'MOTHER', by.y = 'mother', all.y = T)

# stats
(table(df_meds$medication)); 19/1117# 1.5% of the women used SSRI during pregnancy - no MAO users
table(df_meds$medication, df_meds$clus) # majority was in the high increasing group (18.1% used SSRI)
table(df_meds$clus)
