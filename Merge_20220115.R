######################################
# DIFFERENTIATING THE MATERNAL BRAIN #
######################################

# MERGE SCRIPT 

# Authors: I.K. Schuurmans, R. Zou, H. El Marroun, ..., A. Luik. 
# 2020-03-03

# NOTE: files mothers and father are kept seperately untill analysis as their ID numbers have overlap!

#----------------------------------------------------------
# LIBRARY
#----------------------------------------------------------

# Open library
library(foreign)
library(magrittr)
library(dplyr)
library(xlsx)

setwd("O:/Medewerkers/042647 Schuurmans, I/Project_3_DEP_MRI")

#----------------------------------------------------------
# FUNCTIONS
#----------------------------------------------------------

### Read in data quickly
readquick <- function(path){ # only works for spss
  dataframe <- read.spss(path, use.value.labels = FALSE, to.data.frame = T)
  names(dataframe) <- tolower(names(dataframe))
  return(dataframe)
}

#----------------------------------------------------------
# BASE: GENERAL INFO
#----------------------------------------------------------

# GENERAL
general <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/CHILD-ALLGENERALDATA_12112020.sav")

## ethnicity

### father
general$ethni_f <- rep(NA, nrow(general))
for (i in 1:nrow(general)){
  if (is.na(general$ethnfv2[i])){
    general$ethni_f[i] <- NA}
  else{
    if (general$ethnfv2[i] == 2|general$ethnfv2[i] == 3|general$ethnfv2[i] == 4|
        general$ethnfv2[i] == 5|general$ethnfv2[i] == 6|general$ethnfv2[i] == 7|
        general$ethnfv2[i] == 200|general$ethnfv2[i] == 400|general$ethnfv2[i] == 600){
      general$ethni_f[i] <- 1}
    else {
      if (general$ethnfv2[i] == 1|general$ethnfv2[i] == 300| general$ethnfv2[i] == 500|
          general$ethnfv2[i] == 700|general$ethnfv2[i] == 800){
        general$ethni_f[i] <- 0}}}}

general$ethni_p <- rep(NA, nrow(general))
for (i in 1:nrow(general)){
  if (is.na(general$ethnp[i])){
    general$ethni_p[i] <- NA}
  else{
    if (general$ethnp[i] == 2|general$ethnp[i] == 3|general$ethnp[i] == 4|
        general$ethnp[i] == 5|general$ethnp[i] == 6|general$ethnp[i] == 7|
        general$ethnp[i] == 200|general$ethnp[i] == 400|general$ethnp[i] == 600){
      general$ethni_p[i] <- 1}
    else {
      if (general$ethnp[i] == 1|general$ethnp[i] == 300| general$ethnp[i] == 500|
          general$ethnp[i] == 700|general$ethnp[i] == 800){
        general$ethni_p[i] <- 0}}}}

general$ethni_f <- ifelse(general$partbf == 1, general$ethni_f, general$ethni_p)

### mother
general$ethni <- rep(NA, nrow(general))
for (i in 1:nrow(general)){
  if (is.na(general$ethnmv2[i])){
    general$ethni[i] <- NA}
  else{
    if (general$ethnmv2[i] == 2|general$ethnmv2[i] == 3|general$ethnmv2[i] == 4|
        general$ethnmv2[i] == 5|general$ethnmv2[i] == 6|general$ethnmv2[i] == 7|
        general$ethnmv2[i] == 200|general$ethnmv2[i] == 400|general$ethnmv2[i] == 600){
      general$ethni[i] <- 1}
    else {
      if (general$ethnmv2[i] == 1|general$ethnmv2[i] == 300| general$ethnmv2[i] == 500|
          general$ethnmv2[i] == 700|general$ethnmv2[i] == 800){
        general$ethni[i] <- 0}}}}

## education
general$edu_intermediate <- ifelse(general$educm == 3 | general$educm == 4, 1, 0)
general$edu_intermediate_partner <- ifelse(general$educp == 3 | general$educp == 4, 1, 0)
general$edu_high <- ifelse(general$educm > 4, 1, 0)
general$edu_high_partner <- ifelse(general$educp > 4, 1, 0)

## income
general$income_dich <- ifelse(general$income < 7, 1, 0)

# get confounders
conf <- general[,c('idc', 'mother', 'partner', 'ethni_f', 'ethni', 'edu_intermediate', 'edu_intermediate_partner', 'edu_high', 'edu_high_partner')]

# fix doubles (as some mothers/partners are repeated because they have multiple children)
conf_m <- conf %>% # dataframe 
  filter(!is.na(conf$mother)) %>% # only take the not missing and usable moms
  group_by(mother) %>% # group by moms
  filter(edu_high == max(edu_high)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(mother)) 

conf_p <- conf %>% # dataframe 
  filter(!is.na(conf$partner)) %>% # only take the not missing and usable moms
  group_by(partner) %>% # group by moms
  filter(edu_high_partner == max(edu_high_partner)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(partner)) 

#----------------------------------------------------------
# OUTCOME: BRAIN
#----------------------------------------------------------

## read in core MRI data
mri_core <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/Ouder MRI/MRIparent_MRIcore_data_02072021.sav")
levels(mri_core$id_dicom) <- trimws(levels(mri_core$id_dicom))
mri_core$id_dicom <- as.character(mri_core$id_dicom)

## read in structural data (t1)
freesurfer1 <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/Ouder MRI/ParentMRI_aseg_stats_20210601.sav")
freesurfer2 <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/Ouder MRI/ParentMRI_tbv_20210601.sav")
freesurfer <- merge(freesurfer1, freesurfer2, by = c('id'), all = T)
levels(freesurfer$id) <- trimws(levels(freesurfer$id))
freesurfer$id <- as.character(freesurfer$id)

# get measures
freesurfer$amygdala <- as.numeric(scale((freesurfer$left_amygdala_vol + freesurfer$right_amygdala_vol)/2))
freesurfer$hippocampus <- as.numeric(scale((freesurfer$left_hippocampus_vol + freesurfer$right_hippocampus_vol)/2))
freesurfer$thalamus <- as.numeric(scale((freesurfer$left_thalamus_proper_vol + freesurfer$right_thalamus_proper_vol)/2))
freesurfer$caudate <- as.numeric(scale((freesurfer$left_caudate_vol + freesurfer$right_caudate_vol)/2))
freesurfer$putamen <- as.numeric(scale((freesurfer$left_putamen_vol + freesurfer$right_putamen_vol)/2))
freesurfer$pallidum <- as.numeric(scale((freesurfer$left_pallidum_vol + freesurfer$right_pallidum_vol)/2))
freesurfer$accumbens <- as.numeric(scale((freesurfer$left_accumbens_area_vol + freesurfer$right_accumbens_area_vol)/2))
freesurfer$totalgrayvol <- as.numeric(scale(freesurfer$totalgrayvol))
freesurfer$cerebralwhitemattervol <- as.numeric(scale(freesurfer$cerebralwhitemattervol))

# select outcome measures
freesurfer <- freesurfer[,c('id','totalgrayvol', 'cerebralwhitemattervol', "amygdala", 
                            "hippocampus", "thalamus", "caudate", "putamen","pallidum", 
                            "accumbens", 'etiv', 'wm_hypointensities_vol')]

## read CVD markers (flair and t2* infarcts + cmb)
cvd <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/Ouder MRI/MRIparent_CMB_infarcten_31052021.sav")

# merge mri data
df1 <- merge(mri_core, freesurfer, by.x = 'id_dicom', by.y = 'id', all = T)
df3 <- merge(df1, cvd[,c(1:2, 5:10)], by = c('mother', 'partner'), all = T)

#----------------------------------------------------------

# PREDICTOR: LONGITUDINAL DEPRESSION

## read in data
bsipre <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/prenatal/GR1003-BSI D1_22112016.sav")
bsipref <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/depression/GR1004-BSI G1_22112016.sav")
agepre <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/20211126_invulleeftijd_GR1003.sav")
source("O:/Medewerkers/042647 Schuurmans, I/Project_1_ELS_CM/R Code/SyntaxBSI_20191230.R")
names(BSI_totalscore) <- tolower(names(BSI_totalscore))
bsinine <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/GR1081_E1-3_5-6_08082016.sav")

## merge
da <- merge(general[,c("idc", "idm", "mother","partner", "age_m_v2","age_bf_v2","age_p_v2", 'mardich', 'income_dich')], bsipre, all = T, by = 'idm')
db <- merge(da, bsipref, all = T, by = 'idm')
dc <- merge(db, agepre, all = T, by  = 'idm')
dd <- merge(dc, BSI_totalscore, all = T, by = 'idc')
df <- merge(dd, bsinine, all = T, by = 'idc')

## rename depression (mother)
df$mp <- df$dep
df$m3 <- df$dbsi3m
df$m9 <- df$dbsi9m

## rename depression (father)
df$fp <- df$dep_p
df$f3 <- df$dbsi3f
df$f9 <- df$dbsi9f

## number of moments filled in
df$mmis <- rowSums(is.na(data.frame(df$mp, df$m3, df$m9))) 
df$fmis <- rowSums(is.na(data.frame(df$fp, df$f3, df$f9))) 

## age mom or dad intake
df$mage <- df$age_m_v2
df$fage <- ifelse(is.na(df$age_bf_v2), df$age_p_v2, df$age_bf_v2)

## rename age 
df$agepre <- df$leeftijdmoeder_invullengr1003  
df$agepre_dad <- df$fage + (df$leeftijdmoeder_invullengr1003-df$age_m_v2)
df$age3y <- df$bsi.age_gr1065
df$age3y_dad <- df$bsi.age_gr1065
df$age9y <- df$bsi.agemothergr1081
df$age9y_dad <- df$bsi.agefathergr1083

## extract only the necessary variable, make different dataset for man and woman
dm <- df[,c('mother', 'mp', 'm3', 'm9', 'mmis', 'agepre', 'age3y', 'age9y', 'mage', 'mardich', 'income_dich')]
dp <- df[,c('partner', 'fp', 'f3', 'f9','fmis', 'agepre_dad', 'age3y_dad', 'age9y_dad', 'fage', 'mardich', 'income_dich')]

# fix doubles (as some mothers/partners are repeated because they have multiple children)
dm <- dm %>% # dataframe 
  filter(!is.na(dm$mother)) %>% # only take the not missing and usable moms
  group_by(mother) %>% # group by moms
  filter(mmis == min(mmis)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(mother)) 

dp <- dp %>% # dataframe 
  filter(!is.na(dp$partner)) %>% # only take the not missing and usable moms
  group_by(partner) %>% # group by moms
  filter(fmis == min(fmis)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(partner)) 

#----------------------------------------------------------

# VISIT DATA DEPRESSION

# read in the data
last <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/Ouder MRI/MRIparent_interview_03052021.sav")
last$depression <- last$cesdscorewgt

# divide for mother and father
last_mother <- last[!is.na(last$mother),c('mother', 'depression')]
last_partner <- last[!is.na(last$partner),c('partner', 'depression')]

# remove doubles
last_mother <- last_mother %>% # dataframe 
  filter(!is.na(last_mother$mother)) %>% # only take the not missing and usable moms
  group_by(mother) %>% # group by moms
  filter(depression == max(depression)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(mother)) 

last_partner <- last_partner %>% # dataframe 
  filter(!is.na(last_partner$partner)) %>% # only take the not missing and usable moms
  group_by(partner) %>% # group by moms
  filter(depression == max(depression)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(partner)) 

# merge with earlier depression data
dm <- merge(dm, last_mother, by = 'mother', all = T)
dp <- merge(dp, last_partner, by = 'partner', all.x = T)

#----------------------------------------------------------

# CONFOUNDER: SMOKING

# smoking mother
smoking_mother <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/GR1001-F1-10_22112016.sav")
smoking_mother <- merge(general[,c(1:3)], smoking_mother, by = 'idm', all.y = T)
smoking_mother$smoking <- as.numeric(smoking_mother$f0700101)
#smoking_mother$smoking[which(smoking_mother$f0700401 == 0)] <- 2
#smoking_mother$smoking <- as.factor(smoking_mother$smoking)

# smoking partner
smoking_partner <- readquick("O:/medewerkers/042647 Schuurmans, I/DATA/GR1001-G1-7_22112016.sav")
smoking_partner <- merge(general[,c(1:2, 4)], smoking_partner, by = 'idm', all.y = T)
smoking_partner$smoking <- ifelse(smoking_partner$g0100101 == 2, NA, smoking_partner$g0100101)

# fix doubles (as some mothers/partners are repeated because they have multiple children)
smoking_mother <- smoking_mother %>% # dataframe 
  filter(!is.na(smoking_mother$mother)) %>% # only take the not missing and usable moms
  group_by(mother) %>% # group by moms
  filter(smoking == max(smoking)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(mother)) 

smoking_partner <- smoking_partner %>% # dataframe 
  filter(!is.na(smoking_partner$partner)) %>% # only take the not missing and usable moms
  group_by(partner) %>% # group by moms
  filter(smoking == max(smoking)) %>% # take by duplicated scans the one with the higher quality
  filter(!duplicated(partner)) 

#----------------------------------------------------------
# FLOWCHART

# merge depression, brain and smoking
dfmother <- merge(df3[df3$type==2,], dm, all.x = T, by = 'mother')
dfmother <- merge(dfmother, smoking_mother[,c('mother', 'smoking')], all.x = T, by = 'mother')
dfmother <- merge(dfmother, conf_m[,c('mother', 'ethni', 'edu_intermediate', 'edu_high')], all.x = T, by = 'mother')

dfpartner <- merge(df3[df3$type==3,], dp, all.x = T, by = 'partner')
dfpartner <- merge(dfpartner, smoking_partner[,c('partner', 'smoking')], all.x = T, by = 'partner')
dfpartner <- merge(dfpartner, conf_p[,c('partner', 'ethni_f', 'edu_intermediate_partner', 'edu_high_partner')], all.x = T, by = 'partner')

names(dfpartner) <- names(dfmother)
df <- rbind(dfmother, dfpartner)

# block 1
block1 <- nrow(df)

# block 2: missings depression
df2 <- df[which(df$mmis<2),]
block2 <- nrow(df2)

# ADD THIS STEP TO MS
# block 3: missing neuroimaging (t1, wml, flair, t2*)
df3a <- df2[which(df2$t1_available==1),]
df3 <- df3a[which(df3a$t2_available==1),]
block3 <- nrow(df3)

# block 4: ifs
df4 <- df3[which(df3$exclude_incidental == 0),]
block4 <- nrow(df4)

# block 5: quality (t1, wml, flair, t2*)
df5a <- df4[which(df4$exclude_quality_t1==0),]
df5 <- df5a[which(df5a$exclude_quality_t2==0),]
block5 <- nrow(df5)

# flow
c(block1, block2-block1, block2, block3-block2, block3, block5-block4, block4, block4-block3, block5)

# final dataset
df <- df5
saveRDS(df, 'data_depr_project.rds')

#----------------------------------------------------------
# check ethnicities
table(general[which(df$idparent[df$type==2] %in% general$mother), "ethnmv2"])
table(general[which(df$idparent[df$type==3] %in% general$mother), "ethnfv2"]) #CHEXCK
# non-western: Indonesian (2), Cape verdian (3), Morrocan (4), Dutch Antilles (5), 
# Surinamese (6), Turkish (7), African (200), Asian, non-western (400), African, non-western (600)
# western: Dutch (1), American, western (300), Asian, western (500), 
# European (700), Oceanie (800)
  
# check mean ages
ages <- df[c("idparent","type", "agepre", "age3y", "age9y")]
names(ages) <- c("idparent","type", "agepre_1", "age3y_1", "age9y_1")
mean(ages$agepre_1, na.rm = T); sd(ages$agepre_1, na.rm = T)
mean(ages$age3y_1, na.rm = T); sd(ages$age3y_1, na.rm = T)
mean(ages$age9y_1, na.rm = T); sd(ages$age9y_1, na.rm = T)

# merge with df because we wanted all ages but do now want to rerun depression trajs
summary(ages)
df <- readRDS("df.rds")
df_1 <- merge(df, ages, by = c('idparent', 'type'), all = T)
saveRDS(df_1, 'df_1.rds')
