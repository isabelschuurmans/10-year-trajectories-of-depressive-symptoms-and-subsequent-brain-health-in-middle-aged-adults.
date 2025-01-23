# NON RESPONSE ANALYSES OUDER MRI
# Isabel Schuurmans, 2021-11-15

# open library
library(haven)
library(xlsx)

# function
nonresponse <- function(ids_study, ids_full, dataset, variables){
  
  # convert dataset to dataframe
  dataset <- as.data.frame(dataset)[]
  
  # make indicator variable for included and excluded participants
  dataset$indicator <- ids_full %in% ids_study
  
  # get empty table
  table_stats <- matrix(NA, nrow = (length(variables)), ncol = 3)
  
  # compare on the variables
  
  for (i in 1:length(variables)){
    
    # t-test continuous outcomes
    if (is.numeric(dataset[,variables[i]])){
      
      table_stats[i,1] <- as.character(paste0(round(mean(dataset[dataset$indicator == T,variables[i]], na.rm=T),2), ' ± ', round(sd(dataset[dataset$indicator == T,variables[i]], na.rm=T), 2), ')', sep = ''))
      table_stats[i,2] <- as.character(paste0(round(mean(dataset[dataset$indicator == F,variables[i]], na.rm=T),2), ' ± ', round(sd(dataset[dataset$indicator == F,variables[i]], na.rm=T), 2), ')', sep = ''))
      table_stats[i,3] <- as.character(paste0('t (df) = ', as.numeric(round(t.test(dataset[dataset$indicator == T,variables[i]], dataset[dataset$indicator == F,variables[i]])$statistic, 2)), ' (',
                                              round(t.test(dataset[dataset$indicator == T,variables[i]], dataset[dataset$indicator == F,variables[i]])$parameter, 1), '), p = ',
                                              round(t.test(dataset[dataset$indicator == T,variables[i]], dataset[dataset$indicator == F,variables[i]])$p.value, 3)))
      
    } 
    
    # chi-square for categorical outcomes 
    else {
      
      table_stats[i,1] <- as.character(paste0(round(as.numeric(table(dataset[dataset$indicator == T,variables[i]])[1])),' (', 
                                              round(as.numeric(prop.table(table(dataset[dataset$indicator == T,variables[i]]))[1])*100,1), ')'))
      table_stats[i,2] <- as.character(paste0(round(as.numeric(table(dataset[dataset$indicator == F,variables[i]])[1])),' (', 
                                              round(as.numeric(prop.table(table(dataset[dataset$indicator == F,variables[i]]))[1])*100,1), ')'))
      table_stats[i,3] <- as.character(paste0('X2 (df) = ', round(chisq.test(dataset[,variables[i]], dataset$indicator)$statistic, 2), ' (',
                                              round(chisq.test(dataset[,variables[i]], dataset$indicator)$parameter, 1), '), p = ',
                                              round(chisq.test(dataset[,variables[i]], dataset$indicator)$p.value, 3)))
      
    }
    
  }
  
  # clean up stats
  colnames(table_stats) <- c('M study', 'M full', 'test stats')
  rownames(table_stats) <- variables
  
  return(table_stats)
  
}

# read in data
tmp <- readRDS("O:/medewerkers/042647 Schuurmans, I/Project_1_ELS_CM/tmp.rds")
gen <- read_sav("O:/medewerkers/042647 Schuurmans, I/DATA/CHILD-ALLGENERALDATA_12112020.sav")
mri <- readRDS("O:/Medewerkers/042647 Schuurmans, I/Project_3_DEP_MRI/data_depr_project.rds")
tmp <- merge(tmp, gen, by.x = 'idc', by.y = 'IDC')

# recode smoking 
tmp$smoking <- ifelse(tmp$smoking ==3, 2, tmp$smoking)

# recode ethnicity
tmp$ethnicity_women <- tmp$ETHNMv2
tmp$ethnicity_women[tmp$ETHNMv2 == 1] <- 0
tmp$ethnicity_women[tmp$ETHNMv2 == 300] <- 0
tmp$ethnicity_women[tmp$ETHNMv2 == 500] <- 0
tmp$ethnicity_women[tmp$ETHNMv2 == 700] <- 0
tmp$ethnicity_women[tmp$ETHNMv2 == 800] <- 0
tmp$ethnicity_women[tmp$ETHNMv2 == 2] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 3] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 4] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 5] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 6] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 7] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 200] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 400] <- 1
tmp$ethnicity_women[tmp$ETHNMv2 == 600] <- 1
tmp$ethnicity_women <- as.numeric(tmp$ethnicity_women)

tmp$ethnicity_men <- tmp$ETHNFv2
tmp$ethnicity_men[tmp$ETHNFv2 == 1] <- 0
tmp$ethnicity_men[tmp$ETHNFv2 == 300] <- 0
tmp$ethnicity_men[tmp$ETHNFv2 == 500] <- 0
tmp$ethnicity_men[tmp$ETHNFv2 == 700] <- 0
tmp$ethnicity_men[tmp$ETHNFv2 == 800] <- 0
tmp$ethnicity_men[tmp$ETHNFv2 == 2] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 3] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 4] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 5] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 6] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 7] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 200] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 400] <- 1
tmp$ethnicity_men[tmp$ETHNFv2 == 600] <- 1
tmp$ethnicity_men <- as.numeric(tmp$ethnicity_men)

# recode  edu
tmp$education_women3 <- ifelse(tmp$EDUCM3 < 3, 0, 1)
tmp$education_men3 <- ifelse(tmp$EDUCP3 < 3, 0, 1)

# recode mari
tmp$mari <- as.numeric(tmp$mari)

# age men
tmp$age_inclusion_men <- ifelse(!is.na(tmp$PARTBF) & !is.na(tmp$AGE_BF_v2), tmp$AGE_BF_v2, 
                                ifelse(!is.na(tmp$PARTBF) & is.na(tmp$AGE_BF_v2), tmp$AGE_P_v2, 
                                       ifelse(is.na(tmp$AGE_BF_v2), tmp$AGE_P_v2, tmp$AGE_BF_v2)))


# make van tmp a parent data
tmp1a <- tmp[,c('MOTHER', 'gender', 'gestbir', 'weight','parity','age_m_v2', 'ethnicity_women','BMI_0','education_women3', 'mari','dep','dbsi3m', 'visit13','smoking')]
tmp1a$parentID <- paste0(tmp1a$MOTHER, 'm')
tmp1b <- tmp[,c('PARTNER', 'gender', 'gestbir', 'weight','parity','age_inclusion_men', 'ethnicity_men','BMI_P','education_men3', 'mari','dep_p','dbsi3f', 'visit13')]
tmp1b$smoking <- rep(NA,nrow(tmp1b))
tmp1b$parentID <- paste0(tmp1b$PARTNER, 'P')

# bind
names(tmp1b) <- names(tmp1a) <- c('ID', 'gender', 'gestbir', 'weight','parity','age_inclusion', 'ethnicity','BMI','education', 'marital_status','dep_pregn','dep_3', 'visit13','smoking','parentID')
tmp2 <- rbind(tmp1a, tmp1b)

# change ids
mri$parentID <- rep(NA, nrow(mri))
for (i in 1:nrow(mri)) {
  mri$parentID[i] <- ifelse(mri$type[i]==2, paste0(mri$idparent[i], 'm', sep = ''), paste0(mri$idparent[i], 'p', sep = ''))}

# convert to correct data types GENDER
tmp2[,c('gender','ethnicity', 'education', 'marital_status', 'smoking')] <- 
  apply(tmp2[,c('gender','ethnicity', 'education', 'marital_status', 'smoking')], 2, as.factor)

tmp2$marital_status[tmp2$marital_status == "NaN"] <- NA

# get ids that participated in @13
ids_13 <- tmp2[which(tmp2$visit13==1),]

# run non-response
tab <- nonresponse(ids_study = mri$parentID, ids_full = ids_13$parentID, dataset = ids_13, 
                   variables = c('age_inclusion', 'ethnicity','BMI','education', 'marital_status','dep_pregn','dep_3', 'gender', 'gestbir', 'weight','parity'))

write.xlsx(tab, 'O:/Medewerkers/042647 Schuurmans, I/Project_3_DEP_MRI/Tables/TableS1.xlsx')



