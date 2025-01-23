# functions depression project
# function to get stats
get_table1 <- function(df){
  
  table <- data.frame(rep(NA, ncol(df)+1)) 
  
  table[1,1] <- paste0(nrow(df), ' (', round(nrow(df)/1676*100,1), ')', collapse = '')
  
  for (i in 1:ncol(df)){
    
    # get n & % for factor data
    if (is.factor(df[,i]) | is.character(df[,i])){
      
      n <- length(df[which(df[,i]==1),i])
      pt <- round(n/nrow(df)*100,1)
      table[(i+1),1] <- paste0(n, ' (', pt, ')', collapse = '')
      
    }
    
    else {
      
      m <- round(mean(df[,i], na.rm=T),2)
      sd <- round(sd(df[,i], na.rm=T),2)
      table[(i+1),1] <- paste0(m, ' ± ', sd, collapse = '')
      
    }
    
  }
  
  return(table)
  
}

# change from two to one column  
ci_interval <- function(table){
  
  for (i in 1:nrow(table)){
    table[i,2] <- (paste(table[i,2], table[i,3], sep = ', '))
  }
  
  table <- table[,c(1:2,4)]
  return(table)
  
}

# coefficient table
coefficienttable <- function(outcomes, pred, confounders, analysis = 'lm', neff, interactor, MI = FALSE, implist){
  
  # make an empty table 
  outcomes <- as.data.frame(outcomes)
  confounders <- as.data.frame(confounders)
  table_1 <- data.frame(matrix(NA, length(outcomes)*4, 4))
  
  if (analysis == 'lm'){
    
    if (MI){
      
      for (i in 1:ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking)),
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv)))
        
        # save in table
        position <- i*4-2
        table_1[position:(position+2),] <- data.frame(round(summary(pool(mod))[2:4,2],2), 
                                                      round(summary(pool(mod))[2:4,2]-1.96*summary(pool(mod))[2:4,3],2),
                                                      round(summary(pool(mod))[2:4,2]+1.96*summary(pool(mod))[2:4,3],2),
                                                      summary(pool(mod))[2:4,6])
        
      }  
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        x <- lm(outcomes[,i] ~ ., data = totalpred) %>% summary()
        
        # save in table
        position <- i*4-2
        table_1[position:(position+2),] <- data.frame(round(x$coefficients[2:4],2),
                                                      round(x$coefficients[2:4]-1.96*x$coefficients[2:4,2],2),
                                                      round(x$coefficients[2:4]+1.96*x$coefficients[2:4,2],2),
                                                      x$coefficients[2:4,4])
      }
    }
  }  
  
  # add row and column names
  names(table_1) <- c('beta', 'CI', 'CI ul', 'p')
  names_new <- rep(NA,4*length(outcomes))
  abcd <- c('A', 'B', 'C', 'D')
  
  for (k in 1:length(outcomes)){
    
    for (l in 1:length(abcd)){
      names_new[(k*4 - 4 + l)] <- paste0(names(outcomes[k]), sep = '_', abcd[l])
    }
    
  }
  
  if (analysis == 'glm'){
    
    if (MI){
      
      for (i in ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking, family = 'binomial')),
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv, family = 'binomial')))
        
        # save in table
        position <- i*4-2
        table_1[position:(position+2),] <- data.frame(round(exp(summary(pool(mod))[2:4,2]),2), 
                                                      round(exp(summary(pool(mod))[2:4,2]-1.96*summary(pool(mod))[2:4,3]),2),
                                                      round(exp(summary(pool(mod))[2:4,2]+1.96*summary(pool(mod))[2:4,3]),2),
                                                      summary(pool(mod))[2:4,6])
      }  
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        x <- glm(outcomes[,i] ~ ., data = totalpred, family = 'binomial') %>% summary()
        
        # save in table
        position <- i*4-2
        table_1[position:(position+2),] <- data.frame(round(exp(x$coefficients[2:4]),2),
                                                      round(exp(x$coefficients[2:4]-1.96*x$coefficients[2:4,2]),2),
                                                      round(exp(x$coefficients[2:4]+1.96*x$coefficients[2:4,2]),2),
                                                      x$coefficients[2:4,4])
      }
    }
  }  
  
  # add row and column names
  names(table_1) <- c('Odd', 'CI', 'CI ul', 'p')
  names_new <- rep(NA,4*length(outcomes))
  abcd <- c('A', 'B', 'C', 'D')
  
  for (k in 1:length(outcomes)){
    
    for (l in 1:length(abcd)){
      names_new[(k*4 - 4 + l)] <- paste0(names(outcomes[k]), sep = '_', abcd[l])
    }
    
  }
  
  row.names(table_1) <- names_new
  
  # add confidence interval
  table_2 <- ci_interval(table_1)
  
  # round p
  table_2$p <- round(table_2$p, 3)

  return(table_2)
  
}

# new coefficient table
bse_table <- function(outcomes, pred, confounders, analysis = 'lm', neff, interactor, MI = FALSE, implist){
  
  # make an empty table 
  outcomes <- as.data.frame(outcomes)
  confounders <- as.data.frame(confounders)
  table_1 <- data.frame(matrix(NA, 3, length(outcomes)))
  
  if (analysis == 'lm'){
    
    if (MI){
      
      for (i in 1:ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking)),
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv)))
        
        # save in table
        x<-as.data.frame(ifelse(summary(pool(mod))[2,6] < .05 | summary(pool(mod))[3,6] < .05 | summary(pool(mod))[4,6] < .05,
                                                  data.frame(paste0(round(summary(pool(mod))[2:4,2],2), '* (', 
                                                                    round(summary(pool(mod))[2:4,2] - 1.96*summary(pool(mod))[2:4,3],2), ', ',
                                                                    round(summary(pool(mod))[2:4,2] + 1.96*summary(pool(mod))[2:4,3],2),')')),
                                                  data.frame(paste0(round(summary(pool(mod))[2:4,2],2), ' (', 
                                                                    round(summary(pool(mod))[2:4,2] - 1.96*summary(pool(mod))[2:4,3],2), ', ',
                                                                    round(summary(pool(mod))[2:4,2] + 1.96*summary(pool(mod))[2:4,3],2),')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1:3,i] <- x$obs
          
      }  
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        mod <- lm(outcomes[,i] ~ ., data = totalpred) %>% summary()
        
        # save in table
        x<-as.data.frame(ifelse(mod$coefficients[2,4] < .05 | mod$coefficients[3,4] < .05 | mod$coefficients[4,4] < .05,
                                data.frame(paste0(round(mod$coefficients[2:4],2), '* (', 
                                                  round(mod$coefficients[2:4]-1.96*mod$coefficients[2:4,2],2), ', ',
                                                  round(mod$coefficients[2:4]+1.96*mod$coefficients[2:4,2],2), ')')),
                                data.frame(paste0(round(mod$coefficients[2:4],2), ' (', 
                                                  round(mod$coefficients[2:4]-1.96*mod$coefficients[2:4,2],2), ', ',
                                                  round(mod$coefficients[2:4]+1.96*mod$coefficients[2:4,2],2), ')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1:3,i] <- x$obs
        
        
      }
    }
  }  
  
  if (analysis == 'glm'){
    
    if (MI){
      
      for (i in ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking, family = 'binomial')),
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv, family = 'binomial')))
        
        # save in table
        x<-as.data.frame(ifelse(summary(pool(mod))[2,6] < .05 | summary(pool(mod))[3,6] < .05 | summary(pool(mod))[4,6] < .05,
                                data.frame(paste0(round(exp(summary(pool(mod))[2:4,2]),2), '* (', 
                                                  round(exp(summary(pool(mod))[2:4,2] - 1.96*summary(pool(mod))[2:4,3]),2), ', ',
                                                  round(exp(summary(pool(mod))[2:4,2] + 1.96*summary(pool(mod))[2:4,3]),2), ')')),
                                data.frame(paste0(round(exp(summary(pool(mod))[2:4,2]),2), ' (',
                                                  round(exp(summary(pool(mod))[2:4,2] - 1.96*summary(pool(mod))[2:4,3]),2), ', ',
                                                  round(exp(summary(pool(mod))[2:4,2] + 1.96*summary(pool(mod))[2:4,3]),2), ')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1:3,i] <- x$obs
        
      }
      
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        mod <- glm(outcomes[,i] ~ ., data = totalpred, family = 'binomial') %>% summary()
        
        # save in table
        x<-as.data.frame(ifelse(mod$coefficients[2,4] < .05 | mod$coefficients[3,4] < .05 | mod$coefficients[4,4] < .05,
                                data.frame(paste0(round(exp(mod$coefficients[2:4]),2), '* (', 
                                                  round(exp(mod$coefficients[2:4] - 1.96*mod$coefficients[2:4,2]),2), ', ',
                                                  round(exp(mod$coefficients[2:4] + 1.96*mod$coefficients[2:4,2]),2), ')')),
                                data.frame(paste0(round(exp(mod$coefficients[2:4]),2), ' (', 
                                                  round(exp(mod$coefficients[2:4] - 1.96*mod$coefficients[2:4,2]),2), ', ',
                                                  round(exp(mod$coefficients[2:4] + 1.96*mod$coefficients[2:4,2]),2), ')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1:3,i] <- x$obs
      }
    }
  }  
  
  
  # add row and column names
  names(table_1) <- names(outcomes)
  rownames(table_1) <- c('low increasing', 'decreasing', 'high increasing')
  
  return(table_1)
  
}



bse_table_depr <- function(outcomes, pred, confounders, analysis = 'lm', neff, interactor, MI = FALSE, implist){
  
  # make an empty table 
  outcomes <- as.data.frame(outcomes)
  table_1 <- data.frame(matrix(NA, 1, length(outcomes)))
  
  if (analysis == 'lm'){
    
    if (MI){
      
      for (i in 1:ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking)),
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv)))
        
        # save in table
        x<-as.data.frame(ifelse(summary(pool(mod))[2,6] < .05,
                                data.frame(paste0(round(summary(pool(mod))[2,2],2), '* (', 
                                                  round(summary(pool(mod))[2,2] - 1.96*summary(pool(mod))[2,3],2), ', ',
                                                  round(summary(pool(mod))[2,2] + 1.96*summary(pool(mod))[2,3],2),')')),
                                data.frame(paste0(round(summary(pool(mod))[2,2],2), ' (', 
                                                  round(summary(pool(mod))[2,2] - 1.96*summary(pool(mod))[2,3],2), ', ',
                                                  round(summary(pool(mod))[2,2] + 1.96*summary(pool(mod))[2,3],2),')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1,i] <- x$obs
        
      }  
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        mod <- lm(outcomes[,i] ~ ., data = totalpred) %>% summary()
        
        # save in table
        x<-as.data.frame(ifelse(mod$coefficients[2,4] < .05,
                                data.frame(paste0(round(mod$coefficients[2],2), '* (', 
                                                  round(mod$coefficients[2]-1.96*mod$coefficients[2,2],2), ', ',
                                                  round(mod$coefficients[2]+1.96*mod$coefficients[2,2],2), ')')),
                                data.frame(paste0(round(mod$coefficients[2],2), ' (', 
                                                  round(mod$coefficients[2]-1.96*mod$coefficients[2,2],2), ', ',
                                                  round(mod$coefficients[2]+1.96*mod$coefficients[2,2],2), ')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1,i] <- x$obs
        
        
      }
    }
  }  
  
  if (analysis == 'glm'){
    
    if (MI){
      
      for (i in ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking, family = 'binomial')),
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv, family = 'binomial')))
        
        # save in table
        x<-as.data.frame(ifelse(summary(pool(mod))[2,6] < .05,
                                data.frame(paste0(round(exp(summary(pool(mod))[2,2]),2), '* (', 
                                                  round(exp(summary(pool(mod))[2,2] - 1.96*summary(pool(mod))[2,3]),2), ', ',
                                                  round(exp(summary(pool(mod))[2,2] + 1.96*summary(pool(mod))[2,3]),2), ')')),
                                data.frame(paste0(round(exp(summary(pool(mod))[2,2]),2), ' (',
                                                  round(exp(summary(pool(mod))[2,2] - 1.96*summary(pool(mod))[2,3]),2), ', ',
                                                  round(exp(summary(pool(mod))[2,2] + 1.96*summary(pool(mod))[2,3]),2), ')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1,i] <- x$obs
        
      }
      
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        mod <- glm(outcomes[,i] ~ ., data = totalpred, family = 'binomial') %>% summary()
        
        # save in table
        x<-as.data.frame(ifelse(mod$coefficients[2,4] < .05,
                                data.frame(paste0(round(exp(mod$coefficients[2]),2), '* (', 
                                                  round(exp(mod$coefficients[2] - 1.96*mod$coefficients[2,2]),2), ', ',
                                                  round(exp(mod$coefficients[2] + 1.96*mod$coefficients[2,2]),2), ')')),
                                data.frame(paste0(round(exp(mod$coefficients[2]),2), ' (', 
                                                  round(exp(mod$coefficients[2] - 1.96*mod$coefficients[2,2]),2), ', ',
                                                  round(exp(mod$coefficients[2] + 1.96*mod$coefficients[2,2]),2), ')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1,i] <- x$obs
      }
    }
  }  
  
  
  # add row and column names
  names(table_1) <- names(outcomes)
  rownames(table_1) <- c('depression')
  
  return(table_1)
  
}


bse_table_sexint <- function(outcomes, pred, confounders, analysis = 'lm', neff, moderator, MI = FALSE, implist){
  
  # make an empty table 
  outcomes <- as.data.frame(outcomes)
  table_1 <- data.frame(matrix(NA, 7, length(outcomes)))
  tot <- ncol(as.data.frame(confounders))
  
  if (analysis == 'lm'){
    
    if (MI){
      
      for (i in 1:ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, lm(outcomes[,i] ~ pred * moderator + ageparent + ethni + edu_intermediate + edu_high + smoking)),
               mod <- with(implist, lm(outcomes[,i] ~ pred * moderator + ageparent + ethni + edu_intermediate + edu_high + smoking + etiv)))
        
        # save in table
        x<-as.data.frame(ifelse(summary(pool(mod))[tot+6,6] < .05 | summary(pool(mod))[tot+7,6] < .05 | summary(pool(mod))[tot+8,6] < .05,
                                data.frame(paste0(round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2],2), '* (', 
                                                  round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2] - 1.96*summary(pool(mod))[c(2:5, (tot+6):(tot+8)),3],2), ', ',
                                                  round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2] + 1.96*summary(pool(mod))[c(2:5, (tot+6):(tot+8)),3],2),')')),
                                data.frame(paste0(round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2],2), ' (', 
                                                  round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2] - 1.96*summary(pool(mod))[c(2:5, (tot+6):(tot+8)),3],2), ', ',
                                                  round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2] + 1.96*summary(pool(mod))[c(2:5, (tot+6):(tot+8)),3],2),')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1:7,i] <- x$obs
        
      }  
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        ifelse(tot == 1, 
               mod <- summary(lm(outcomes[,i] ~ pred * moderator + confounders)),
               mod <- summary(lm(outcomes[,i] ~ pred * moderator + confounders[,1] + confounders[,2])))
        
        # save in table
        x<-as.data.frame(ifelse(mod$coefficients[tot+6,4] < .05 | mod$coefficients[tot+7,4] < .05 | mod$coefficients[tot+8,4] < .05,
                                data.frame(paste0(round(mod$coefficients[c(2:5, (tot+6):(tot+8))],2), '* (', 
                                                  round(mod$coefficients[c(2:5, (tot+6):(tot+8))]-1.96*mod$coefficients[c(2:5, (tot+6):(tot+8)),2],2), ', ',
                                                  round(mod$coefficients[c(2:5, (tot+6):(tot+8))]+1.96*mod$coefficients[c(2:5, (tot+6):(tot+8)),2],2), ')')),
                                data.frame(paste0(round(mod$coefficients[c(2:5, (tot+6):(tot+8))],2), ' (', 
                                                  round(mod$coefficients[c(2:5, (tot+6):(tot+8))]-1.96*mod$coefficients[c(2:5, (tot+6):(tot+8)),2],2), ', ',
                                                  round(mod$coefficients[c(2:5, (tot+6):(tot+8))]+1.96*mod$coefficients[c(2:5, (tot+6):(tot+8)),2],2), ')')))[[1]])
        x$obs <- as.character(x[,1])
        table_1[1:7,i] <- x$obs
        
        
      }
    }
  }  
  

  
  # add row and column names
  names(table_1) <- names(outcomes)
  rownames(table_1) <- c('low increasing', 'decreasing', 'high increasing', 'sex', 'sex x low increasing', 'sex x decreasing', 'sex x high increasing')
  
  return(table_1)
  
}




# get adjusted mean differences, test-statistic & p-value






adjustedmeans <- function(outcomes, pred, confounders, analysis = 'lm', neff, interactor, MI = FALSE, implist){
  
  library(emmeans)
  
  if (analysis == 'lm'){
    
    # make an empty table 
    outcomes <- as.data.frame(outcomes)
    table_1 <- data.frame(matrix(NA, length(outcomes)*4, 4))
    
    if (MI){
      
      for (i in 1:ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               x <- with(implist, lm(outcomes[,i] ~ clus + ageparent + type + ethni + edu_intermediate + edu_high + smoking)),
               x <- with(implist, lm(outcomes[,i] ~ clus + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv)))
        
        # get emmeans
        mm <- as.data.frame(emmeans(x, specs = pairwise ~ clus)$emmeans)
        mm_p <- as.data.frame(emmeans(x, specs = pairwise ~ clus)$contrasts)
        
        # save in table 
        position <- i*4-3
        table_1[(position):(position+3),] <- data.frame(round(mm[,2],2),
                                                        round(mm[,5],2),
                                                        round(mm[,6],2),
                                                        round(c(NA, mm_p[1:3,6]),3))
        
        }
      }  
    
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        x <- lm(outcomes[,i] ~ ., data = totalpred) 
        
        # get emmeans
        mm <- as.data.frame(emmeans(x, specs = pairwise ~ pred)$emmeans)
        mm_p <- as.data.frame(emmeans(x, specs = pairwise ~ pred)$contrasts)
        
        # save in table 
        position <- i*4-3
        table_1[(position):(position+3),] <- data.frame(round(mm[,2],2),
                                                        round(mm[,5],2),
                                                        round(mm[,6],2),
                                                        round(c(NA, mm_p[1:3,6]),3))
      }
    }
  }  
  
  # add row and column names
  names(table_1) <- c('Mean', 'LL', 'UL', 'p')
  names_new <- rep(NA,4*length(outcomes))
  abcd <- c('A', 'B', 'C', 'D')
  
  for (k in 1:length(outcomes)){
    
    for (l in 1:length(abcd)){
      names_new[(k*4 - 4 + l)] <- paste0(names(outcomes[k]), sep = '_', abcd[l])
    }
    
  }
  
  row.names(table_1) <- names_new
  
  table_1$outcome <- rep(names(outcomes), each = 4)
  table_1$group <- rep(abcd, length(outcomes))
  
  table_2 <- table_1
  
  return(table_2)
  
}


coefficienttable_sexint <- function(outcomes, pred, moderator, confounders, neff, MI = FALSE, implist){
  
  # make an empty table 
  outcomes <- as.data.frame(outcomes)
  table_1 <- data.frame(matrix(NA, length(outcomes)*7, 4))
  tot <- ncol(as.data.frame(confounders))
  
  if (MI){
    
    for (i in 1:ncol(outcomes)){
      
      ifelse(confounders == 'global', 
             mod <- with(implist, lm(outcomes[,i] ~ pred * moderator + ageparent + ethni + edu_intermediate + edu_high + smoking)),
             mod <- with(implist, lm(outcomes[,i] ~ pred * moderator + ageparent + ethni + edu_intermediate + edu_high + smoking + etiv)))
      
      
      
      # save in table
      position <- i*7-6
      table_1[position:(position+6),] <- data.frame(round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2],2), 
                                                    round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2]-1.96*summary(pool(mod))[c(2:5, (tot+6):(tot+8)),3],2),
                                                    round(summary(pool(mod))[c(2:5, (tot+6):(tot+8)),2]+1.96*summary(pool(mod))[c(2:5, (tot+6):(tot+8)),3],2),
                                                    summary(pool(mod))[c(2:5, (tot+6):(tot+8)),6])
      
    }  
  }
  
  else { 
    
    # regessions
    for (i in 1:ncol(outcomes)){
      
      # make predictor sections
      ifelse(tot == 1, 
             x <- summary(lm(outcomes[,i] ~ pred * moderator + confounders)),
             x <- summary(lm(outcomes[,i] ~ pred * moderator + confounders[,1] + confounders[,2])))
      
      
      # save in table
      position <- i*7-6
      table_1[position:(position+6),] <- data.frame(round(x$coefficients[c(2:5, (tot+6):(tot+8))],2),
                                                    round(x$coefficients[c(2:5, (tot+6):(tot+8))]-1.96*x$coefficients[c(2:5, (tot+6):(tot+8)),2],2),
                                                    round(x$coefficients[c(2:5, (tot+6):(tot+8))]+1.96*x$coefficients[c(2:5, (tot+6):(tot+8)),2],2),
                                                    x$coefficients[c(2:5, (tot+6):(tot+8)),4])
    }
  }
  
  
  # add row and column names
  names(table_1) <- c('beta', 'CI', 'CI ul', 'p')
  names_new <- rep(NA,7*length(outcomes))
  names_int <- c('B', 'C', 'D', 'type', 'sex_B', 'sex_C', 'sex_D')
  
  for (k in 1:length(outcomes)){
    
    for (l in 1:length(names_int)){
      names_new[(k*7 - 7 + l)] <- paste0(names(outcomes[k]), sep = '_', names_int[l])
    }}
  
  
  row.names(table_1) <- names_new
  
  # add confidence interval
  table_2 <- ci_interval(table_1)
  
  # round p
  #table_2$p <- round(table_2$p, 3)

  return(table_2)
  
}

coefficienttable_depr <- function(outcomes, pred, confounders, analysis = 'lm', neff, interactor, MI = FALSE, implist){
  
  if (analysis == 'lm'){
    
    # make an empty table 
    outcomes <- as.data.frame(outcomes)
    table_1 <- data.frame(matrix(NA, length(outcomes), 4))
    
    if (MI){
      
      for (i in 1:ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking)),
               mod <- with(implist, lm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv)))
        
        # save in table
        table_1[i,] <- data.frame(round(summary(pool(mod))[2,2],2), 
                                  round(summary(pool(mod))[2,2]-1.96*summary(pool(mod))[2,3],2),
                                  round(summary(pool(mod))[2,2]+1.96*summary(pool(mod))[2,3],2),
                                  summary(pool(mod))[2,6])
        
      }  
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        x <- lm(outcomes[,i] ~ ., data = totalpred) %>% summary()
        
        # save in table
        table_1[i,] <- data.frame(round(x$coefficients[2],2),
                                  round(x$coefficients[2]-1.96*x$coefficients[2,2],2),
                                  round(x$coefficients[2]+1.96*x$coefficients[2,2],2),
                                  x$coefficients[2,4])
      }
    }
  }  
  
  if (analysis == 'glm'){
    
    # make an empty table 
    outcomes <- as.data.frame(outcomes)
    table_1 <- data.frame(matrix(NA, length(outcomes), 4))
    
    if (MI){
      
      for (i in 1:ncol(outcomes)){
        
        # run regression
        ifelse(confounders == 'global', 
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking, family = 'binomial')),
               mod <- with(implist, glm(outcomes[,i] ~ pred + ageparent + type + ethni + edu_intermediate + edu_high + smoking + etiv, family = 'binomial')))
        
        # save in table
        table_1[i,] <- data.frame(round(summary(pool(mod))[2,2],2), 
                                  round(summary(pool(mod))[2,2]-1.96*summary(pool(mod))[2,3],2),
                                  round(summary(pool(mod))[2,2]+1.96*summary(pool(mod))[2,3],2),
                                  summary(pool(mod))[2,6])
      }  
    }
    
    else { 
      
      # regessions
      for (i in 1:ncol(outcomes)){
        
        # make predictor sections
        totalpred <- cbind(pred, confounders)
        
        # run regression
        x <- glm(outcomes[,i] ~ ., data = totalpred, family = 'binomial') %>% summary()
        
        # save in table
        table_1[i,] <- data.frame(round(exp(x$coefficients[2]),2),
                                  round(exp(x$coefficients[2]-1.96*x$coefficients[2,2]),2),
                                  round(exp(x$coefficients[2]+1.96*x$coefficients[2,2]),2),
                                  x$coefficients[2,4])
      }
    }
  }  
  
  # add row and column names
  names(table_1) <- c('Odd/Beta', 'CI', 'CI ul', 'p')
  row.names(table_1) <- names(outcomes)
  
  # add confidence interval
  table_2 <- ci_interval(table_1)
  
  return(table_2)
  
}

