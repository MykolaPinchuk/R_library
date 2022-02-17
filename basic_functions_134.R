### This is a collection of functions to work with CRSP and Compustat. Here are all useful functions I produced since August 2018. ###

### Versions ###
# in version 26 I replaced ME with me.
# between v.77 and v.80 I deleted bunch of obsolete functions.
# there were bf88 and bf88t. I forgot the difference and continued with 88t. This is 88t successor. similarly 97__ vs 97.
# version 93 uses the latest (06/30/19) fama_mac and deletes its previous versions.
# in v_117 I deleted some old breakpoints and unisorters.
# in v_122 I deleted some obsolete functions.
# around Feb-March 2021 I have created v129, but smth went wrong and it was lost forever.

# for simple grouped sum pracma::accumarray vastly outperforms data.table for small vectors. on accumarray vs aggregate():
# aggregate has higher fixed cost and lower marginal cost. if vector in below 500, use accumarray. if above 1000, use aggregate.
# is above 5000, dt is comparable to aggregate.

# on prd: prd variable depends on the market i want to study. futures market will have different prd than equity market.
# therefore, i should try to save all input files without prd variables, with "date" only.
# then, depending on the type of analysis i need, i should load respective prd_map and join it with data.

## WARNING: dplyr:: group_by %>% mutate()   is not what it seems. need to use group_by %>% do() to replicate for loop. 
## some code here uses group_by %>% mutate() and may be wrong.
## all modern functions hopefully use data.table.

## WARNING: beware of floating point issues. two largest numbers in q=seq(-0.75, -0.15, 0.15) are not -0.15 and -0.3.
## they are very close, but not equal to -0.15 and -0.3. while q==-0.6 suceeds, q==-0.3 fails.

## whenever data.table join gives the non-uniqueness warning, using allow.cartesian=TRUE is a mistake!
## it means I have some duplicates, so I have to solve that first.

## to report nice t-stats in famamac, need: 
## https://stackoverflow.com/questions/42995868/r-stargazer-package-eliminate-t-label-from-reported-test-statistics
## trace(stargazer:::.stargazer.wrap, edit = T)

# roll_regress can spit out the error smth like "dchdd failed with a code 1...". It is not debuggable. The most frequent reason is constant regressor, often it is RET=0 or RET=-RF. 
# can delete such permno-ymprd via "permno_ymprd_zeroRET15d.csv".
# good piece on dt joins: https://medium.com/analytics-vidhya/r-data-table-joins-48f00b46ce29

# Coding tips:
# Use Ctrl+Alt+B to run evth before the cursor.
# use ctrl+Shift+C to comment out a chunk of code.
# use Ctrl+Alt+T to run code section-by-section with sections, divided by ####

# Meta-programming tips:
# I tend to start implementing smth too early w/o thinking enough about the plan.
# at the margin, think a bit longer on how you plan implementing smth before you start coding.
# in the end, it will save more time and you will be less likely to be stuck with ugly legacy code too cumbersome too rewrite.
# use sections to organise the code, they are useful.
# try to write code within each section to be reusable, i.e., not dependent ob objects, modified in place.
# to do that, you can fix the objects at the end of section (e.g., crsp0 <- copy(crsp)) and then use fixed objects at the beginning of the next section.




library("data.table")
library("lubridate")
library("dplyr")
library("Rcpp")
library("RcppArmadillo")
#options("scipen"=100, "digits"=4)


### Column SD ###
colSd <- function(x, nao=T)sqrt(rowMeans((t(x)-colMeans(x, na.rm=nao))^2, na.rm=nao)*((dim(x)[1])/(dim(x)[1]-1)))
#END


### Lagger ###
lag <- function(x,l=1){
  output <- c(rep(NA, l), x[1:(length(x)-l)])
  return (output)
}
#END

### Matrix Lagger ###
m_lag <- function(x, l=1, na=T){
  if(na==T){
    tempdf <- data.frame(matrix(NA, l, ncol(x)))
    colnames(tempdf) <- colnames(x)
    output <- rbind(tempdf, x[1:(nrow(x)-l),])
  }
  if(na==F){
    tempdf <- data.frame(matrix(0, l, ncol(x)))
    colnames(tempdf) <- colnames(x)
    output <- rbind(tempdf, x[1:(nrow(x)-l),])
  }
  colnames(output) <- paste0("l", l, colnames(output))
  return (output)
}
#END


############################
### FAST AUTOCORRELATION ###
############################
# added around early May 2020

acr <- function(x, l){
  sum((x[(-l):-1]-mean(x))*(x[-length(x):-(length(x)+l-1)]-mean(x)))/sum((x-mean(x))^2)
}

# END



#####################
### CPCRSP LAGGER ###
#####################
# traditional way: just lag the whole variable, do couple of checks
# dplyr way: group by permno, enumerate by prd, fill the holes with left_join, then just use "lag()" to the vector of interest
# dplyr way is more elegant, but is probably slower. template for dplyr way is in "...59t_cpcrsp_lagger"

cpcrsp_lagger <- function(data, period, variable, lag){
  
  ##cpst <- cpcrsp_lagger(cpst, "year", "at", 1)
  
  data <- data[order(data$PERMNO, data[,period]),]
  
  temp <- m_lag(data[,c("PERMNO", period, variable)], lag)
  cnames <- colnames(temp)
  output <- cbind(data, temp)
  output <- output[(1+lag):nrow(output),]

  output[!((output[, "PERMNO"] == output[, cnames[1]]) & (output[, period] == output[, cnames[2]]+lag)), cnames[3]] <- NA
  output[, cnames[1:2]] <- list(NULL)
  
  return(output)
  
}

#END





### Fast %in% ###
library("fastmatch")
`%fin%` <- function(x, table) {
  stopifnot(require(fastmatch))
  fmatch(x, table, nomatch = 0L) > 0L
}
#END



### Summary statistic ###
# added on 05/11/20
# corrected on 01/06/21

summstat <- function(x, digits = 4){
  
  output <- data.frame(matrix(-88, ncol(x), 11))
  output <- data.frame(t(round(rbind(apply(x, 2, min, na.rm = TRUE), apply(x, 2, quantile, probs=c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99), na.rm = TRUE), apply(x, 2, max, na.rm=TRUE), apply(x, 2, mean, na.rm = TRUE), apply(x, 2, sd, na.rm = TRUE)), digits)))
  setnames(output, names(output), c("Min", "p1", "p10", "p25", "Median", "p75", "p90", "p99", "Max", "Mean", "SD"))
  
  return(output)
  
}

# END




### P-values for correlations ###
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z <- round(z, 5)
  z
}
#END



####################
### ROLLING MEAN ###
####################


rolling_mean <- function(data, variables, periods){
  
  ##ret_t_ <- rolling_mean(ret_t, c("eRETtp35", "mkt"), 24)
  
  data <- data[order(data$year, data$month),]
  data[, paste(variables, "rm", sep="")] <- rep(NA, nrow(data))
  
  if(length(variables)==1){
    for(i in 1:(nrow(data)-periods)){
      data[periods+i, paste(variables, "rm", sep="")] <- mean(data[i:(i+periods-1),variables])
    }
  }
  else{
    for(i in 1:(nrow(data)-periods)){
      data[periods+i, paste(variables, "rm", sep="")] <- colMeans(data[i:(i+periods-1),variables])
    }
  }
  return(data)
  
}

#END


#########################
### ROLLING PREDICTOR ###
#########################
# predictor using ols via rolling\expanding window
# it should be general, ie it should use 1-dimensional period variable (original or reduced)
# added around March 2019

roll_predict <- function(data, variable, ols_predictors, skipprd, type = "rolling", rollprd){
  
  ##Eg: dpdata_roll30 <- roll_predict(dpdata, "epr", "ldp", 10, type = "rolling", rollprd = 30)
  
  data[, paste(variable, "_pred_mean", sep="")] <- rep(NA, nrow(data))
  data[, paste(variable, "_pred_ols", sep="")] <- rep(NA, nrow(data))
  
  if(type == "expanding"){
    
    for (i in (skipprd+1):nrow(data)){
      data[, paste(variable, "_pred_mean", sep="")][i] <- mean(data[, variable][1:(i-1)], na.rm=T)
      data[, paste(variable, "_pred_ols", sep="")][i] <- (lm.fit(cbind(rep(1,(i-1)),data[, ols_predictors][1:(i-1)]), data[, variable][1:(i-1)])$coefficients)%*%c(1,data[, ols_predictors][i])
    }
    
  }
  
  if(type == "rolling"){
    
    for (i in (skipprd+1):nrow(data)){
      data[, paste(variable, "_pred_mean", sep="")][i] <- mean(data[, variable][(max(1, i-rollprd)):(i-1)], na.rm=T)
      data[, paste(variable, "_pred_ols", sep="")][i] <- (lm.fit(cbind(rep(1,(i - (max(1, i-rollprd)))),data[, ols_predictors][(max(1, i-rollprd)):(i-1)]), data[, variable][(max(1, i-rollprd)):(i-1)])$coefficients)%*%c(1,data[, ols_predictors][i])
    }
    
  }
  
  
  return(data)
  
}

#END



#######################################################
### CRSP CLEANER FROM PERMNOS WITH FEW OBSERVATIONS ###
#######################################################
# general function, applicable for either monthly or daily returns
# added on 11/18, v.42
# as far as I remember, it is needed to guarantee correct operation of make_momentum

clean_crsp_by_obs_old <- function(data, N){
  
  #crspc <- clean_crsp_by_obs(crspc, 45), i.e., kill permnos with less than 45 days of returns
  
  crsp_summ <- data %>%
    group_by(PERMNO) %>%
    dplyr::summarize(n = n()) %>%
    as.data.frame()
  
  permnos_todie <- crsp_summ$PERMNO[crsp_summ$n < N]
  
  crsp <- data[!(data$PERMNO %fin% permnos_todie),]
  
  return(crsp)
}

# END



# added on 12/14/20 to replace dplyr syntax with dt.

clean_crsp_by_obs <- function(data, n){
  
  #crspc <- clean_crsp_by_obs(crspc, 45), i.e., kill permnos with less than 45 days of returns
  
  crsp_summ <- data[,.N,by=PERMNO]
  crsp <- data[PERMNO %in% crsp_summ[N>=n,PERMNO]]

  return(crsp)
}

# END





############################################################################
### THE FUNCTION TO PRODUCE NICE OUTPUT AFTER average_returns_calculator ###
############################################################################
# It will create asterisks and square brackets
# usually, there n ntile sorts and long-short portfolio, so nportfolios = # of sorts + 1
# created on 11/22/18
# col_skip added on 11/23/18

make_stars <- function(data, nportfolios, col_skip = 2){
  
  ast_th <- qnorm(c(0.95, 0.975, 0.995))
  temptsts <- data[seq(2, nrow(data), 2), (col_skip+1):ncol(data)]
  output <- data
  output <- output[seq(1, nrow(output)-1, 2), (col_skip+1):ncol(output)]
  
  for (i in seq(1, nrow(output), 1)){
    output[i,] <- sprintf("%0.2f", output[i,])
  }
  
  for (j in 1:(nportfolios)){
    output[((abs(temptsts[,j])< 100500)&(abs(temptsts[,j])>ast_th[3])),j] <- paste0(output[((abs(temptsts[,j])< 100500)&(abs(temptsts[,j])>ast_th[3])),j], "as3")
    output[((abs(temptsts[,j])< ast_th[3])&(abs(temptsts[,j])>ast_th[2])),j] <- paste0(output[((abs(temptsts[,j])< ast_th[3])&(abs(temptsts[,j])>ast_th[2])),j], "as2")
    output[((abs(temptsts[,j])< ast_th[2])&(abs(temptsts[,j])>ast_th[1])),j] <- paste0(output[((abs(temptsts[,j])< ast_th[2])&(abs(temptsts[,j])>ast_th[1])),j], "as1")
  }
  
  for (i in seq(1, nrow(temptsts), 1)){
    temptsts[i,] <- sprintf("%0.2f", temptsts[i,])
    temptsts[i,] <- paste0("[", temptsts[i,], "]")
  }
  
  output$index <- seq(1, (nrow(output)*2-1), 2)
  temptsts$index <- seq(2, (nrow(temptsts)*2), 2)
  
  results <- rbind(output, temptsts)
  results <- results[order(results$index),]
  results$index <- NULL
  results <- cbind(data[,1:col_skip], results)
  
  return(results)
  
}

# END





############################################################################
### THE FUNCTION TO PRODUCE NICE OUTPUT AFTER average_returns_calculator ###
############################################################################
# It will create asterisks and square brackets
# usually, there n ntile sorts and long-short portfolio, so nportfolios = # of sorts + 1
# created on 10/14/19

make_stars_19 <- function(data, col_skip = 0, brackets_exist = FALSE){
  
  ## make_stars_19(copy(table4), 1, TRUE)
  ## make_stars_19(copy(table3))
  
  ast_th <- qnorm(c(0.95, 0.975, 0.995))
  temptsts <- data[seq(2, nrow(data), 2), (col_skip+1):ncol(data)]
  output <- data
  output <- output[seq(1, nrow(output)-1, 2), (col_skip+1):ncol(output)]
  nportfolios <- ncol(data)-col_skip
  
  if(col_skip>0){
    ornames <- names(data)[1:col_skip]
  }
  
  if(brackets_exist == TRUE){
    temptsts <- gsub("\\[", "", temptsts)
    temptsts <- gsub("\\]", "", temptsts)
    temptsts <- data.frame(t(as.numeric(temptsts)))
  }
  
  for (i in seq(1, nrow(output), 1)){
    output[i,] <- sprintf("%0.2f", output[i,])
  }

  for (j in 1:(nportfolios)){
    output[((abs(temptsts[,j])< 100500)&(abs(temptsts[,j])>ast_th[3])),j] <- paste0(output[((abs(temptsts[,j])< 100500)&(abs(temptsts[,j])>ast_th[3])),j], "as3")
    output[((abs(temptsts[,j])< ast_th[3])&(abs(temptsts[,j])>ast_th[2])),j] <- paste0(output[((abs(temptsts[,j])< ast_th[3])&(abs(temptsts[,j])>ast_th[2])),j], "as2")
    output[((abs(temptsts[,j])< ast_th[2])&(abs(temptsts[,j])>ast_th[1])),j] <- paste0(output[((abs(temptsts[,j])< ast_th[2])&(abs(temptsts[,j])>ast_th[1])),j], "as1")
  }
  
  for (i in seq(1, nrow(temptsts), 1)){
    temptsts[i,] <- sprintf("%0.2f", temptsts[i,])
    temptsts[i,] <- paste0("[", temptsts[i,], "]")
  }
  
  output$index <- seq(1, (nrow(output)*2-1), 2)
  temptsts$index <- seq(2, (nrow(temptsts)*2), 2)
  names(temptsts) <- names(output)
  
  results <- rbind(output, temptsts)
  results <- results[order(results$index),]
  results$index <- NULL
  if (col_skip>0){
    results <- cbind(data[,1:col_skip], results)
    names(results)[1:col_skip] <- ornames
  }
  
  return(results)
  
}

# END


##################################################
### The function to add stars in summary stats ###
##################################################

# added on 3/7/2020 for Yixin`s RA.
# data can be data.frame with row.names. 


make_stars_yc <- function(data){
  
  ## table1ew <- make_stars_yc(table1ew)
  
  ast_th <- qnorm(c(0.95, 0.975, 0.995))
  
  rnames <- row.names(data)
  setDT(data)
  
  data[,LS:=sprintf("%0.3f", LS)]
  data[] <- lapply(data, as.numeric)
  data[,stars:=(abs(Tstat)>=ast_th[3])+(abs(Tstat)>=ast_th[2])+(abs(Tstat)>=ast_th[1])]
  data[,LS:=paste0(LS, paste0("as", stars))]
  data <- data.frame(data[,stars:=NULL])
  rownames(data) <- rnames
  
  return(data)
  
}

#END




##########################################################
### THE FUNCTION TO PRODUCE NICE OUTPUT AFTER BJS_regs ###
##########################################################
# It will create asterisks and square brackets
# It is intended to be used on BJS_regs output, i.e. 4x_ df with returns and t-stats below. 
# created on 09/03/19

make_stars_BJS <- function(data, nportf = 2){
  
  ##make_stars_BJS(table3)
  
  output1 <- data[1:2,-1]
  output2 <- data.frame()
  ast_th <- qnorm(c(0.95, 0.975, 0.995))
  temp <- gsub("]", "", output1[2,])
  temp <- gsub("\\[ ", "", temp)
  temp <- as.numeric(temp)
  
  for (j in 1:(ncol(output1))){
    output1[((abs(temp[j])< 100500)&(abs(temp[j])>ast_th[3])),j][1] <- paste0(output1[((abs(temp[j])< 100500)&(abs(temp[j])>ast_th[3])),j][1], "as3")
    output1[((abs(temp[j])< ast_th[3])&(abs(temp[j])>ast_th[2])),j][1] <- paste0(output1[((abs(temp[j])< ast_th[3])&(abs(temp[j])>ast_th[2])),j][1], "as2")
    output1[((abs(temp[j])< ast_th[2])&(abs(temp[j])>ast_th[1])),j][1] <- paste0(output1[((abs(temp[j])< ast_th[2])&(abs(temp[j])>ast_th[1])),j][1], "as1")
  }
  
  if(nportf==2){
    
    output2 <- data[3:4,-1]
    temp <- gsub("]", "", output2[2,])
    temp <- gsub("\\[ ", "", temp)
    temp <- as.numeric(temp)
    
    for (j in 1:(ncol(output2))){
      output2[((abs(temp[j])< 100500)&(abs(temp[j])>ast_th[3])),j][1] <- paste0(output2[((abs(temp[j])< 100500)&(abs(temp[j])>ast_th[3])),j][1], "as3")
      output2[((abs(temp[j])< ast_th[3])&(abs(temp[j])>ast_th[2])),j][1] <- paste0(output2[((abs(temp[j])< ast_th[3])&(abs(temp[j])>ast_th[2])),j][1], "as2")
      output2[((abs(temp[j])< ast_th[2])&(abs(temp[j])>ast_th[1])),j][1] <- paste0(output2[((abs(temp[j])< ast_th[2])&(abs(temp[j])>ast_th[1])),j][1], "as1")
    }
  }
  
  output <- cbind(data[,1], rbind(output1, output2))
  names(output) <- names(data)
  
  return(output)
  
}

#END


###############################
### DEPENDENT DOUBLE SORTER ###
###############################
# [On time-invariant variable]
# probably obsolete

# REQUIRES: NOTHING
# ARGUMENTS:
# data - data.frame
# nportf - scalar, number of portfolios to sort into
# variable1 - a variable, on which the data is already sorted. is of the form "decilesmth"
# variable2 - a variable, on which we want to sort conditionally on sorting on variable1
# OUTPUT: data_sorted, a dataframe with added column, indicating the portfolio the firm belongs to


dependent_double_sorter_inv <- function(data, variable1, variable2, nportf){
  
  ##crsp_sb <- dependent_double_sorter(crspn, "decileME", "beta", 10)
  
  intoutput <- matrix(NA, nrow <- nrow(data)+1, ncol = ncol(data)+1)
  colnames(intoutput) <- c(colnames(data), paste("decilecond", variable2, sep = ""))
  
  i <- 1
  for (yr in unique(data$year)){
    if (sum(((data$year == yr)&(data$month == 6))) == 0){
      next
    }
    tempdata <- data[(data$year == yr), ]
    for (d in 1:10){
      ttempdata <- tempdata[tempdata[, variable1] == d,]
      ttempdata[, paste("decilecond", variable2, sep = "")] <- ntile(ttempdata[, variable2], nportf)
      old_i <- i+1
      i <- i + nrow(ttempdata) 
      intoutput[old_i:i,] <- data.matrix(ttempdata)
    }
    
  }
  

  crspd_ <- na.omit(intoutput)
  crspn7 <- data.frame(crspd_)
  crspn7$month <- 7
  crspn8 <- data.frame(crspd_)
  crspn8$month <- 8
  crspn9 <- data.frame(crspd_)
  crspn9$month <- 9
  crspn10 <- data.frame(crspd_)
  crspn10$month <- 10
  crspn11 <- data.frame(crspd_)
  crspn11$month <- 11
  crspn12 <- data.frame(crspd_)
  crspn12$month <- 12
  
  crspd_[,"year"] <- crspd_[,"year"]+1
  crspn1 <- data.frame(crspd_)
  crspn1$month <- 1
  crspn2 <- data.frame(crspd_)
  crspn2$month <- 2
  crspn3 <- data.frame(crspd_)
  crspn3$month <- 3
  crspn4 <- data.frame(crspd_)
  crspn4$month <- 4
  crspn5 <- data.frame(crspd_)
  crspn5$month <- 5
  crspn6 <- data.frame(crspd_)
  crspn6$month <- 6
  
  data_sorted <- rbind(crspn1, crspn2, crspn3, crspn4, crspn5, crspn6, crspn7, crspn8, crspn9, crspn10, crspn11, crspn12)
  return(data_sorted)
}

#END


###############################
### DEPENDENT DOUBLE SORTER ###
###############################
# [On time-variant variable]
# probably obsolete

dependent_double_sorter <- function(data, variable1, variable2, nportf){
  
  ##crsp_sb <- dependent_double_sorter(crspn, "decileME", "beta", 10)
  
  data <- data.matrix(data)
  
  data_sorted <- matrix(NA, nrow <- nrow(data)+1, ncol = ncol(data)+1)
  colnames(data_sorted) <- c(colnames(data), paste("decilecond", variable2, sep = ""))
  

  i <- 1
  for (yr in unique(data[,"year"])){
    for (mn in unique(data[,"month"])){
      for (d in 1:10){
        if (sum(((data[,"year"] == yr)&(data[,"month"] == mn)&(data[, variable1] == d))) == 0){
          next
        }
        tempdata <- data[(data[,"year"]==yr)&(data[,"month"]==mn)&(data[, variable1]==d),]
        tempdata <- cbind(tempdata, ntile(tempdata[, variable2], nportf))
        colnames(tempdata)[ncol(data)+1] <- paste("decilecond", variable2, sep = "")
        old_i <- i+1
        i <- i + nrow(tempdata) 
        data_sorted[old_i:i,] <- tempdata
        
      }
      
    }

  }
  
  data_sorted <- na.omit(data.frame(data_sorted))
  return(data_sorted)
}

#END




####################################################
### DAILY FAMA-MACBETH REGRESSION WITH SIMPLE SE ###
####################################################

# REQUIRES: nothing
# ARGUMENTS:
# data - data.table
# regressors <- vector of the colnames to be used as regressors
# depvar <- string, colname of dependent variable
# OUTPUT: df with 4 rows, reporting key statistics. To report results nicely, run dummy olses, use them as arguments for stargazer, specifying coef and se from this output df.

# it modifies input DT, so use copy(data) as an argument
# added on 06/30/19
# 08/18/19: .SDcols makes the function run for 20 sec instead of 5 sec...

fama_macbeth_simple_se_19_d <- function(data, regressors, depvar){
  
  ##fama_macbeth_simple_se_19_d(copy(crsp), c("beta", "size"), "RET")
  
  n_regressors <- length(regressors)
  regressors1 <- c("int", regressors)
  fm_results <- matrix(NA, nrow = 3, ncol = length(regressors)+1)
  rownames(fm_results) <- c("coefmean", "se", "tstat")
  data[,int:=1]
  
  temp <- data[, as.list(lm.fit(as.matrix(.SD[,1:(n_regressors+1)]), .SD[[n_regressors+2]])$coefficients), .SDcols = c(regressors1, depvar), by = prd]
  
  temp <- na.omit(temp)
  names(temp)[3:(length(regressors)+2)] <- regressors
  colnames(fm_results) <- c("Intercept", regressors)
  
  fm_results[1,] <- colMeans(temp)[2:ncol(temp)]
  fm_results[2,] <- (colSd(temp)[2:ncol(temp)])/(sqrt(nrow(temp)))
  fm_results[3,] <- fm_results[1,]/fm_results[2,]
  
  return(fm_results)
  
}

# END



### the version with R2 ###

# added on 06/30/19

fama_macbeth_simple_se_19_d_r2 <- function(data, regressors, depvar){
  
  ##fama_macbeth_simple_se_19_d(copy(crsp), c("beta", "size"), "RET")
  
  n <- length(regressors)
  regressors1 <- c("int", regressors)
  fm_results <- matrix(NA, nrow = 4, ncol = length(regressors)+1)
  rownames(fm_results) <- c("coefmean", "se", "tstat", "r2")
  data[,int:=1]
  
  temp <- data[, as.list(c(lm.fit(as.matrix(.SD[,1:(n+1)]), .SD[[n+2]])$coefficients, (((.SD[[n+2]]-mean(.SD[[n+2]]))^2)%*%(rep(1,.N)) - ((.SD[[n+2]]-lm.fit(as.matrix(.SD[,1:(n+1)]), .SD[[n+2]])$fitted.values)^2)%*%(rep(1,.N)))/ ((.SD[[n+2]]-mean(.SD[[n+2]]))^2)%*%(rep(1,.N)) )), .SDcols = c(regressors1, depvar), by = prd]
  
  temp <- na.omit(temp)
  names(temp)[3:(length(regressors)+3)] <- c(regressors, "r2")
  colnames(fm_results) <- c("Intercept", regressors)
  
  fm_results[1,] <- colMeans(temp)[2:(ncol(temp)-1)]
  fm_results[2,] <- (colSd(temp)[2:(ncol(temp)-1)])/(sqrt(nrow(temp)))
  fm_results[3,] <- fm_results[1,]/fm_results[2,]
  fm_results[4,] <- colMeans(temp)[ncol(temp)]
  
  return(fm_results)
  
}

# END



#########################################################
### DAILY FAMA-MACBETH REGRESSION WITH NEWEY-WEST SEs ###
#########################################################

# requires sandwitch.
# lags are calculated according to BEM14.
# added on 06/30/19

fama_macbeth_nw_se_19_d_r2 <- function(data, regressors, depvar){
  
  ##fama_macbeth_simple_se_19_d(copy(crsp), c("beta_bw", "size"), "RET")
  
  n <- length(regressors)
  regressors1 <- c("int", regressors)
  fm_results <- matrix(NA, nrow = 4, ncol = length(regressors)+1)
  rownames(fm_results) <- c("coefmean", "se", "tstat", "r2")
  data[,int:=1]
  
  temp <- data[, as.list(c(lm.fit(as.matrix(.SD[,1:(n+1)]), .SD[[n+2]])$coefficients, (((.SD[[n+2]]-mean(.SD[[n+2]]))^2)%*%(rep(1,.N)) - ((.SD[[n+2]]-lm.fit(as.matrix(.SD[,1:(n+1)]), .SD[[n+2]])$fitted.values)^2)%*%(rep(1,.N)))/ ((.SD[[n+2]]-mean(.SD[[n+2]]))^2)%*%(rep(1,.N)) )), .SDcols = c(regressors1, depvar), by = prd]
  
  temp <- na.omit(temp)
  names(temp)[3:(length(regressors)+3)] <- c(regressors, "r2")
  colnames(fm_results) <- c("Intercept", regressors)
  
  lagg <- ceiling((nrow(temp)/100)^(2/9)*4)
  for (i in 1:(n+1)){
    fm_results[2,i] <- sqrt(NeweyWest(lm(temp[[i+1]]~1), lag=lagg))
  }
  
  fm_results[1,] <- colMeans(temp)[2:(ncol(temp)-1)]
  fm_results[3,] <- fm_results[1,]/fm_results[2,]
  fm_results[4,] <- colMeans(temp)[ncol(temp)]
  
  return(fm_results)
  
}

# END




################################
### PORTFOLIOS FROM UNISORTS ###
################################

# the latest version to parameterize weighting variable and to get rid of _nt thing
# added on 03/10/19
# even given data.table for the heaviest calculations, this function is very slow...

characteristic_portfolios <- function(data, variable1, ret_variable, sort, weighting_type="value", grouping_var, weighting_var){
  
  ##unisorted_portfolios_gen(returns, "cmom11", 'RET', 5, "equal", "trdprd")
  ##characteristic_portfolios(crsp, "ind", 'RET', 49, "value", "prd", "laggedme")
  # it uses data.table syntax to make aggregation faster
  
  portfolio_list <- vector("list", length = prod(sort))
  
  data <- data[!is.na(data[, variable1]),]
  
  count <- 1
  for(i in 1:sort){
    
    name <- paste(variable1, ret_variable, i, sep="")
    tempdata <- data[(data[,variable1] == i),]
    tempdata <- data.table(tempdata)
    if (weighting_type == "value"){
      tempret <- (tempdata[, sum(eval(as.name(weighting_var))*eval(as.name(ret_variable)), na.rm=T), by=grouping_var])[,2]/(tempdata[,sum(eval(as.name(weighting_var)), na.rm=T), by=grouping_var][,2])
    } else if(weighting_type == "equal"){
      tempret <- (tempdata[, mean(eval(as.name(ret_variable)), na.rm=TRUE), by=grouping_var])[,2]
    } else {
      stop('"weighting_type" not correctly specified')
    }
    tempdata <- cbind(tempdata[,sum(eval(as.name(ret_variable)), na.rm=T), by=grouping_var][,1],  tempret)
    colnames(tempdata)[2] <- name
    portfolio_list[[count]] <- data.frame(tempdata)
    count <- count+1
    
  }
  
  portfolios <- Reduce(function(...) full_join(..., by=grouping_var), portfolio_list)
  portfolios <- portfolios[!is.na(portfolios[, grouping_var]),]
  portfolios <- portfolios[order(portfolios$prd),]
  
  return(portfolios)
  
}

#END



#####################
### DAILY VERSION ###
#####################

# this is not an ultimate portfolio generator, we can optimize the code further by leveraging on data.table more
# it will lose speed for very fine sorts
# added on 05/29/19
# need to check
# updated on 06/03/19 to enable other LHS
# updated on 06/07/20 to fix some error in R 3.6.3


characteristic_portfolios_d <- function(data, variable, ret_variable, weighting_type="value", weighting_var){
  
  ##characteristic_portfolios_d(crsp, "beta_nt_10", 'RET', "value", "lme")
  
  setorder(data, prd)
  varname <- gsub( "_nt_.*$", "", variable)
  setnames(data, variable, "varbl")
  setnames(data, ret_variable, "RETT")
  setnames(data, weighting_var, "wght")
  sort <- max(data[,varbl], na.rm = TRUE)
  data <- data[!is.na(data[, varbl])]
  portfolio_list <- vector("list", length = prod(sort))
  
  for(i in 1:sort){
    
    tdata <- data[varbl==i]
    if (weighting_type == "value"){
      tempret <- tdata[,.(pRET=sum(wght*RETT, na.rm=TRUE)/sum(wght, na.rm=TRUE)), by=prd]
    } else if(weighting_type == "equal"){
      tempret <- tdata[,.(pRET=mean(RETT, na.rm=TRUE)), by=prd]
    } else {
      stop('"weighting_type" not correctly specified')
    }
    
    setnames(tempret, "pRET", paste(varname, i, sep="_"))
    portfolio_list[[i]] <- tempret
    
  }
  
  portfolios <- Reduce(function(...) full_join(..., by="prd"), portfolio_list)
  setDT(portfolios)
  portfolios <- portfolios[!is.na(prd)]
  setorder(portfolios, prd)
  
  return(portfolios)
  
}

#END




#################################
### DAILY/MONTHLY BREAKPOINTS ###
#################################

# as long as period variable is called "prd", this function is agnostic of frequency.
# for 40M dCRSP it takes less than 6 sec
# added on 05/28/19
# updated on 08/31/19 to use 30-40-30 for tercile sort.

make_bpoints_d <- function(data, variable, sort){
  
  # breaks_beta_10 <- make_bpoints_d(crsp, "beta", 10)
  
  setnames(data, variable, "varbl")
  cnames <- c("prd", paste("br", 0:(sort), sep=""))
  if (sort!=3){
    temp <- data[, .(q_v = quantile(varbl, seq(0, 1, 1/sort), na.rm = T)), by=prd]
  }
  else{
    temp <- data[, .(q_v = quantile(varbl, c(0, 0.3, 0.7, 1), na.rm = T)), by=prd]
  }
  temp[, q_n:=rep(1:(sort+1),length(unique(prd)))]
  temp <- dcast.data.table(temp, prd~q_n, value.var = "q_v")
  
  names(temp) <- cnames
  temp[,br0 := min(br0, na.rm = TRUE)-10]
  temp[,ncol(temp):=max(data[,varbl], na.rm = TRUE)+10]
  setnames(data, "varbl", variable)
  temp <- temp[complete.cases(temp),]
  
  return(temp)
  
}

# END


##################################################
### DAILY BREAKPOINTS WITH MONTHLY REBALANCING ###
##################################################

# breakpoints are fixed as of the first days of the month
# this version modifies the object, passed as the first argument so copy it!
# this function is likely to spit out the warning that some columns do not exist, which is fine.
# added on 06/29/19

make_bpoints_d_monthlyreb <- function(data, variable, sort){
  
  # breaks_mom122_10m <- make_bpoints_d_monthlyreb(copy(crsp), "mom122", 10)
  
  setnames(data, variable, "varbl")
  cnames <- c("prd", paste("br", 0:(sort), sep=""))
  data[,c("ymprd"):=list(NULL)]
  
  prd_map <- fread("daily_prd_map.csv")
  prd_map[, date:=as.Date(date)]
  prd_map[,c("year", "month", "day"):=list(year(date), month(date), day(date))]
  prd_map[,date:=NULL]
  
  data <- prd_map[data, on="prd"]
  setorder(data, prd)
  invisible(gc)
  
  temp <- data[order(day),head(day,1),by=ymprd]
  setnames(temp, "V1", "day")
  temp <- temp[data, on=c("ymprd", "day"), nomatch=0]
  temp <- unique(temp, by=c("PERMNO", "prd"))
  invisible(gc())
  temp <- data[, .(q_v = quantile(varbl, seq(0, 1, 1/sort), na.rm = T)), by=ymprd]
  temp[, q_n:=rep(1:(sort+1),length(unique(ymprd)))]
  temp <- dcast.data.table(temp, ymprd~q_n, value.var = "q_v")
  names(temp) <- c("ymprd", cnames[2:length(cnames)])
  temp[,br0 := min(br0, na.rm = TRUE)-10]
  temp[,ncol(temp):=max(data[,varbl], na.rm = TRUE)+10]
  temp <- temp[complete.cases(temp),]
  temp <- temp[prd_map[,.(prd, ymprd)], on="ymprd", nomatch=0]
  temp[,ymprd:=NULL]
  setcolorder(temp, c("prd", names(temp)[1:(ncol(temp)-1)]))
  
  return(temp)
  
}

# END



######################
### ANNUAL VERSION ###
######################

# to be applied to annual data
# added on 03/10/2020

make_bpoints_a <- function(data, variable, sort){
  
  # breaks_beta_10 <- make_bpoints_a(cpst, "op", 10)
  
  setnames(data, variable, "varbl")
  cnames <- c("year", paste("br", 0:(sort), sep=""))
  if (sort!=3){
    temp <- data[, .(q_v = quantile(varbl, seq(0, 1, 1/sort), na.rm = T)), by=year]
  }
  else{
    temp <- data[, .(q_v = quantile(varbl, c(0, 0.3, 0.7, 1), na.rm = T)), by=year]
  }
  temp[, q_n:=rep(1:(sort+1),length(unique(year)))]
  temp <- dcast.data.table(temp, year~q_n, value.var = "q_v")
  
  names(temp) <- cnames
  temp[,br0 := min(br0, na.rm = TRUE)-10]
  temp[,ncol(temp):=max(data[,varbl], na.rm = TRUE)+10]
  setnames(data, "varbl", variable)
  temp <- temp[complete.cases(temp),]
  
  return(temp)
  
}

# END




###############################
### DAILY/MONTHLY UNISORTER ###
###############################

# added on 05/28/19

unisorter_var_d <- function(data, variable, bks){
  
  ##crsp <- unisorter_var_d(crsp, "beta", breaks_beta_10)
  
  setnames(data, variable, "varbl")
  bks <- bks[, data.frame(embed(unlist(.SD),2)[,2:1]), by=prd]
  bks[, ntile := seq_len(.N), by=prd]
  data[bks, on=c("prd"="prd","varbl>=X1","varbl<=X2"), ntile := i.ntile]
  setnames(data, "varbl", variable)
  setnames(data, "ntile", paste(variable, "nt", max(data[,ntile], na.rm = TRUE), sep="_"))
  
  return(data)
  
}

#END



########################
### ANNUAL UNISORTER ###
########################

#to be applied to annual data
# added on 03/10/2020

unisorter_var_a <- function(data, variable, bks){
  
  ##crsp <- unisorter_var_a(cpst, "op", breaks_10)
  
  setnames(data, variable, "varbl")
  bks <- bks[, data.frame(embed(unlist(.SD),2)[,2:1]), by=year]
  bks[, ntile := seq_len(.N), by=year]
  data[bks, on=c("year"="year","varbl>=X1","varbl<=X2"), ntile := i.ntile]
  setnames(data, "varbl", variable)
  setnames(data, "ntile", paste(variable, "nt", max(data[,ntile], na.rm = TRUE), sep="_"))
  
  return(data)
  
}

#END




####################################
### PORTFOLIOS FROM DOUBLE SORTS ###
####################################


doublesorted_portfolios <- function(data, variable1, variable2, ret_variable, sort){
  
  ##portfolios_mebm <- doublesorted_portfolios(cpcrsp, "me", "bm", "RET", c(2,3))
  # it uses data.table syntax to make aggregation faster
  
  variable1full <- paste("nt", variable1, sort[1], sep='_')
  variable2full <- paste("nt", variable2, sort[2], sep='_')
  portfolio_list <- vector("list", length = prod(sort))
  
  count <- 1
  for(i in 1:sort[1]){
    for (j in 1:sort[2]){
      name <- paste(variable1, i, variable2, j, ret_variable, sep="")
      tempdata <- data[(data[,variable1full] == i)&(data[,variable2full] == j),]
      tempdata <- data.table(tempdata)
      tempret <- (tempdata[, mean(me*eval(as.name(ret_variable))), by='year,month'])[,3]/(tempdata[,mean(me), by='year,month'][,3])
      tempdata <- cbind(tempdata[,mean(eval(as.name(ret_variable))), by='year,month'][,1:2],  tempret)
      colnames(tempdata)[3] <- name
      portfolio_list[[count]] <- data.frame(tempdata)
      count <- count+1
    }
  }
  
tempfactors <- Reduce(function(...) full_join(..., by=c("year", "month")), portfolio_list)
return(tempfactors)
  
}

#END


# upgdared fully dt version, added on 06/12/19:

doublesorted_portfolios_dt <- function(data, variable1, variable2, ret_variable="RET", sort){
  
  ##portfolios_mebm <- doublesorted_portfolios_dt(copy(crsp), "me", "bm", "RET", c(2,3))
  # it uses data.table syntax to make aggregation faster
  
  variable1full <- paste(variable1, "nt", sort[1], sep='_')
  variable2full <- paste(variable2, "nt", sort[2], sep='_')
  portfolio_list <- vector("list", length = prod(sort))
  
  count <- 1
  for(i in 1:sort[1]){
    for (j in 1:sort[2]){
      name <- paste(variable1, i, variable2, j, ret_variable, sep="")
      tempdata <- data[eval(as.name(variable1full)) == i][eval(as.name(variable2full)) == j]
      tempret <- tempdata[, .(name=mean(lme*eval(as.name(ret_variable)), na.rm = TRUE)/mean(lme, na.rm = TRUE)), by=prd]
      setnames(tempret, "name", name)
      portfolio_list[[count]] <- tempret
      count <- count+1
    }
  }
  
  tempfactors <- Reduce(function(...) full_join(..., by=c("prd")), portfolio_list)
  return(data.table(tempfactors))
  
}

#END



####################################
### PORTFOLIOS FROM TRIPLE SORTS ###
####################################


triplesorted_portfolios <- function(data, variable1, variable2, variable3, ret_variable, sort){
  
  ##portfolios_mebm <- doublesorted_portfolios(cpcrsp, "me", "bm", "RET", c(2,3))
  # it uses data.table syntax to make aggregation faster
  
  variable1full <- paste("nt", variable1, sort[1], sep='_')
  variable2full <- paste("nt", variable2, sort[2], sep='_')
  variable3full <- paste("nt", variable3, sort[3], sep='_')
  portfolio_list <- vector("list", length = prod(sort))
  
  count <- 1
  for(i in 1:sort[1]){
    for (j in 1:sort[2]){
      for (k in 1:sort[3]){
        
        name <- paste(variable1, i, variable2, j, variable3, k, ret_variable, sep="")
        tempdata <- data[(data[,variable1full] == i)&(data[,variable2full] == j)&(data[,variable3full] == k),]
        tempdata <- data.table(tempdata)
        tempret <- (tempdata[, mean(me*eval(as.name(ret_variable))), by='year,month'])[,3]/(tempdata[,mean(me), by='year,month'][,3])
        tempdata <- cbind(tempdata[,mean(eval(as.name(ret_variable))), by='year,month'][,1:2],  tempret)
        colnames(tempdata)[3] <- name
        portfolio_list[[count]] <- data.frame(tempdata)
        count <- count+1
        
      }
      
    }
  }
  
  tempfactors <- Reduce(function(...) full_join(..., by=c("year", "month")), portfolio_list)
  return(tempfactors)
  
}

#END




# new version, added on 03/22/19
# prd instead of year-month

BJS_regression_unisorted_19 <- function(testrets, ret_var, factrets, regressors){
  
  ##temptable <- BJS_regression_unisorted(treturns5, "eRET_vw", factors6, c("Mkt.RF", "SMB", "HML"))
  
  var1 <- colnames(testrets)[2]
  var1p <- length(unique(testrets[,var1]))
  
  retlist <- list()
  for(i in 1:var1p){
    tempdata <- testrets[(testrets[, var1]==i),]
    tempdata <- inner_join(tempdata, factrets, by=c("prd"))
    retlist[[(length(retlist)+1)]] <- tempdata
  }
  
  
  listlg <- 5 + 2*(length(regressors))
  
  bjs_table <- data.frame(matrix(NA, nrow = var1p, ncol = listlg-1))
  
  if (listlg == 7){
    colnames(bjs_table) <- c(var1, "c_int", "t_int", paste("c_", regressors[1], sep=""), paste("t_", regressors[1], sep=""), "R^2")
    
    for (i in 1:length(retlist)){
      df <- retlist[[i]]
      model <- lm(get(ret_var)~get(regressors), df)
      bjs_table[i,] <- c(df[1, var1], model$coefficients[1], summary(model)[["coefficients"]][,"t value"]["(Intercept)"], model$coefficients[2], summary(model)[["coefficients"]][,"t value"][2], summary(model)$adj.r.squared)    
    }
    return(bjs_table)
  }
  
  if (listlg == 9){
    colnames(bjs_table) <- c(var1, "c_int", "t_int", paste("c_", regressors[1], sep=""), paste("t_", regressors[1], sep=""), paste("c_", regressors[2], sep=""), paste("t_", regressors[2], sep=""), "R^2")
    
    for (i in 1:length(retlist)){
      df <- retlist[[i]]
      model <- lm(get(ret_var)~get(regressors[1])+get(regressors[2]), df)
      bjs_table[i,] <- c(df[1, var1], model$coefficients[1], summary(model)[["coefficients"]][,"t value"]["(Intercept)"], model$coefficients[2], summary(model)[["coefficients"]][,"t value"][2], model$coefficients[3], summary(model)[["coefficients"]][,"t value"][3], summary(model)$adj.r.squared)    
    }
    return(bjs_table)
  }
  
  if (listlg == 11){
    colnames(bjs_table) <- c(var1, "c_int", "t_int", paste("c_", regressors[1], sep=""), paste("t_", regressors[1], sep=""), paste("c_", regressors[2], sep=""), paste("t_", regressors[2], sep=""), paste("c_", regressors[3], sep=""), paste("t_", regressors[3], sep=""), "R^2")
    
    for (i in 1:length(retlist)){
      df <- retlist[[i]]
      model <- lm(get(ret_var)~get(regressors[1])+get(regressors[2])+get(regressors[3]), df)
      bjs_table[i,] <- c(df[1, var1], model$coefficients[1], summary(model)[["coefficients"]][,"t value"]["(Intercept)"], model$coefficients[2], summary(model)[["coefficients"]][,"t value"][2], model$coefficients[3], summary(model)[["coefficients"]][,"t value"][3], model$coefficients[4], summary(model)[["coefficients"]][,"t value"][4], summary(model)$adj.r.squared)    
    }
    return(bjs_table)
    
  }
  
  if (listlg == 13){
    colnames(bjs_table) <- c(var1, "c_int", "t_int", paste("c_", regressors[1], sep=""), paste("t_", regressors[1], sep=""), paste("c_", regressors[2], sep=""), paste("t_", regressors[2], sep=""), paste("c_", regressors[3], sep=""), paste("t_", regressors[3], sep=""), paste("c_", regressors[4], sep=""), paste("t_", regressors[4], sep=""), "R^2")
    
    for (i in 1:length(retlist)){
      df <- retlist[[i]]
      model <- lm(get(ret_var)~get(regressors[1])+get(regressors[2])+get(regressors[3])+get(regressors[4]), df)
      bjs_table[i,] <- c(df[1, var1], model$coefficients[1], summary(model)[["coefficients"]][,"t value"]["(Intercept)"], model$coefficients[2], summary(model)[["coefficients"]][,"t value"][2], model$coefficients[3], summary(model)[["coefficients"]][,"t value"][3], model$coefficients[4], summary(model)[["coefficients"]][,"t value"][4], model$coefficients[5], summary(model)[["coefficients"]][,"t value"][5], summary(model)$adj.r.squared)    
    }
    return(bjs_table)
    
  }
  
  if (listlg == 15){
    colnames(bjs_table) <- c(var1, "c_int", "t_int", paste("c_", regressors[1], sep=""), paste("t_", regressors[1], sep=""), paste("c_", regressors[2], sep=""), paste("t_", regressors[2], sep=""), paste("c_", regressors[3], sep=""), paste("t_", regressors[3], sep=""), paste("c_", regressors[4], sep=""), paste("t_", regressors[4], sep=""), paste("c_", regressors[5], sep=""), paste("t_", regressors[5], sep=""), "R^2")
    
    for (i in 1:length(retlist)){
      df <- retlist[[i]]
      model <- lm(get(ret_var)~get(regressors[1])+get(regressors[2])+get(regressors[3])+get(regressors[4])+get(regressors[5]), df)
      bjs_table[i,] <- c(df[1, var1], model$coefficients[1], summary(model)[["coefficients"]][,"t value"]["(Intercept)"], model$coefficients[2], summary(model)[["coefficients"]][,"t value"][2], model$coefficients[3], summary(model)[["coefficients"]][,"t value"][3], model$coefficients[4], summary(model)[["coefficients"]][,"t value"][4], model$coefficients[5], summary(model)[["coefficients"]][,"t value"][5], model$coefficients[6], summary(model)[["coefficients"]][,"t value"][6], summary(model)$adj.r.squared)    
    }
    return(bjs_table)
    
  }
  
  if (listlg == 17){
    colnames(bjs_table) <- c(var1, "c_int", "t_int", paste("c_", regressors[1], sep=""), paste("t_", regressors[1], sep=""), paste("c_", regressors[2], sep=""), paste("t_", regressors[2], sep=""), paste("c_", regressors[3], sep=""), paste("t_", regressors[3], sep=""), paste("c_", regressors[4], sep=""), paste("t_", regressors[4], sep=""), paste("c_", regressors[5], sep=""), paste("t_", regressors[5], sep=""), paste("c_", regressors[6], sep=""), paste("t_", regressors[6], sep=""), "R^2")
    
    for (i in 1:length(retlist)){
      df <- retlist[[i]]
      model <- lm(get(ret_var)~get(regressors[1])+get(regressors[2])+get(regressors[3])+get(regressors[4])+get(regressors[5])+get(regressors[6]), df)
      bjs_table[i,] <- c(df[1, var1], model$coefficients[1], summary(model)[["coefficients"]][,"t value"]["(Intercept)"], model$coefficients[2], summary(model)[["coefficients"]][,"t value"][2], model$coefficients[3], summary(model)[["coefficients"]][,"t value"][3], model$coefficients[4], summary(model)[["coefficients"]][,"t value"][4], model$coefficients[5], summary(model)[["coefficients"]][,"t value"][5], model$coefficients[6], summary(model)[["coefficients"]][,"t value"][6], model$coefficients[7], summary(model)[["coefficients"]][,"t value"][7], summary(model)$adj.r.squared)    
    }
    return(bjs_table)
    
  }
  
  if (!(listlg %in% c(7, 9, 11, 13, 15, 17))){
    stop("too many factors")
  }
  
  
}

#END




################################
### UNIVERSAL BJS REGRESSION ###
################################

# this is more efficient version: period-agnostic, allowing for any number of parameters and producing RNM-style output
# notice that due to string formatting, needed to get [] for t-stat, I can not parameterize digits in stargazer. need to set up digits here.
# I assume that the last column in "data" is L/S returns
# added on 06/01/2019

BJS_regs <- function(data, factors, long_output = TRUE, digits_ = 3){
  
  ##BJS_regs(ptfs_vw, factors, TRUE, 2)
  
  setcolorder(factors, c(names(factors)[names(factors)!="prd"], "prd"))
  regdata <- factors[data, on="prd", nomatch=0]
  n_portf <- ncol(data)-1
  n_factors <- ncol(factors)-1
  
  if(long_output == TRUE){
    
    output <- data.frame(matrix(NA, n_portf*2, 4 + n_factors))
    names(output) <- c(c("Ntile", "Ret", "Alpha"), names(factors)[1:n_factors], "adjR2")
    
    for (i in 1:n_portf){
      
      ols_eq <- as.formula(paste(names(regdata)[n_factors+1+i], paste(names(factors)[1:n_factors], collapse = " + "), sep = " ~ "))
      tempmod <- summary(lm(ols_eq, regdata))
      output[(2*i-1):(2*i), 1] <- i
      output[(2*i-1):(2*i), 2] <- c(regdata[,mean(eval(as.name(names(regdata)[n_factors+1+i])), na.rm = TRUE)], sqrt(nrow(regdata))*regdata[,mean(eval(as.name(names(regdata)[n_factors+1+i])), na.rm = TRUE)/sd(eval(as.name(names(regdata)[n_factors+1+i])), na.rm = TRUE)])   
      output[2*i-1, 3:(ncol(output)-1)] <- (tempmod$coefficients)[,1]
      output[2*i, 3:(ncol(output)-1)] <- (tempmod$coefficients)[,3]
      output[(2*i-1):(2*i), ncol(output)] <- tempmod$adj.r.squared
      
    }
    
    output[(2*n_portf-1):(2*n_portf), 1] <- "LS"
    output[,-1] <- format(round(output[,-1], digits_), digits_)
    
    for(i in seq(2,2*n_portf,2)){
      output[i,2:(n_factors+3)] <- lapply(output[i,2:(n_factors+3)], function(x) paste("[ ", x, "]", sep=""))
      output[i,1] <- ""
      output[i,(n_factors+4)] <- ""
    }
    
  }
  
  else{
    
    output <- data.frame(matrix(NA, 2, n_portf+2))
    names(output) <- c(c("Model", "Stat"), names(data)[-1])
    
    for (i in 1:n_portf){
      tempmod <- summary(lm(get(names(regdata)[4+i])~EMKT+SMB+HML, regdata))
      output[1, (i+2)] <- (tempmod$coefficients)[1,1]
      output[2, (i+2)] <- (tempmod$coefficients)[1,3]
    }
    
    output[1:2,2] <- c("Alpha", "T_stat")
    output[2,1] <- ""
    output[1:2, 3:(n_portf+2)] <- format(round(output[1:2, 3:(n_portf+2)], digits_), digits_)
    output[2, 3:(n_portf+2)] <- lapply(output[2, 3:(n_portf+2)], function(x) paste("[ ", x, "]", sep=""))
    
  }
  
  return(output)
  
}

#END


###########################
### FAST DATE FORMATTER ###
###########################

# create joiner and join it to character dates whenever needed
# added on 04/05/2021

y <- seq(as.Date('1900-01-01'), Sys.Date(), by = 'day')
id_date <- data.table(id = as.character(y), fdate = as.Date(y), key = 'id')
setnames(id_date, "id", "date")

#set.seed(128)
#x <- as.character(Sys.Date()-sample(40000, 1e6, TRUE))
##system.time(date4 <- id_date[setDT(list(id = x)), on='id', date])

# use id_date as a joiner

# END



###########################################
### CONVERTER OF CPST INTO MONTHLY DATA ###
###########################################

# added on 04/03/2019
# updated on 04/04/2019
# data must have prd, year and month variables
# WARNING: if you use this function to produce ff-style lags, smth like cpst_prd <- make_cpst_prd_19(copy(cpst), 128), 
#          then you can get duplicated values of accounting variables in July. to fix that, do smth like  > cpst_prd <- cpst_prd[!duplicated(cpst_prd, by=c("PERMNO", "prd"), fromLast=TRUE)].

make_cpst_prd <- function(data, ignore_fmonth = TRUE, lag = 6){
  
  ## cpst_prd <- make_cpst_prd(cpst, ignore_fmonth = FALSE, lag=6)
  ## lag=6 means that fyear1962 data is matched to returns starting from July 1963
  
  df <- data.frame(matrix(NA, (2019-1961)*12, 1))
  df[,1] <- rep(1962:2019, each=12)
  colnames(df) <- "year"
  df$month <- rep(1:12, 2019-1961)
  df$prd_j <- df$year*12 + df$month - 23500 + 12 + lag
  
  if (ignore_fmonth == FALSE){
    data <- data[data$month == 12,]
  }

  joiner <- df[, c("year", "prd_j")]
  data <- left_join(joiner, data, by = c("year")) 
  setnames(data, "prd", "prd_acc")
  setnames(data, "prd_j", "prd")
  data[,c("year", "month")] <- list(NULL)
  
  return(data)
  
}

# END


# new dt version of this function
# WARNING: the output will have around 80M rows, so cpst should have only the variables you really need
# added on 06/28/19
# update on 09/16/20: "daily_prd_map" is overwritten, there may be discrepancies in date format.


make_cpst_prd_19 <- function(data, lag = 126){
  
  ## cpst_prd <- make_cpst_prd(copy(cpst), lag=126)
  ## lag=126 means that fyear1962 data is matched to returns starting from July 1963
  lag <- lag+1
  
  joiner <- fread("daily_prd_map.csv")
  joiner[,date:=mdy(date)]
  joiner[,year:=year(date)]
  years <- joiner[,.(year, prd)]
  setnames(years, c("year", "prd"), c("lyear", "lprd"))
  setorder(joiner, prd)
  joiner[, lprd:=c(rep(NA, 252+lag), prd[1:(.N-as.integer(252+lag))])]
  joiner <- years[joiner, on="lprd", nomatch=0]
  joiner <- joiner[,.(prd, lyear)]
  setnames(joiner, "lyear", "year")

  data <- data[joiner, on="year", allow.cartesian=TRUE]
  data <- data[!is.na(PERMNO)]
  setorder(data, PERMNO, prd)
  setnames(data, "year", "cpst_year")
  
  return(data)
  
}

# END


# map for joining cpst with monthly crsp. 
# notice that here I update cpst data across firms monthly.

make_cpst_prd_mon_monupd <- function(data){
  
  ##cpst_prd <- make_cpst_prd_mon_monupd(copy(cpst))
  
  joiner <- data.table(rel_prd = rep(45:740, each=12), prd = rep(0:11, 741-45))
  joiner[,prd:=rel_prd+prd]
  
  data <- joiner[data, on="rel_prd", allow.cartesian = TRUE]
  setorder(data, PERMNO, rel_prd)
  setkey(data, PERMNO, prd)
  dups <- duplicated(data, fromLast = TRUE, by=key(data))
  data <- data[!dups]
  
  return(data)
  
}

# END


###############
### TRIMMER ###
###############

# input is assumed to be data.table
# added on 06/03/19
# not sure whether it is useful

trim <- function(data, variable, prob){
  setnames(data, variable, "varbl")
  data <- data[varbl>quantile(varbl, prob[1])][varbl<quantile(varbl, prob[2])]
  setnames(data, "varbl", variable)
  return(data)
}

#END


trim1 <- function(data, variable, prob){
  setnames(data, variable, "varbl")
  values <- data[,varbl]>quantile(data[,varbl], prob[1])&data[,varbl]<quantile(data[,varbl], prob[2])
  setnames(data, "varbl", variable)
  return(values)
}

#END

# these functions are not yet nicely applicable



####################################
### ### ### CRSP CLEANER ### ### ###
####################################

# start_year > 1962

crsp_cleaner <- function(data, start_year=1973, end_year=2018, minPRC=0){
  
  ##
  
  crsp <- fread(data)
  invisible(gc())
  crsp <- data.frame(crsp)
  crsp <- crsp[3:nrow(crsp), ]
  crsp$RET <- as.numeric(crsp$RET)*100
  crsp$vwretd <- as.numeric(crsp$vwretd)*100
  crsp$ewretd <- as.numeric(crsp$ewretd)*100
  crsp$sprtrn <- as.numeric(crsp$sprtrn)*100
  crsp$PRC <- abs(crsp$PRC)
  crsp <- crsp[crsp$PRC>minPRC,]
  crsp$date <- dmy(crsp$date)
  crsp$year <- year(crsp$date)
  crsp$month <- month(crsp$date)
  crsp <- crsp[(crsp$year>=start_year)&(crsp$year<=end_year),]
  dim(crsp)
  length(unique(crsp$PERMNO))
  
  crsp <- crsp[((crsp$SHRCD == 10)|(crsp$SHRCD == 11)),]
  crsp <- crsp[((crsp$EXCHCD == 1)|(crsp$EXCHCD == 2)|(crsp$EXCHCD == 3)),]
  dim(crsp)
  
  crsp <- crsp[,c("PERMNO", "year", "month", "RET",  "vwretd", "ewretd", "sprtrn", "PRC", "SHROUT")]
  crsp <- na.omit(crsp)
  crsp <- distinct(crsp)
  dim(crsp)
  crsp$prd <- (crsp$year*12+crsp$month)-min(crsp$year)*12
  
  rfr <- read.csv("rfr.csv")
  rfr$date <- mdy(rfr$MCALDT)
  rfr$MCALDT <- NULL
  rfr$year <- year(rfr$date)
  rfr$month <- month(rfr$date)
  crsp$me <- crsp$PRC*crsp$SHROUT/1000
  crsp <- crsp[!is.na(crsp$me),]
  crsp <- crsp[order(crsp$prd),]
  crsp <- crsp[,c("PERMNO", "year", "month", "prd", "RET",  "vwretd", "ewretd", "sprtrn", "me")]
  
  crsp <- inner_join(crsp, rfr, by = c("year", "month"))
  invisible(gc())
  crsp <- crsp[,c("PERMNO", "year", "month", "prd", "RET",  "vwretd", "ewretd", "sprtrn", "TMYTM", "me")]
  crsp <- distinct(crsp)
  crsp <- na.omit(crsp)
  crsp <- crsp[order(crsp$PERMNO, crsp$year, crsp$month),]
  dim(crsp)
  crsp$rfr <- crsp$TMYTM/(12)
  crsp$TMYTM <- NULL
  crsp$eRET <- crsp$RET - crsp$rfr
  crsp$mkt <- crsp$vwretd - crsp$rfr
  
  # adding momentum
  
  mom72 <- fread("crsp_63_mom72.csv")
  mom72 <- data.frame(mom72)
  mom72 <- mom72[,c("PERMNO", "year", "month", "mom72RET")]
  crsp <- inner_join(crsp, mom72, by=c("PERMNO", "year", "month"))
  
  mom122 <- fread("crsp_63_mom122.csv")
  mom122 <- data.frame(mom122)
  mom122 <- mom122[,c("PERMNO", "year", "month", "mom122RET")]
  crsp <- inner_join(crsp, mom122, by=c("PERMNO", "year", "month"))
  
  mom127 <- fread("crsp_63_mom127.csv")
  mom127 <- data.frame(mom127)
  mom127 <- mom127[,c("PERMNO", "year", "month", "mom127RET")]
  crsp <- inner_join(crsp, mom127, by=c("PERMNO", "year", "month"))
  
  rm("mom72", "mom122", "mom127", "rfr")
  invisible((gc))
  
  return(crsp)
  
}

#END




##############################
### ULTIMATE MAKE_MOMENTUM ###
##############################

# added on 05/03-04/2019

# This should be the most efficient make_momentum possible w/o using C.
# minrets is the minimum number of nonmissing RET in data. Usually set it to 121 or 252. if too small, you can get an error (usually) or sometimes just more than order of magnitude longer execution time.
# b_f_60tt contains complete temporary function with speed becnhmarks and complete comparison with old versions


# in the future, want to do smth similar to speed up rolling regressions. can use roll_regres from here: 
# https://stackoverflow.com/questions/12139934/is-there-a-fast-way-to-run-a-rolling-regression-inside-data-table

# restricting the data to three variables speeds things by the factor of two, can implement it later.
# implementing this is not trivial, since data.table uses assignment by reference, i.e., deleting all other columns will delete mom too.

# this function works well for relatively short windows (up to 20-50), but for longer windows need the different function without hole-filling.


make_momentum_19f <- function(data, lag1, lag2, minrets){
  
  ##temp21n <- make_momentum_19f(crsp, 21, 1, 121)
  
  momentum_variable <- paste0("mom", lag1, lag2)
  
  fdata <- data[, CJ(prd=min(prd):max(prd)), by = "PERMNO"]
  fdata <- data[fdata, on=c("PERMNO", "prd")]
  # this fills holes. order of magnitude speed gain vs dplyr/tidyverse.

  temp <- fdata[, .(nna=sum(!is.na(RET))) , by="PERMNO"]
  permnos <- temp[nna<minrets,PERMNO]
  fdata <- fdata[!(PERMNO %fin% permnos)]
  # this excludes firms with few nnms RET

  # The following code implements very fast rolling product, suggested here: https://stackoverflow.com/questions/30603999/how-do-i-take-a-rolling-product-using-data-table
  N <- as.integer(lag1)
  fdata <- fdata[, mom := Reduce(`*`, shift(RET, as.integer(lag2):(N), type = "lag")), by = "PERMNO"]
  # this produces more than 2 orders of magnitude speed gain over probably the most effecient dplyr solution
  
  colnames(fdata)[colnames(fdata)=="mom"] <- momentum_variable
  newdata <- fdata[!is.na(fdata$RET),]
  
  return (newdata)
}

#END


# functional programming in data.table is tricky. Maybe clean solution for this problem does not even exist.
# useful trick to avoid using actual functional programming, while achieving exactly the same results: 
# just rename the column to the fixed name before and after doing computations.
# ideally, I should write the general function now, but i am a bit worried about making RET RET_index and then renaming it back - 
# there is some chance to mess things up. but mostly i am just too lazy.

make_momentum_19f_index <- function(data, lag1, lag2, varname, minrets){
  
  ##temp21n <- make_momentum_19f(crsp, 21, 1, 121)
  
  setnames(data, varname, "RET_index")
  
  momentum_variable <- paste0("mom_", varname, lag1, lag2)
  
  fdata <- data[, CJ(prd=min(prd):max(prd)), by = "PERMNO"]
  fdata <- data[fdata, on=c("PERMNO", "prd")]
  # this fills holes. order of magnitude speed gain vs dplyr/tidyverse.
  
  temp <- fdata[, .(nna=sum(!is.na(RET_index))) , by="PERMNO"]
  permnos <- temp[nna<minrets,PERMNO]
  fdata <- fdata[!(PERMNO %fin% permnos)]
  # this excludes firms with few nnms RET
  
  # The following code implements very fast rolling product, suggested here: https://stackoverflow.com/questions/30603999/how-do-i-take-a-rolling-product-using-data-table
  N <- as.integer(lag1)
  fdata <- fdata[, mom := Reduce(`*`, shift(RET_index, as.integer(lag2):(N), type = "lag")), by = "PERMNO"]
  # this produces more than 2 orders of magnitude speed gain over probably the most effecient dplyr solution
  
  colnames(fdata)[colnames(fdata)=="mom"] <- momentum_variable
  newdata <- fdata[!is.na(fdata$RET_index),]
  setnames(newdata, "RET_index", varname)
  
  return (newdata)
}

#END


### This momentum calculator allows for some number of holes, no larger than minrets-(lag1-lag2).
### It is very convenient for long windows. 
### it is the fastest momentum_calculator as of 08/21/19.
# modified to be consistent with make_momentum_19f_index on 09/01/19.
# modified to include window_scalar on 10/01/19.

make_long_momentum <- function(data, lag1, lag2, minrets, varname, window_scalar = 1.1){
  
  ##tcrsp <- make_long_momentum(copy(crsp), 252*3, 0, 252*3, "RET")
  
  setnames(data, varname, "RET_index")
  momentum_variable <- paste0("mom_", varname, lag1, lag2)
  setorder(data, PERMNO, prd)
  
  temp <- data[, .(nna=sum(!is.na(RET_index))) , by="PERMNO"]
  permnos <- temp[nna<minrets,PERMNO]
  data <- data[!(PERMNO %fin% permnos)]
  
  data <- data[, mom := Reduce(`*`, shift(RET_index, as.integer(lag2):as.integer(lag1), type = "lag")), by = "PERMNO"]
  
  # below I exclude the observations where the window, used to calculate momentum, was too long
  data[, lprd:=c(rep(NA, lag1), prd[1:(.N-as.integer(lag1))]), by=PERMNO]
  data[, window:=prd-lprd]
  data <- data[window<=(lag1+1-lag2)*window_scalar]
  data[,c("lprd", "window"):=list(NULL)]
  
  colnames(data)[colnames(data)=="mom"] <- momentum_variable
  newdata <- data[!is.na(data$RET_index),]
  setnames(newdata, "RET_index", varname)
  setnames(data, "RET_index", varname)
  
  return (newdata)
}

#END



#############################
### DAILY BETA CALCULATOR ###
#############################

# this will be fast universal beta calculator. It was developed in script "daily_beta_calc_11".
# window is always rolling
# the output inludes only observations with proper beta
# whole daily CRSP will run in less than 1 min
# on beta estimation: Alan in FIN418, l4s7 claims that we use ols with intercept to estimate beta.
# added on 05/27/19
# updated on 06/02/19 to use only past data to calculate beta. negligible effect.
# updated on 06/03/19 to do the same for r2
# updated on 06/04/19 to correct lagger
# as of 06/05 dependent and independent variables are hardcoded in beta and roll_sd functions. parameterizing them should be trivial.
# on holes: I estimate betas using chronologically-ordered data with rows = window and then see the length of the window in days. this length should not exceed max_window.

make_beta <- function(data, window=252, clean_window=378, max_window=278, r2 = FALSE, only_past = TRUE){
  
  ##crsp <- make_beta(crsp)
  
  window <- as.integer(window)
  data <- clean_crsp_by_obs(data, clean_window)
  setkey(data, PERMNO, prd)
  
  if(r2 == FALSE){
    data[, beta := roll_regres.fit(x = cbind(1, .SD[["vwretd"]]), y = .SD[["RET"]], do_downdates = TRUE, width = window)$coefs[, 2], by = PERMNO]
  }
  else{
    data[, c("beta", "r2") := list(roll_regres.fit(x = cbind(1, .SD[["vwretd"]]), y = .SD[["RET"]], do_downdates = TRUE, width = window, do_compute = c("r.squareds"))$coefs[, 2],
                                   roll_regres.fit(x = cbind(1, .SD[["vwretd"]]), y = .SD[["RET"]], do_downdates = TRUE, width = window, do_compute = c("r.squareds"))$r.squareds), by = PERMNO]
  }
  
  invisible(gc())
  
  tdata <- data[,.(PERMNO, prd, beta)]
  if(only_past == TRUE){
    tdata[, beta:=c(NA, beta[-.N]), by=PERMNO]
    if (r2 == TRUE){
      data[, r2:=c(NA, r2[-.N]), by=PERMNO]
    }
  }
  rows <- which(!is.na(tdata$beta))
  tdata[rows, windoww := tdata[rows, prd]-tdata[rows-(window-1),prd]]
  tdata <- tdata[windoww<max_window]
  tdata$windoww <- NULL
  
  data[,beta:=NULL]
  data <- tdata[data,on=c("PERMNO", "prd"), nomatch=0]
  data <- unique(data)
  
  return(data)
  
}

#END


#dt <- data.table(PERMNO=c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2), prd=c(1,2,3,7,8,9,10,11,12,13,14,15, 1,2,11,12,13,14,15,16,17,18,19,20), RET = c(1,2,4,1.5,4,2,6,3,9,-2,5,0,1,11,0,-4,2,7,1,4,7,1,1,3), vwretd=c(1,2,3,1,5,0,-1,1,2,-1,3,1,0,-2,1,3,1,4,1,2,7,0,2,1))
#make_beta(dt,3,3,3)


### Beta calculator with parameterized dependent and independent variables ###

# added on 06/06/19
# parameterized depvar on 06/18/19
# edited on 12/15/20

make_beta_param <- function(data, variable="vwretd", retvar = "RET", window=252, clean_window=378, max_window=378, r2 = FALSE, only_past = TRUE){
  
  ##crsp <- make_beta(crsp)
  
  window <- as.integer(window)
  data <- clean_crsp_by_obs(data, clean_window)
  setkey(data, PERMNO, prd)
  
  if(r2 == FALSE){
    data[, beta := roll_regres.fit(x = cbind(1, .SD[[variable]]), y = .SD[[retvar]], do_downdates = TRUE, width = window)$coefs[, 2], by = PERMNO]
  }
  else{
    data[, c("beta", "r2") := list(roll_regres.fit(x = cbind(1, .SD[[variable]]), y = .SD[[retvar]], do_downdates = TRUE, width = window, do_compute = c("r.squareds"))$coefs[, 2],
                                   roll_regres.fit(x = cbind(1, .SD[[variable]]), y = .SD[[retvar]], do_downdates = TRUE, width = window, do_compute = c("r.squareds"))$r.squareds), by = PERMNO]
  }
  
  invisible(gc())
  
  tdata <- data[,.(PERMNO, prd, beta)]
  setorder(tdata, PERMNO, prd)
  if(only_past == TRUE){
    tdata[, beta:=c(NA, beta[-.N]), by=PERMNO]
    if (r2 == TRUE){
      data[, r2:=c(NA, r2[-.N]), by=PERMNO]
    }
  }
  rows <- which(!is.na(tdata$beta))
  tdata[rows, windoww := tdata[rows, prd]-tdata[rows-(window-1),prd]]
  tdata <- tdata[windoww<max_window]
  tdata[,windoww := NULL]
  
  data[,beta:=NULL]
  data <- tdata[data,on=c("PERMNO", "prd"), nomatch=0]
  data <- unique(data)
  
  return(data)
  
}

#END


### Beta calculator with parameterized dependent and independent variables and controls ###

# added on 07/15/19
# data should include the controls

make_beta_param_controls <- function(data, variable, retvar = "RET", window=252, clean_window=378, max_window=278, r2 = FALSE, only_past = TRUE, controls){
  
  ##crsp <- make_beta_param_controls(copy(crsp0), variable="q97_pu", window=24, clean_window=24, max_window=32, r2=FALSE, controls=c("vwretd", "SMB", "HML"))
  
  window <- as.integer(window)
  data <- clean_crsp_by_obs(data, clean_window)
  setkey(data, PERMNO, prd)
  data[,int:=1]
  n_regressors <- length(c(variable, controls))+1
  
  
  if(r2 == FALSE){
    data[, beta := roll_regres.fit(x = as.matrix(.SD[,1:n_regressors]), y = .SD[[n_regressors+1]], do_downdates = TRUE, width = window)$coefs[, 2], .SDcols = c("int", variable, controls, retvar), by = PERMNO]
  } else{
    data[, c("beta", "r2") := list(roll_regres.fit(x = as.matrix(.SD[,1:n_regressors]), y = .SD[[n_regressors+1]], do_downdates = TRUE, width = window)$coefs[, 2],
                                   roll_regres.fit(x = as.matrix(.SD[,1:n_regressors]), y = .SD[[n_regressors+1]], do_downdates = TRUE, width = window, do_compute = c("r.squareds"))$r.squareds),
         .SDcols = c("int", variable, controls, retvar), by = PERMNO]
  }

  invisible(gc())
  
  tdata <- data[,.(PERMNO, prd, beta)]
  if(only_past == TRUE){
    tdata[, beta:=c(NA, beta[-.N]), by=PERMNO]
    if (r2 == TRUE){
      data[, r2:=c(NA, r2[-.N]), by=PERMNO]
    }
  }
  rows <- which(!is.na(tdata$beta))
  tdata[rows, windoww := tdata[rows, prd]-tdata[rows-(window-1),prd]]
  tdata <- tdata[windoww<max_window]
  tdata$windoww <- NULL
  
  data[,c("beta", "int"):=list(NULL)]
  data <- tdata[data,on=c("PERMNO", "prd"), nomatch=0]
  data <- unique(data)
  
  return(data)
  
}

#END



#######################
### IVOL CALCULATOR ###
#######################

# based on make_beta, it calculated iVOL as sd of residuals from BJS regression on factors (capm, ff3 or ff5) rolling over 1 month.
# it is coded in seemingly stupid way, but more elegant solution (.SDcols, similar to fama_mac) is thrice slower.
# minret is the minimum number of nonhissing days in a month. period is already shifted forward after the application of this function.
# beware: it is slow. It runs for 5-8 minutes on daily crsp with 4-5 variables.
# added on 07/02/19
# updated on 12/24/19

make_ivol <- function(data, model="ff3", minret=15){
  
  ff <- unique(fread("FF7_daily.csv", select = c("prd", "EMKT", "SMB", "HML", "RMW", "CMA", "Mom", "ST_Rev")))
  ff[,int:=1]
  
  if(model=="capm"){
    regressors <- c("int", "EMKT")
    data <- ff[, c("prd", regressors), with=FALSE][data, on="prd", nomatch=0]
    temp <- data[,sd(lm.fit(cbind(int, EMKT), RET)$residuals),by=c("PERMNO", "ymprd")]
    
  }
  if(model=="ff3"){
    regressors <- c("int", "EMKT", "SMB", "HML")
    data <- ff[, c("prd", regressors), with=FALSE][data, on="prd", nomatch=0]
    temp <- data[,sd(lm.fit(cbind(int, EMKT, SMB, HML), RET)$residuals),by=c("PERMNO", "ymprd")]
    
  }
  if(model=="ff5"){
    regressors <- c("int", "EMKT", "SMB", "HML", "RMW", "CMA")
    data <- ff[, c("prd", regressors), with=FALSE][data, on="prd", nomatch=0]
    data <- data[!is.na(RMW)]
    temp <- data[,sd(lm.fit(cbind(int, EMKT, SMB, HML, RMW, CMA), RET)$residuals),by=c("PERMNO", "ymprd")]
  }
  
  setnames(temp, "V1", "ivol")
  goodobs <- data[,.N, by=c("PERMNO", "ymprd")]
  goodobs <- goodobs[N>=minret]
  temp <- goodobs[,.(PERMNO, ymprd)][temp, on=c("PERMNO", "ymprd"), nomatch=0]
  temp[,ymprd:=ymprd+1]
  
  return(temp)

}

# END



# updated on 06/16/20:

make_ivol_20 <- function(data, model="ff3", minret=15){
  
  ff <- unique(fread("dFF7_19.csv", select = c("prd", "EMKT", "SMB", "HML", "RMW", "CMA", "Mom", "ST_Rev")))
  ff[,int:=1]
  
  if(model=="capm"){
    regressors <- c("int", "EMKT")
    data <- ff[, c("prd", regressors), with=FALSE][data, on="prd", nomatch=0]
    temp <- data[,sd(lm.fit(cbind(int, EMKT), RET)$residuals),by=c("PERMNO", "ymprd")]
    
  }
  if(model=="ff3"){
    regressors <- c("int", "EMKT", "SMB", "HML")
    data <- ff[, c("prd", regressors), with=FALSE][data, on="prd", nomatch=0]
    temp <- data[,sd(lm.fit(cbind(int, EMKT, SMB, HML), RET)$residuals),by=c("PERMNO", "ymprd")]
    
  }
  if(model=="ff5"){
    regressors <- c("int", "EMKT", "SMB", "HML", "RMW", "CMA")
    data <- ff[, c("prd", regressors), with=FALSE][data, on="prd", nomatch=0]
    data <- data[!is.na(RMW)]
    temp <- data[,sd(lm.fit(cbind(int, EMKT, SMB, HML, RMW, CMA), RET)$residuals),by=c("PERMNO", "ymprd")]
  }
  
  setnames(temp, "V1", "ivol")
  goodobs <- data[,.N, by=c("PERMNO", "ymprd")]
  goodobs <- goodobs[N>=minret]
  temp <- goodobs[,.(PERMNO, ymprd)][temp, on=c("PERMNO", "ymprd"), nomatch=0]
  temp[,ymprd:=ymprd+1]
  
  return(temp)
  
}

# END


####################
### ROLLING MEAN ###
####################

# rolling mean by PERMNO, rolling over prd.
# this function is developed from make_beta above
# added on 06/06/19

make_rollingmean <- function(data, variable, window, clean_window, max_window, only_past = TRUE){
  
  ##temp4 <- make_rollingmean(temp3, "ac", 16, 16, 252*3, TRUE)
  
  window <- as.integer(window)
  data <- clean_crsp_by_obs(data, clean_window)
  setkey(data, PERMNO, prd)
  
  data[, rmean := rollmean(eval(as.name(variable)), k = window, fill = NA, align = "right"), by = PERMNO]
  invisible(gc())
  
  tdata <- data[,.(PERMNO, prd, rmean)]
  if(only_past == TRUE){
    tdata[, rmean:=c(NA, rmean[-.N]), by=PERMNO]
  }
  rows <- which(!is.na(tdata$rmean))
  tdata[rows, windoww := tdata[rows, prd]-tdata[rows-(window-1),prd]]
  tdata <- tdata[windoww<max_window]
  tdata$windoww <- NULL
  
  data[,rmean:=NULL]
  data <- tdata[data,on=c("PERMNO", "prd"), nomatch=0]
  data <- unique(data)
  
  return(data)
  
}

#END



#####################################
### ROLLING VOLATILITY CALCULATOR ###
#####################################

# added on 06/05/19

make_rollsd <- function(data, variable = "RET", window=252, clean_window=378, max_window=278, only_past = TRUE){
  
  ##crsp2 <- make_rollsd(copy(crsp), "RET", window = 48, clean_window = 48, max_window = 3*252)
  
  window <- as.integer(window)
  data <- clean_crsp_by_obs(data, clean_window)
  setkey(data, PERMNO, prd)

  data[, roll_sd := roll_sd(eval(as.name(variable)), window, fill=NA, align="right"), by="PERMNO"]
  invisible(gc())
  
  tdata <- data[,.(PERMNO, prd, roll_sd)]
  if(only_past == TRUE){
    tdata[, roll_sd:=c(NA, roll_sd[-.N]), by=PERMNO]
  }
  rows <- which(!is.na(tdata$roll_sd))
  tdata[rows, windoww := tdata[rows, prd]-tdata[rows-(window-1),prd]]
  tdata <- tdata[windoww<max_window]
  tdata$windoww <- NULL
  data[,roll_sd:=NULL]
  data <- tdata[data,on=c("PERMNO", "prd"), nomatch=0]
  data <- unique(data)
  
  return(data)
  
}

#END




###########################################################################
### FUNCTION, CALCULATING SAMPLE MEAN RETURNS AND THEIR STANDARD ERRORS ###
###########################################################################
# [for now - regular standard errors]
# can implement NW se from sandwich via NeweyWest(), but it will be better to write my own NW se estimator: NeweyWest(model1, lag = 6, adjust = T)
# assume the first two columns are year, month

average_returns_calculator <- function(data, type="t_stat"){
  
  output <- matrix(NA, ncol = ncol(data), nrow = 2)
  
  if(type=="se"){
    for (i in 1:ncol(data)){
      tempmodel <- lm(data[,i]~1)
      output[,i] <- coef(summary(tempmodel))[1:2]
    }
    
    output <- data.frame(output)
    rownames(output) <- c("Mean", "SE")
  }
  

  else{
    for (i in 1:ncol(data)){
      tempmodel <- lm(data[,i]~1)
      output[,i] <- coef(summary(tempmodel))[1:2]
    }
    
    output <- data.frame(output)
    output[2,] <- output[1,]/output[2,]
    rownames(output) <- c("Mean", "T_stat")
  }
  
  return(output)
  
}

#END



######################
### SUE CALCULATOR ###
######################
# 09/28/18

sue_calculator <- function(data, minannounc=6, mininnov=4){
  
  ## cpst_sue <- sue_calculator(cpst)
  
  output <- data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)+3))
  colnames(output) <- c(colnames(data), "qprd", "sue", "ue")
  
  j <- 1
  
  for (firm in unique(data$gvkey)){
    
    tempdata <- data[data$gvkey == firm,]
    tempdata$qprd <- tempdata$fyearq*4 + tempdata$fqtr - 7900
    tempdata <- tempdata[order(tempdata$qprd),]
    tempdata$SUE <- NA
    if(nrow(tempdata)<6){
      next
    }
    tempdata$l4eps <- c(NA, NA, NA, NA, tempdata$epspxq[1:(length(tempdata$epspxq)-4)])
    tempdata$UE <- tempdata$epspxq - tempdata$l4eps
    for(i in 6:nrow(tempdata)){
      if(sum(!is.na(tempdata$epspxq[(max(1,(i-8))):i]))<minannounc){
        next
      }
      if(sum(!is.na(tempdata$UE[(max(1,(i-8))):i]))<mininnov){
        next
      }
      tempdata$SUE[i] <- (tempdata$UE[i])/(sd(tempdata$UE[(max(1,(i-7))):i], na.rm = T))
    }
    tempdata$l4eps <- NULL
    
    
    output[j:(j+nrow(tempdata)-1),] <- tempdata
    j <- j+nrow(tempdata)
    
  }
  
  output <- output[!is.na(output$year),]
  return(output)
  
}

#END




#############################
### FF 49 INDUSRTY SORTER ###
#############################
# cwd should contain ff49_table
# updated dt version, added on 12/4/19

ff49_industries_sorter <- function(data, industryname = "SICCD"){
  
  industries <- fread("FF49_ind.csv")
  setnames(industries, "SICCD", industryname)
  data <- industries[data, on=industryname]

  return(data)
  
}

#END




#####################################
### Correlation matrix with stars ###
#####################################
# stolen from: https://gist.github.com/anhqle/5855936

corstarsl <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .01, "***", ifelse(p < .05, "** ", ifelse(p < .1, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars 
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}
#END



#################################
### SAMPLE RESTRICTER BY MCAP ###
#################################

# for every month, it restricts the sample to PERMNOs with inflation-adjusted ME above _ as of the first day of the month
# lower-frequency restriction is to mitigate obvious bias towards firms with high returns in the previous period.
# added on 06/06/19
# as of v.96, 08/09 I am not sure whether this is correct, see line 15.

restrict_crsp <- function(data, mcap=1000000){
  
  data <- unique(data, by=c("PERMNO", "prd"))
  cal_map <- fread("daily_calendar_map.csv")
  cal_map$cdate <- as.Date(cal_map$cdate)
  cal_map <- unique(cal_map, by="prd")
  data <- cal_map[data, on="prd", nomatch=0]
  invisible(gc())
  data[, c("year", "month", "day") := list(year(cdate), month(cdate), day(cdate))]
  data[, ymprd := year*12+month-23500]
  data[,c("year", "month"):=NULL]
  setorder(data, ymprd)
  invisible(gc())
  
  firstdays <- data[order(day),head(day,1),by=ymprd]
  # should it be by=c("PERMNO", "ymprd) ???
  setnames(firstdays, "V1", "day")
  firstdays <- firstdays[data, on=c("ymprd", "day"), nomatch=0]
  firstdays <- unique(firstdays, by=c("PERMNO", "prd"))
  invisible(gc())
  
  cpi <- fread("CPI47.csv")
  cpi <- cpi[,.(date=mdy(DATE), cpi=CPIAUCSL)]
  cpi[,c("year", "month"):=list(year(date), month(date))]
  cpi[,ymprd:=year*12+month-23500]
  cpi <- unique(cpi)
  firstdays <- cpi[,.(ymprd, cpi)][firstdays, on="ymprd"]
  # lets set cpi as of Dec18, i.e., 252.72 as a numeraire
  firstdays[,cpi_adj:=252.72/cpi]
  
  firstdays[, me_adj:=me*cpi_adj]
  firstdays <- firstdays[me_adj>mcap]
  invisible(gc())
  permnos <- firstdays[,.(PERMNO, ymprd)]
  data <- permnos[data, on=c("PERMNO", "ymprd"), nomatch = 0]
  data <- unique(data, by=c("PERMNO", "prd"))
  data <- data[!((PERMNO==18040)&(ymprd==59))]
  # excluding the typo (or an anomaly) in CRSP
  invisible(gc())
  
  return(data)
  
}

#END



### at monthly frequency (annual rebalancing, i.e. it restricts the sample to PERMNO-years with inflation-adjusted ME above _ as of the first month of the year) ###

# added on 08/09/19
# mcap is in millions here

restrict_crsp_monthly <- function(data, mcap=20){
  
  firstdays <- data[order(prd),head(prd,1),by=c("PERMNO", "year")]
  setnames(firstdays, "V1", "prd")
  firstdays <- firstdays[data, on=c("PERMNO", "year", "prd"), nomatch=0]
  firstdays <- unique(firstdays, by=c("PERMNO", "prd"))
  invisible(gc())
  
  cpi <- fread("CPI47.csv")
  cpi <- cpi[,.(date=mdy(DATE), cpi=CPIAUCSL)]
  cpi[,c("year", "month"):=list(year(date), month(date))]
  cpi[,ymprd:=year*12+month-23500]
  cpi <- unique(cpi)
  setnames(cpi, "ymprd", "prd")
  firstdays <- cpi[,.(prd, cpi)][firstdays, on="prd"]
  # lets set cpi as of Dec18, i.e., 252.72 as a numeraire
  firstdays[,cpi_adj:=252.72/cpi]
  
  firstdays[, me_adj:=me*cpi_adj]
  firstdays <- firstdays[me_adj>mcap]
  invisible(gc())
  permnos <- firstdays[,.(PERMNO, year)]
  data <- permnos[data, on=c("PERMNO", "year"), nomatch = 0]
  data <- unique(data, by=c("PERMNO", "prd"))
  data <- data[!((PERMNO==18040)&(prd==59))]
  # excluding the typo (or an anomaly) in CRSP
  invisible(gc())
  
  return(data)
  
}

#END



# added on 01/06/20
# generalized version of the function above: restricts the sample on any variable.
# WARNING: to avoid bias, use lagged variable to restrict on.

restrict_crsp_monthly_20 <- function(data, variable, minvalue, cpiadjust = FALSE){
  
  firstdays <- data[order(prd),head(prd,1),by=c("PERMNO", "year")]
  setnames(firstdays, "V1", "prd")
  firstdays <- firstdays[data, on=c("PERMNO", "year", "prd"), nomatch=0]
  firstdays <- unique(firstdays, by=c("PERMNO", "prd"))
  invisible(gc())
  
  cpi <- fread("CPI47.csv")
  cpi <- cpi[,.(date=mdy(DATE), cpi=CPIAUCSL)]
  cpi[,c("year", "month"):=list(year(date), month(date))]
  cpi[,ymprd:=year*12+month-23500]
  cpi <- unique(cpi)
  setnames(cpi, "ymprd", "prd")
  firstdays <- cpi[,.(prd, cpi)][firstdays, on="prd"]
  # lets set cpi as of Dec18, i.e., 252.72 as a numeraire
  firstdays[,cpi_adj:=252.72/cpi]
  
  if(cpiadjust == TRUE){
    firstdays[, variable_adj:=get(variable)*cpi_adj]
  }
  else {
    firstdays[, variable_adj:=get(variable)]
  }
  firstdays <- firstdays[variable_adj>minvalue]
  invisible(gc())
  permnos <- firstdays[,.(PERMNO, year)]
  data <- permnos[data, on=c("PERMNO", "year"), nomatch = 0]
  data <- unique(data, by=c("PERMNO", "prd"))
  data <- data[!((PERMNO==18040)&(prd==59))]
  # excluding the typo (or an anomaly) in CRSP
  invisible(gc())
  
  return(data)
  
}

#END



#####################################
### PAST TS VOLATILITY CALCULATOR ###
#####################################

# added on 08/21/19

make_tsvol <- function(data, lagg = 126, max_window){
  
  ## tcrsp <- make_tsvol(copy(crsp[,.(PERMNO, prd, RET)]), lagg=126, max_window=138)
  
  setorder(data, PERMNO, prd)
  data[, roll_sd := roll_sd(RET, lagg, fill=NA, align="right"), by="PERMNO"]
  # this thing is slow, but can not find better solution. and it is not a signal (includes current period). 
  data[, lwprd:=shift(prd, lagg-1, "lead"), by="PERMNO"]
  data[, window:=prd-lwprd, by="PERMNO"]
  data[window>max_window, roll_sd:=NA]
  data[, roll_sd:=shift(roll_sd, 1, "lead"), by="PERMNO"]
  # so now roll_sd is a proper signal, known at the beginning of the period
  data[,c("lwprd", "window"):=list(NULL)]
  
  return(data)
  
}

#END




indnames <- c("Agric", "Food", "Soda", "Beer", "Smoke", "Toys", "Fun", "Books", "Hshld", "Clths",
              "Hlth", "MedEq", "Drugs", "Chems", "Rubbr", "Txtls", "BldMt", "Cnstr", "Steel", "FabPr",
              "Mach", "ElcEq", "Autos", "Aero", "Ships", "Guns", "Gold", "Mines", "Coal", "Oil",
              "Util", "Telcm", "PerSv", "BusSv", "Hardw", "Softw", "Chips", "LabEq", "Paper", "Boxes",
              "Trans", "Whlsl", "Rtail", "Meals", "Banks", "Insur", "RlEst", "Fin", "Other")

ff49inds <- data.table(indnumber=1:49, indnames)




################################
### BETA ADJUSTER (WELCH 19) ###
################################

# the code is copypaste from Ivo`s paper. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3371240
# added on 12/12/19.

wins.rel <- function( r, rmin, rmax ) {
  rlower <- ifelse( (rmin<rmax), rmin, rmax )
  rupper <- ifelse( (rmin<rmax), rmax, rmin )
  ifelse( r<rlower, rlower, ifelse(r>rupper, rupper, r) )
}

# END




#######################
### Mode calculator ###
#######################

# copypaste from https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# END

### ### END OF BASIC FUNCTIONS ### ### 

