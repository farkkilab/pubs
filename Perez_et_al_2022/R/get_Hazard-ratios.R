#######################################################################################################################################
############################################# Function to make binarization of data and get HazardR  ##########################################
#######################################################################################################################################

#The next function first take as input the data tables with the survival information and ovaHRDscar values and Telli2016 scars
#The input tables (one for each cohort) have common names
#For ovaHRDscar levels the column name correspond to: HRDsum
#For Telli2016 HRD-AIs leves the column name correspond to: pHRDsum
#Other important columns names are:
####BRCAness: BRCA1/2mutation or deletion status (1 positive, 0 negative)

getHRs <- function(dataf, formulaHR ='Surv(OS.time, OS)~', datalabel="PanCanAtlas", returndata=FALSE, BRCAness.comp=TRUE){
  dataf$ovaHRDscar <- ifelse(dataf$HRDsum >= 54,1,0) #Cutoff For ovaHRDscar
  dataf$Telli2016.54 <- ifelse(dataf$pHRDsum >= 54,1,0) #For Telli2016, cutoff of 54
  dataf$Telli2016 <- ifelse(dataf$pHRDsum >= 42,1,0) #For Telli2016
  dataf$Takaya2020 <- ifelse(dataf$pHRDsum >= 63,1,0) #For Tacaya2020

  #Generating dataframe 2, without BRCAmut/del samples
  dataf2 <- dataf[dataf$BRCAness == 0,]
  #Change to positive status if there is not BRCAness
  dataf2$ovaHRDscar.BRCAwt <- ifelse(dataf2$ovaHRDscar == 1 & dataf2$BRCAness == 0,1,0)
  dataf2$Telli2016.BRCAwt <- ifelse(dataf2$Telli2016 == 1 & dataf2$BRCAness == 0, 1,0)
  dataf2$Takaya2020.BRCAwt <- ifelse(dataf2$Takaya2020 == 1 & dataf2$BRCAness == 0, 1,0)
  dataf2$Telli2016.54.BRCAwt <- ifelse(dataf2$Telli2016.54 == 1 & dataf2$BRCAness == 0,1,0) #For Telli2016, cutoff of 54

  ### Getting the hazard ratios for each the next variables/columns in dataf(2) ####
  if(BRCAness.comp){
    selectedvars <- c("BRCAness", "Telli2016", "Telli2016.54","Takaya2020", "ovaHRDscar")
  }else{
    selectedvars <- c("Telli2016", "Telli2016.54", "Takaya2020", "ovaHRDscar")
  }
  vars2 <- c("Telli2016.BRCAwt", "Telli2016.54.BRCAwt","Takaya2020.BRCAwt", "ovaHRDscar.BRCAwt")
  surv <- Cox_regresion_variables(dataf, variables =selectedvars, formula =formulaHR)
  surv2 <- Cox_regresion_variables(dataf2, variables =vars2, formula =formulaHR)


  names(surv) <- c("beta", "mean", "lower", "upper", "wald.test","p.value")
  surv <- surv[c(2:4,6)]
  surv$CI <- paste0(surv[,1], "(CI:",surv[,2], "-", surv[,3], ")")

  names(surv2) <- c("beta", "mean", "lower", "upper", "wald.test","p.value")
  surv2 <- surv2[c(2:4,6)]
  surv2$CI <- paste0(surv2[,1], "(CI:",surv2[,2], "-",surv2[,3], ")")

  survOSmix <- rbind(surv, surv2)
  nvalues <- c(NA, NA, NA) #This is important for plotting
  survOSmix <- rbind(nvalues, survOSmix)

  ## Count the number of samples selected according to each criteria, make a matrix with that info ####
  count_selected <- sapply(selectedvars, function(x){length(which(dataf[,x] == 1))})
  prop_selected <- paste(as.character(round(count_selected * 100 /nrow(dataf))),"%", sep = "")
  count_selected2 <- sapply(vars2, function(x){length(which(dataf2[,x] == 1))})
  prop_selected2 <- paste(as.character(round(count_selected2 * 100 /nrow(dataf2))), "%", sep="")

  data_counts_selected <- data.frame(names=c(selectedvars,  vars2),
                                     selected=c(count_selected, count_selected2),
                                     proportions=c(prop_selected, prop_selected2))

  column_labels <- c("Criteria", datalabel, "n selected", "prop selected")
  data_counts_selected <- rbind(column_labels, data_counts_selected)

  data_counts_selected_matrix <- as.matrix(data_counts_selected)

  #Return the hazzard ratios by default, but if returndata==1, then return the dataf table
  result <- list(survOSmix, data_counts_selected_matrix)
  if (returndata){
    return(dataf)
  }else{
    return(result)
  }
}

