#'Functions to use in the validation samples


#This function take a data frame, with values for each sample. A patient could have several samples.
#This function average the indicated values indicated per patient and return the resultant dataframe
avg_by_sample <- function(info, patient_column = 1, sample_column = 0, columns_to_avg = c(6,7)){
  mean_info_patient <- NULL
  #Select unique patient (or ID)
  for(j in unique(info[,patient_column])){
    patient.info <- info[info[,patient_column] %in% j,]
    if (sample_column){
      #Select by unique sub sample ID, condition or group
      for(h in unique(patient.info[,sample_column])){
        patient.info.sample <- patient.info[patient.info[,sample_column] %in% h,]
        mean_vals <- NULL
        for (colm in columns_to_avg){
          mean_vals <- c(mean_vals, mean(patient.info.sample[,colm]))
        }
        patient.info.sample[1,c(columns_to_avg)] <- mean_vals
        mean_info_patient <- rbind(mean_info_patient, patient.info.sample[1,])
      }
    }else{
      #If not sub sample ID, then average just by patient ID
      mean_vals <- NULL
      for (colm in columns_to_avg){
        mean_vals <- c(mean_vals, mean(patient.info[,colm]))
      }
      #Select just the first row
      patient.info[1,c(columns_to_avg)]  <-  mean_vals
      mean_info_patient <- rbind(mean_info_patient, patient.info[1,])
    }
  }
  return(mean_info_patient)
}


#This function take a data frame, with values for each sample. A patient could have several samples.
#This function average the indicated values indicated per patient and return the resultant dataframe


