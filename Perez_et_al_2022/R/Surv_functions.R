library(survival)
require("survminer")

#Function to create several Kaplan-Meir plots an arrange those in one Figure.
make.merge.survplots <- function(survobject_input, data_inputset, variables=c("myriadHRDhigh", "newHRDhigh", "tacayaHRDhigh"),  break_time=25, xlabplot="PFI time (months)", ylabplot="PFI probability", palette = c("#FA5A41", "#34A0D3"), max.time=FALSE){
    #User can provide an input max time
    if (max.time == FALSE){
      max.time = max(survobject_input[,1])
    }
    univ_formulas <- lapply(variables,
                          function(x) as.formula(paste0("survobject_input ~",  x)))

    fits <- surv_fit(univ_formulas, data = data_inputset)

    splots <- lapply(fits, function(x){
                ggsurvplot(x, data = data_inputset, pval = TRUE, risk.table = TRUE, xlim = c(0,max.time),
                      risk.table.y.text.col = TRUE, tables.y.text = FALSE, palette = palette,
                      break.time.by = break_time, conf.int = TRUE) + xlab(xlabplot) + ylab(ylabplot)
    })

    lapply(fits, function(x){
      pval <- surv_pvalue(x, data = data_inputset, method ="survdiff")
      print(pval)
      mediansurvival <- surv_median(x)
      print(mediansurvival)
    })

    plot <- arrange_ggsurvplots(splots, print = FALSE,
                                ncol = length(variables), nrow = 1, risk.table.height = 0.35)
    return(plot)
}


library(survival)
require("survminer")

#Function to create several Kaplan-Meir plots an arrange those in one Figure.
median.surv.times <- function(survobject_input, data_inputset, variables=c("myriadHRDhigh", "newHRDhigh", "tacayaHRDhigh")){
  univ_formulas <- lapply(variables,
                          function(x) as.formula(paste0("survobject_input ~",  x)))

  fits <- surv_fit(univ_formulas, data = data_inputset)
  medians <- lapply(fits, function(x){
                    l <- surv_median(x)
                    return(l$median)})
  meds <- as.data.frame(medians)
  marker.class <- c(0,1)
  meds <- cbind(meds,marker.class)
  return(meds)
}
