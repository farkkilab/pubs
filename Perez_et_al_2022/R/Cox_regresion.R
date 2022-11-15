###### Function to do Cox_regression analysis of each variable in a data_set ########
Cox_regresion_variables <- function(data, variables, formula='Surv(PFI.time, PFI)~'){

      univ_formulas <- sapply(variables,
                              function(x) as.formula(paste(formula, x)))

      univ_models <- lapply(univ_formulas, function(x){
        return(tryCatch(coxph(x, data = data, method = "exact"), error = function(e) {ret <<- 0}))
        })

      univ_results <- lapply(univ_models,
                             function(x){
                               if (length(x) > 1){
                               x <- summary(x)
                               positive <- x$nevent
                               prop.positive <- (x$nevent)/x$nevent
                               p.value <- signif(x$wald["pvalue"], digits=2)
                               wald.test <- signif(x$wald["test"], digits=2)
                               lastrow <- nrow(x$coef) #When using several covariables in the same regression, take the last one used
                               beta <- signif(x$coef[lastrow, 1], digits=2);#coeficient beta
                               HR <- signif(x$coef[lastrow, 2], digits=2);#exp(beta)
                               HR.confint.lower <- signif(x$conf.int[lastrow,"lower .95"], 2)
                               HR.confint.upper <- signif(x$conf.int[lastrow,"upper .95"],2)
                               #CI <- paste0(HR.confint.lower, "-", HR.confint.upper)
                               #res<-c(beta, HR, CI, wald.test, p.value)
                               #names(res)<-c("beta", "HR", "CI for HR (lower and upper 95%)", "wald.test","p.value")
                               res<-c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, p.value)
                               names(res)<-c("beta", "HR", "lowerCI_HR", "upperCI_HR", "wald.test","p.value")
                               return(res)
                               }else{
                                 res2 <- c(NA, NA, NA, NA, NA, NA)
                                 names(res2)<-c("beta", "HR", "lowerCI_HR", "upperCI_HR", "wald.test","p.value")
                                 return(res2)
                               }
                             })
      res <- t(as.data.frame(univ_results, check.names = FALSE))
      result <- as.data.frame(res)
      return(result)
}

###### Function to do Cox_regression analysis of each variable in a data_set and return values for table ########
Cox.regresion.variables2 <- function(data, variables, formula='Surv(PFI.time, PFI)~'){

   N <- sapply(variables, function(x){length(which(data[,x] == 1))})
   N <- as.data.frame(N)
   Prop <- paste0(round((N[,1] * 100 /nrow(data)),0), "%")

   univ_formulas <- sapply(variables,
                           function(x) as.formula(paste(formula, x)))

   univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data, method = "exact")})

   univ_results <- lapply(univ_models,
                          function(x){
                             x <- summary(x)
                             positive <- x$nevent
                             prop.positive <- paste0(round((x$nevent) * 100 /x$n,0),"%")
                             p.value <- signif(x$wald["pvalue"], digits=2)
                             wald.test <- signif(x$wald["test"], digits=2)
                             lastrow <- nrow(x$coef) #When using several covariables in the same regression, take the last one used
                             beta <- signif(x$coef[lastrow, 1], digits=3);#coeficient beta
                             HR <- signif(x$coef[lastrow, 2], digits=3);#exp(beta)
                             HR.confint.lower <- signif(x$conf.int[lastrow,"lower .95"], 2)
                             HR.confint.upper <- signif(x$conf.int[lastrow,"upper .95"],2)
                             res <- c(HR, p.value)
                             names(res)<-c("HR", "p.value")
                             return(res)
                          })
   res <- t(as.data.frame(univ_results, check.names = FALSE))
   result <- as.data.frame(res)
   result2 <- cbind(N, Prop, result)

   return(result2)
}
