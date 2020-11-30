
if (!require(metafor)){
  install.packages("metafor")
}

library(metafor)

# load skin_result_inter

data <- read.table("/Users/ligk2e/Desktop/meta/new/skin_result_inter.txt",sep=" ")
row_na = which(is.na(data$V15))
data_clean = data[-(row_na),]

f <- function(row){
  T <- row[seq(2,length(row)-1,3)]
  V <- row[seq(3,length(row)-1,3)]
  i <- row[16]
  random <- tryCatch(rma(T,V),error=function(e) NULL)
  t_combine <- tryCatch(as.numeric(random$beta),error=function(e) NULL)
  v_combine <- tryCatch(as.numeric(random$se),error=function(e) NULL)
  p_combine <- tryCatch(as.numeric(random$pval),error=function(e) NULL)
  I_combine <- tryCatch(as.numeric(random$I2),error=function(e) NULL)
  
  new <- c(t_combine,v_combine,p_combine,I_combine,i)
  return(new)
}


common <- read.table('/Users/ligk2e/Desktop/meta/new/skin_common.txt',sep=" ")
com_data <- cbind(data,common)
result <- apply(com_data[1:17535,],1,f)
output <- matrix(unlist(result),ncol=5,byrow = T)

write.table(output,file='/Users/ligk2e/Desktop/meta/new/skin_result_R.txt',sep=" ")



# let's analyze glioblastoma

data <- read.table("/Users/ligk2e/Desktop/meta/new/glioblastoma/brain_result_inter.txt",sep=" ")
common <- read.table('/Users/ligk2e/Desktop/meta/new/glioblastoma/brain_common.txt',sep=" ")


f <- function(row){
  T <- row[seq(2,length(row)-1,3)]
  V <- row[seq(3,length(row)-1,3)]
  i <- row[10]
  random <- tryCatch(rma(T,V),error=function(e) NULL)
  t_combine <- tryCatch(as.numeric(random$beta),error=function(e) NULL)
  v_combine <- tryCatch(as.numeric(random$se),error=function(e) NULL)
  p_combine <- tryCatch(as.numeric(random$pval),error=function(e) NULL)
  I_combine <- tryCatch(as.numeric(random$I2),error=function(e) NULL)
  
  new <- c(t_combine,v_combine,p_combine,I_combine,i)
  return(new)
}



com_data <- cbind(data,common)
result <- apply(com_data[1:19493,],1,f)
output <- matrix(unlist(result),ncol=3,byrow = T)

write.table(output,file='/Users/ligk2e/Desktop/meta/new/skin_result_R.txt',sep=" ")

T <- as.numeric(data[1,])[seq(2,length(data[1,]),3)]
V <- as.numeric(data[1,])[seq(3,length(data[1,]),3)]
random <- rma(T,V)














