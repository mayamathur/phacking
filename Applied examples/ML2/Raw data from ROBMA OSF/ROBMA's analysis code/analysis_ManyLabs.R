### prepare data
# read the datasets
df <- read.csv("Examples - data/data_ManyLabs.csv", comment.char="#")

# test types
table(df$test.method[!duplicated(df$study.id)])


# use t-statistics and DF with (Welch) two sample t-tests
df1 <- df[df$test.method %in% c("Welch Two Sample t-test", " Two Sample t-test"), ]
table(df1$test.method[!duplicated(df1$study.id)])

# join the sample sizes
df1$n <- df1$stat.n1 + df1$stat.n2
# remove unnecessary columns + some renaming
df1 <- df1[, c("study.id", "test.statistic", "n")]
colnames(df1) <- c("name", "t", "n")
df1$name <- as.character(df1$name)

# save the dataset
saveRDS(df1, "Examples - data/data_ManyLabs_t.RDS")


# and Fisher's z transformed correlation coefficents with the rest
df2 <- df[!df$test.method %in% c("Welch Two Sample t-test", " Two Sample t-test"), ]
table(df2$test.method[!duplicated(df2$study.id)])

df2 <- df2[, c("study.id", "ESCI.r", "stat.N")]
colnames(df2) <- c("name", "r", "n")
df2$name <- as.character(df2$name)

df2 <- na.omit(df2)
# save the dataset
saveRDS(df2, "Examples - data/data_ManyLabs_r.RDS")


##### Analysis script ####
df1 <- readRDS("Examples - data/data_ManyLabs_t.RDS")
df2 <- readRDS("Examples - data/data_ManyLabs_r.RDS")

# run RoBMA on each dataset or just load the saved results
load_RoBMA <- TRUE
res        <- list()
library(RoBMA)
for(i in 1:length(unique(df1$name))){
  
  if(load_RoBMA){
    temp_f <- readRDS(file = file.path("Examples - data", paste0("MLt", unique(df1$name)[i])))    
  }else{
    temp_data   <- df1[df1$name == unique(df1$name)[i], ]
    temp_f      <- RoBMA(t = temp_data$t, n = temp_data$n, parallel = TRUE, chains = 2, iter = 5000)
    saveRDS(temp_f, file = file.path("Examples - data", paste0("MLt", unique(df1$name)[i])))
  }
  
  res[[unique(df1$name)[i]]] <- temp_f$RoBMA$BF
}
for(i in 1:length(unique(df2$name))){
  
  if(load_RoBMA){
    temp_f <- readRDS(file = file.path("Examples - data", paste0("MLr", unique(df2$name)[i])))    
  }else{
    temp_data   <- df2[df2$name == unique(df2$name)[i], ]
    temp_f      <- RoBMA(r = temp_data$r, n = temp_data$n, parallel = TRUE, chains = 2, iter = 5000)
    saveRDS(temp_f, file = file.path("Examples - data", paste0("MLr", unique(df1$name)[i])))
  }
  
  res[[unique(df2$name)[i]]] <- temp_f$RoBMA$BF
}

res <- data.frame(do.call(rbind, res))

mean(res$bias < 1/3);sum(res$bias < 1/3) # evidence against bias
mean(res$bias > 3);sum(res$bias > 3)     # evidence for bias


# use test of the different methods
source("Examples - scripts/functions_other_methods.R")


res_others  <- NULL
for(i in 1:length(unique(df1$name))){
  
  temp_data   <- df1[df1$name == unique(df1$name)[i], ]
  temp_data$d <- psych::t2d(temp_data$t, n = temp_data$n)
  temp_data$v <- temp_data$n/(temp_data$n/2)^2 + temp_data$d^2/(2*temp_data$n)
  temp_data$df<- temp_data$n - 2


  temp_res <- suppressWarnings(rbind(
    RMA.est(d = temp_data$d, v = temp_data$v, long=TRUE),
    PETPEESE.est(temp_data$d, temp_data$v, PP.test = "one-sided", long=TRUE, runRMA=FALSE),
    pc_skew(t = temp_data$t, df = temp_data$df, long=TRUE),
    pcurveEst(t = temp_data$t, df = temp_data$df, progress=FALSE, long=TRUE, CI=FALSE),
    SM_twostep.est(d = temp_data$d, v = temp_data$v, p = pt(temp_data$t, temp_data$df, lower.tail = FALSE)),
    SM_threestep.est(d = temp_data$d, v = temp_data$v, p = pt(temp_data$t, temp_data$df, lower.tail = FALSE)),
    TES(t = temp_data$t, df = temp_data$df, N = temp_data$n, p = pt(abs(temp_data$t), temp_data$df, lower.tail = F)*2),
    TIVA(t = temp_data$t, df = temp_data$df),
    Egger(d = temp_data$d, v = temp_data$v),
    IC(t = temp_data$t, p = pt(abs(temp_data$t), temp_data$df, lower.tail = F)*2, N = temp_data$n)
  ))
  
  temp_res$name <- unique(df1$name)[i][1]
  
  res_others <- rbind(res_others, temp_res)
}
for(i in 1:length(unique(df2$name))){
  
  temp_data   <- df2[df2$name == unique(df2$name)[i], ]
  temp_data$d <- psych::r2d(temp_data$r)
  temp_data$v <- temp_data$n/(temp_data$n/2)^2 + temp_data$d^2/(2*temp_data$n)
  temp_data$df<- temp_data$n - 2
  temp_data$t <- psych::d2t(temp_data$d, temp_data$n - 2)
  
  temp_res <- suppressWarnings(rbind(
    RMA.est(d = temp_data$d, v = temp_data$v, long=TRUE),
    PETPEESE.est(temp_data$d, temp_data$v, PP.test = "one-sided", long=TRUE, runRMA=FALSE),
    pc_skew(t = temp_data$t, df = temp_data$df, long=TRUE),
    pcurveEst(t = temp_data$t, df = temp_data$df, progress=FALSE, long=TRUE, CI=FALSE),
    SM_twostep.est(d = temp_data$d, v = temp_data$v, p = pt(temp_data$t, temp_data$df, lower.tail = FALSE)),
    SM_threestep.est(d = temp_data$d, v = temp_data$v, p = pt(temp_data$t, temp_data$df, lower.tail = FALSE)),
    TES(t = temp_data$t, df = temp_data$df, N = temp_data$n, p = pt(abs(temp_data$t), temp_data$df, lower.tail = F)*2),
    TIVA(t = temp_data$t, df = temp_data$df),
    Egger(d = temp_data$d, v = temp_data$v),
    IC(t = temp_data$t, p = pt(abs(temp_data$t), temp_data$df, lower.tail = F)*2, N = temp_data$n)
  ))
  
  temp_res$name <- unique(df2$name)[i][1]
  
  res_others <- rbind(res_others, temp_res)
}




# keep only the bias tests
res_others <- process_results(res_others)
res_others <- res_others[res_others$term == "bias",]
res_others <- res_others[res_others$variable == "p.value",]
res_others$value <- as.numeric(res_others$value)
round(by(res_others$value < .10, res_others$method, mean), 2)
round(by(res_others$value < .10, res_others$method, sum), 2)



### create plot
temp_test <- rbind(
  data.frame(
    test   = ifelse(unlist(res$bias) < 1/3, -1, 
                    ifelse(unlist(res$bias) > 3, 1, 0)),
    method = "RoBMA"
  ),
  data.frame(
    test   = ifelse(as.numeric(res_others$value[res_others$term == "b1" & res_others$variable == "p.value"]) < .10, 1, 0),
    method = res_others$method[res_others$term == "b1" & res_others$variable == "p.value"]
  ),
  data.frame(
    test   = ifelse(as.numeric(res_others$value[res_others$term == "bias" & res_others$variable == "p.value"]) < .10, 1, 0),
    method = res_others$method[res_others$term == "bias" & res_others$variable == "p.value"]
  ),
  data.frame(
    test   = c(-1,0,1),
    method = "zzz"
  )
  
)
temp_test$method <- ifelse(temp_test$method == "TFtest", "TF", temp_test$method)

pdf("Examples - figures/ManyLabs2.pdf", width = 6, height = 4.5)
par(mar = c(3,7,1,2))
temp_test$method <- factor(temp_test$method, levels = c("TES", "Egger", "RI", "TF", "SM2", "SM3", "RoBMA", "zzz"))
temp_table <- table(temp_test$test*-1, temp_test$method)
temp_table <- temp_table[,-ncol(temp_table)]
for(i in 1:ncol(temp_table))temp_table[,i] <- temp_table[,i]/sum(temp_table[,i])
barplot(temp_table, horiz = TRUE, las = 1, xaxt = "n",
        xlab = "", ylab = "", main = "", cex.lab = 2, cex.axis = 2.5, cex.main = 2.5, cex.sub = 2, cex = 2)
axis(1, c(0, .25, .5, .75, 1), c("0", ".25", ".50", ".75", "1"), cex.lab = 2, cex.axis = 2.5, cex.main = 2.5, cex.sub = 2, cex = 2)
abline(v = .10, lwd = 2)
dev.off()

rownames(temp_table) <- c("bias", "no evidence", "no bias")
round(temp_table,3)*100


