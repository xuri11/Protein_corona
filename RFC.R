library(readxl)
library(dplyr)
library(corrplot)
library(randomForest)
library(pROC)
library(rpart)
dat <- read_excel('RandomForest_database.xlsx', sheet = 'RFC')   ##Supporting Data 3
dat1 <- mutate(dat, H_bonding = (dat$P_Ser + dat$P_Thr + dat$P_Asn + dat$P_Gln),
               P_aromatic = (dat$P_Phe + dat$P_Trp + dat$P_Tyr))%>%
  select("H_bonding", "P_Cys", "P_aromatic",
         "PI", "MW", "GRAVY", "Instability_index",
         "Aliphatic_indices")
cor_heat1 <- cor(dat1)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
corrplot(cor_heat1, method = 'pie',
         order = 'hclust',
         type = 'upper',
         tl.col = 'black',
         tl.cex = 1,
         tl.srt = 45,
         col = rev(col2(200)),
         diag = F,
         addCoef.col = 'black')
dat_RFC <- mutate(dat, H_bonding = (dat$P_Ser + dat$P_Thr + dat$P_Asn + dat$P_Gln),
                  P_aromatic = (dat$P_Phe + dat$P_Trp + dat$P_Tyr), 
                  Fold = log(dat$Enrichment))
colnames(dat_RFC)
dat_RFC <- select(dat_RFC, "H_bonding", "P_Cys", "P_aromatic",
                  "PI", "MW", "GRAVY", "Instability_index",
                  "Particle_Size",  "Particle_Charge", "Fold" )
dat_RFC <- dat_RFC[dat_RFC$Fold != 0,]
dat_RFC$Fold[dat_RFC$Fold > 0] <- 1
dat_RFC$Fold[dat_RFC$Fold < 0] <- 0
table(dat_RFC$Fold)
dat_RFC$Fold <- as.factor(dat_RFC$Fold)
dat_RFC$Particle_Charge <- as.factor(dat_RFC$Particle_Charge)
dat_RFC$Particle_Size <- as.factor(dat_RFC$Particle_Size)
str(dat_RFC)
df <- as.data.frame(dat_RFC)%>%
  na.omit()
Fold = function(Z, w, D, seed){
  n = nrow(w)
  d = 1:n
  e = levels(w[,D])
  N = length(e)
  set.seed(seed)
  dd = lapply(1:N, function(i){
    d0 = d[w[,D] == e[i]]
    j = length(d0)
    ZT = rep(1:Z, ceiling(j/Z))[1:j]
    id = cbind(sample(ZT), d0);id})
  mm = lapply(1:Z, function(i){
    u = NULL; for (j in 1:N)
      u = c(u, dd[[j]][dd[[j]][,1] == i,2]);u
  })
  return(mm)}
mm <- Fold(Z = 10, w = df, D = 10, seed = 666)
mm
Z = 10; w = df
E0 = rep(0,Z); E1 = E0
n = length(w$Fold)
for (i in 1:Z) {
  m = mm[[i]]
  n0 = n-length(m); n1 = length(m)
  a = randomForest(Fold~., w[-m,])
  E0[i] = sum(w[-m,10]!=predict(a, w[-m,], type = 'class'))/n0
  E1[i] = sum(w[m,10]!=predict(a, w[m,], type = 'class'))/n1}
mean(E0); mean(E1)

E2 = rep(0,Z); E3 = E2
for (i in 1:Z) {
  m = mm[[i]]
  n0 = n-length(m); n1 = length(m)
  a = rpart(Fold~., w[-m,])
  E2[i] = sum(w[-m,10]!=predict(a, w[-m,], type = 'class'))/n0
  E3[i] = sum(w[m,10]!=predict(a, w[m,], type = 'class'))/n1}
mean(E2); mean(E3)
err1 <- data.frame(E0, E2)
names(err1) <- c('RFC', 'DTC')
err2 <- data.frame(E1, E3)
names(err2) <- c('RFC', 'DTC')
boxplot(err1, ylab = 'Error rate', col = 'lightblue', main = 'Train data')
boxplot(err2, ylab = 'Error rate', col = 'lightblue', main = 'Test data')

train <- w[-mm[[which.min(E1)]],]
test <- w[mm[[which.min(E1)]],]
set.seed(666)
rdf = randomForest(Fold~., train, ntree = 1000)
print(rdf)
plot(rdf)
t <- tuneRF(train[, -10], train[, 10],
            stepFactor = 0.5,
            plot = T,
            ntreeTry = 500,
            trace = T,
            improve = 0.05)
set.seed(123)
rdf = randomForest(Fold~., train, ntree = 500,
                   mtry = 3, importanc = T,
                   proximity = T)
print(rdf)
hist(treesize(rdf), col = 'lightblue', main = 'NO. of the Nodes for the Trees')
varImpPlot(rdf, main = 'Variable importance', type = 1)
pred <- predict(rdf, newdata = test)
roc <- roc(test$Fold, as.numeric(pred))
plot(roc, print.auc=TRUE, auc.polygon=TRUE, 
     grid=c(0.1, 0.2),grid.col=c("green", "red"),
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
     print.thres=TRUE,main='ROC Curve for Random Forest',
     type = 'bars')
