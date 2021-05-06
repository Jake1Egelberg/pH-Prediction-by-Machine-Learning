#Load packages for CART
library(rpart) 
library(rpart.plot) 

#Load packages for RF
library(randomForest) 

#Load packages for stats
library(rattle) 
library(Metrics)
library(car)

#Set working directory to hard drive
setwd("D:\\pHdata\\Outputs")

#Set Seed 
#initial_seed=as.integer(Sys.time()) <- used to generate seed, run once 
#the_seed=initial_seed %% 100000 
the_seed <- 27137 
set.seed(the_seed) 

#Set data source
data <- read.csv("D:\\pHdata\\SamplingLocationData.csv")

#Split data into 70% testing and 30% training sets 
#train takes random sample of number of rows from data 
train<-sample(1:nrow(data),size=ceiling(0.7*nrow(data)),replace=FALSE)  
#sets training and test sets to matrices containing only rows sampled by train  
training_set<-data[train,] #matrix of 70% rows w/ columns 
test_set<-data[-train,] #matrix of 30% rows w/ columns 

#-------------USER INPUTS-------------

dependent.variable <- as.name('ï..pH')
CARTmethod <- "anova" 
RFmethod <- "anova" 

#-------------END USER INPUTS-------------

#-------------MODEL GENERATION-------------

#Build CART model 
CARTmodel <- rpart(ï..pH ~., data=training_set, method=CARTmethod, xval=10, minsplit=20) 

#Prune CART model 
CARTmodel<- prune(CARTmodel, cp=CARTmodel$cptable[which.min(CARTmodel$cptable[,"xerror"]),"CP"]) 

#Set RF max nodes
observations2 <- nrow(training_set)
max.nodes <- round(observations2*0.3,0)

#Build RF 
RFmodel <- randomForest(ï..pH ~.,data=training_set, ntree=1000,importance=TRUE,nodesize=1,maxnodes=max.nodes) #INPUT


#Visualize CART Model as .tiff with high res 
tiff("CARToutput.tiff",width=4,height=4,units="in",res=300)
plot1<-fancyRpartPlot(CARTmodel)
dev.off()

#Visualize MSE as trees added
tiff("RFerror.tiff",width=4,height=4,units="in",res=300)
plot1<-plot(RFmodel)
dev.off()

#-------------MODEL PREDICTIONS-------------

#Gets column num for pH in test dataset
dependentIndex <- grep(dependent.variable,colnames(test_set))
#Creates dataset to make predictions with (missing pH data)
predict_set <- test_set[-dependentIndex]
#Creates list of actual pH values in test set
test_setActual <- test_set[as.character(dependent.variable)]
#converts list to numeric vector
test_setActual <- unlist(test_setActual) 

#CART predict testing set
if(CARTmethod=='anova'){
  type.value <- 'vector' #yields vector of responses
}

#Store CART predictions of test set pH values
CARTpredictions <- predict(CARTmodel,predict_set,type=type.value) 

#CART RRMSE 
if(CARTmethod=="anova"){
  CARTpredictions <- as.numeric(CARTpredictions) #ensures vector elements are numeric 
  CARTrmse <- rmse(test_setActual,CARTpredictions)
  meantest_setActual <- mean(test_setActual)
  CARTrrmse <- CARTrmse/meantest_setActual
  
  CARTRRMSEmatrix <- matrix(nrow=1,ncol=1)
  colnames(CARTRRMSEmatrix)<-c("% RRMSE")
  rownames(CARTRRMSEmatrix)<-c("CART Model")
  CARTRRMSEmatrix[,1]<-c(CARTrrmse*100)
} 

#Store CART variable importance
CARTimp<-CARTmodel$variable.importance

#Store RF predictions of test set pH values 
RFpredictions <- predict(RFmodel,predict_set) 

#RF RRMSE 
if(RFmethod=="anova"){
  RFpredictions <- as.numeric(RFpredictions) #ensures vector elements are numeric 
  RFrmse <- rmse(test_setActual,RFpredictions)
  meantest_setActual <- mean(test_setActual)
  RFrrmse <- RFrmse/meantest_setActual
  
  RFRRMSEmatrix <- matrix(nrow=1,ncol=1)
  colnames(RFRRMSEmatrix)<-c("% RRMSE")
  rownames(RFRRMSEmatrix)<-c("RF Model")
  RFRRMSEmatrix[,1]<-c(RFrrmse*100)
} 

#Store RF model importance
RFimp<-RFmodel$importance

#-------------Correlations to pH-------------
corM <- matrix(nrow=11,ncol=3)
colnames(corM)<-c("rho","p-value","Adjusted p-value")
rownames(corM)<-c(names(data))

for(i in 1:11){
  cor <- cor.test(data[,1],data[,i],method="spearman")
  adjustedP<-p.adjust(cor$p.value,method="bonferroni",n=10) #10 associations b/w variables and pH
  
  corM[i,1]<-round(cor$estimate,digits=4)
  corM[i,2]<-round(cor$p.value,digits=4)
  corM[i,3]<-round(adjustedP,digits=4)
}

#Generating Zhang pH predictions
ZhangpH<-vector(length=30)
ZhangpH[1:30]<-5.43

#-------------Stat tests-------------
tM <- matrix(nrow=4,ncol=2)
rownames(tM)<-c("CART pH v.s. Real pH","RF pH v.s. Real pH", "CART pH v.s. RF pH","Zhang pH v.s. Real pH")
colnames(tM)<-c("Levene's Test p-value","Welch's t-test p-value")

#CART pH v.s. actual pH
CARTpred<-t.test(CARTpredictions,test_setActual,alternative="two.sided",var.equal = FALSE)
tM[1,2]<-round(CARTpred$p.value,digits=4)

CARTpredFactor<-as.factor(CARTpredictions)
CARTlev<-leveneTest(test_setActual~CARTpredFactor)
tM[1,1]<-round(CARTlev$`Pr(>F)`[1],digits=4)

#RF pH v.s. actual pH
RFpred<-t.test(RFpredictions,test_setActual,alternative="two.sided",var.equal = FALSE)
tM[2,2]<-round(RFpred$p.value,digits=4)

RFpredFactor<-as.factor(RFpredictions)
RFlev<-leveneTest(test_setActual~RFpredFactor)
tM[2,1]<-round(RFlev$`Pr(>F)`[1],digits=4)

#CART pH v.s. RF pH
comparison<-t.test(CARTpredictions,RFpredictions,alternative="two.sided",var.equal = FALSE)
tM[3,2]<-round(comparison$p.value,digits=4)

complev<-leveneTest(RFpredictions~CARTpredFactor)
tM[3,1]<-round(complev$`Pr(>F)`[1],digits=4)

#Zhang pH v.s. Real pH
Zhangpred<-t.test(ZhangpH,test_setActual,alternative="two.sided",var.equal = FALSE)
tM[4,2]<-round(Zhangpred$p.value,digits=4)

#ZhangPredFactor<-as.factor(ZhangpH)
#Zhanglev<-leveneTest(test_setActual~ZhangPredFactor)
#tM[4,1]<-round(Zhanglev$`Pr(>F)`[1],digits=4)


#-------------Print outputs-------------
corM
tM
CARTimp
RFimp
CARTRRMSEmatrix
RFRRMSEmatrix


