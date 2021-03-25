#Load packages for CART
library(rpart) 
library(rpart.plot) 

#Load packages for RF
library(randomForest) 

#Load packages for stats
library(rattle) 
library(Metrics)

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
corM <- matrix(nrow=11,ncol=2)
colnames(corM)<-c("r","p-value")
rownames(corM)<-c(names(data))

for(i in 1:11){
  cor <- cor.test(data[,1],data[,i])
  corM[i,1]<-cor$estimate
  corM[i,2]<-cor$p.value
}

#-------------T tests-------------
tM <- matrix(nrow=3,ncol=1)
rownames(tM)<-c("CART pH v.s. Real pH","RF pH v.s. Real pH", "CART pH v.s. RF pH")
colnames(tM)<-c("p-value")

#CART pH v.s. actual pH
CARTpred<-t.test(CARTpredictions,test_setActual,alternative="two.sided")
tM[1,1]<-CARTpred$p.value
#RF pH v.s. actual pH
RFpred<-t.test(RFpredictions,test_setActual,alternative="two.sided")
tM[2,1]<-RFpred$p.value
#CART pH v.s. RF pH
comparison<-t.test(CARTpredictions,RFpredictions,alternative="two.sided")
tM[3,1]<-comparison$p.value

#-------------Print outputs-------------
corM
tM
CARTimp
RFimp
CARTRRMSEmatrix
RFRRMSEmatrix


