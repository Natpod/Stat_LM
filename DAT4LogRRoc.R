# PRACTICAL ASSIGNMENT - LIQUID BIOPSY DATA (DAT4)
# Author: Natalia García Sánchez
# Date : 15/12/2022
# Description: Evaluation of variables in liquid biopsy
# dataset for pancreatic cancer prediction.
# ----------------------------------------------------

#----- Library imports
library("skimr")
library('MASS')
library('dplyr')
library('ROCR')
library('caret')
library('pROC')

#----- Data exploration

# After loading workspace from .Rdata file with DAT4 table,
# a summary of the data will be provided

summary(DAT4)
skim(DAT4)

# write into csv for report

write.csv(summary(DAT4), "sum_DAT4.csv", row.names=FALSE)



#----- Data preprocessing

# process clinical groups as binary qualitative type (terminal-M or non terminal)
DAT4$group <- ifelse(DAT4$group == "M", "terminal", "nonterminal")

#----- 2 by 2 Contingency matrixes 


# (sex - status)
sexcont <- as.matrix(table(DAT4[,c(2, 4)]))
chisq.test(sexcont)$expected
chisq.test(sexcont)

# (group-status)
groupcont <- table(DAT4[,c(7, 4)])
chisq.test(groupcont)$expected
chisq.test(groupcont)

# (CTC - status)
ctccont <- table(DAT4[,c(5, 4)])
chisq.test(ctccont)$expected
chisq.test(ctccont)

# (Mutation - status)
mutcont <- table(DAT4[,c(6, 4)])
chisq.test(mutcont)$expected
chisq.test(mutcont)

# we deem the quantitative variables not representative 
# enough as risk factors for survival

#----------------------------------------------------


# association btw two variables

ICORRR = function(T, alpha = 0.05, deci = 2) {
  a = T[1, 1]
  b = T[1, 2]
  c = T[2, 1]
  d = T[2, 2]
  OR = (a * d)/(b * c)
  RR = (a/(a + b))/(c/(c + d))
  sdOR = sqrt((1/a) + (1/b) + (1/c) + (1/d))
  sdRR = sqrt((1/a) - (1/(a + b)) + (1/b) -
                 (1/(c + d)))
  upperOR = OR * exp(qnorm(1 - alpha/2) *
                      sdOR)
  lowerOR = OR * exp(-qnorm(1 - alpha/2) *
                      sdOR)
  upperRR = RR * exp(qnorm(1 - alpha/2) *
                      sdRR)
  lowerRR = RR * exp(-qnorm(1 - alpha/2) *
                      sdRR)
  sol = matrix(0, ncol = 2, nrow = 2)
  sol[1, 1] = round(lowerOR, deci)
  sol[1, 2] = round(upperOR, deci)
  sol[2, 1] = round(lowerRR, deci)
  sol[2, 2] = round(upperRR, deci)
  colnames(sol) = c("lower", "upper")
  rownames(sol) = c("OR", "RR")
  return(sol)
}

ICORRR(groupcont)
# Simulamos los datos para hacer un pearson
set.seed(1234)

evital <- sample(c("event", "censored"), replace = TRUE,
                  300)
G <- matrix(0, nrow = 300, ncol = 100)
for (i in 1:100) {
  k <- runif(1, min = 0.1, max = 0.9)
  G[, i] <- sample(c("terminal", "nonterminal"), replace = TRUE,
                    size = 300, prob = c(k, 1 - k))
}
sol <- matrix(0, nrow = 100, ncol = 4)
for (i in 1:100) {
  tf <- fisher.test(table(evital, G[, i]))
  sol[i, ] <- c(tf$estimate, tf$conf.int,
                tf$p.value)
}
sol <- round(sol, 3)
colnames(sol) <- c("OR", "lw", "up", "p.value")
rownames(sol) <- paste("G", 1:100, sep = "")
# instancias de datos con pval menor que sg factor
sol[sol[, 4] < 0.05, ]


ICORRR(mutcont)
# Simulamos los datos para hacer un pearson
set.seed(1234)

evital <- sample(c("event", "censored"), replace = TRUE,
                  300)
G <- matrix(0, nrow = 300, ncol = 100)
for (i in 1:100) {
  k <- runif(1, min = 0.1, max = 0.9)
  G[, i] <- sample(c("NEG", "POS"), replace = TRUE,
                    size = 300, prob = c(k, 1 - k))
}
sol <- matrix(0, nrow = 100, ncol = 4)
for (i in 1:100) {
  tf <- fisher.test(table(evital, G[, i]))
  sol[i, ] <- c(tf$estimate, tf$conf.int,
                tf$p.value)
}
sol <- round(sol, 3)
colnames(sol) <- c("OR", "lw", "up", "p.value")
rownames(sol) <- paste("G", 1:100, sep = "")
# instancias de datos con pval menor que sg factor
sol[sol[, 4] < 0.05, ]


# --------- Logistic regression


# modelo lineal var predictor/prediccion

DAT4$sexo = as.factor(DAT4$sex)
DAT4$group = as.factor(DAT4$group)
DAT4$CTC = as.factor(DAT4$CTC)
DAT4$Mutant  = as.factor(DAT4$Mutant)
DAT4$status = ifelse(DAT4$status == "event", 1, 0)

fisher.test(table(DAT4$status, DAT4$sexo))
mel1 = glm(status ~ sexo, data = DAT4,family = binomial)
summary(mel1)
#Transformamos las estimas en ORs
round(exp(cbind(OR = coef(mel1), confint(mel1))),3)


#To select the best model for a data set taking into account the AIC
#criterion, we have the function step(), which has as an argument a tted
#model with all possible variables.
mels = step(glm(status ~ sexo + group + CTC + Mutant + edad, data = DAT4, family = binomial))

# best model is the one that takes into account all variables
mel3 = glm(status ~ sexo + group + CTC + Mutant + edad, data = DAT4, family = binomial)
#Roc curve
predic = predict(mel3, as.data.frame(DAT4), type = "response")
# Melanoma$status01 the reality
pred = prediction(predic, DAT4$status)

# sensitivity and Specificity
perf = performance(pred, "tpr", "fpr")
perf2 = performance(pred, "auc") # Accuracy
auc = unlist(slot(perf2, "y.values"))

#roc area ci
CIROC95 = function(AUC, N1, N2, q = 0.95) {
  1.96=Z(alfa/2)
  z = qnorm(q)
  Q1 = AUC/(2 - AUC)
  Q2 = (2 * AUC^2)/(1 + AUC)
  A = AUC * (1 - AUC)
  B = (N1 - 1) * (Q1 - AUC^2)
  C = (N2 - 1) * (Q2 - AUC^2)
  se = sqrt((A * B * C)/(N1 * N2))
  l = AUC - z * se
  u = AUC + z * se
  return(c(l, u))
  }
x11()
plot(perf)
lines(c(0, 1))

text(auc)

#------------------ 3 fold cross validation
library(caret)
DAT4$status<-as.factor(DAT4$status)
#specify the cross-validation method
ctrl <- trainControl(method = "cv", number = 3)

#fit a regression model and use k-fold CV to evaluate performance
model <- train(status ~  sexo + group + CTC + Mutant + edad, data = DAT4, method = "glm", trControl = ctrl)

#view summary of k-fold CV               
print(model)
predict.model

# Select a parameter setting
selectedIndices <- model$pred$mtry == 2
# Plot:
plot.roc(model$pred$obs[selectedIndices],
         model$pred$M[selectedIndices])





# ----------------
# Partition data and create index matrix of selected values
index <- createDataPartition(DAT4$status, p=.8, list=FALSE, times=1)
# Create test and train data frames
train_df <- DAT4[index,]
test_df <- DAT4[-index,]
nrow(train_df)
nrow(test_df)

# Re-label values of outcome variable for train_df
train_df$status[train_df$status==1] <- "event"
train_df$status[train_df$status==0] <- "censored"

# Re-label values of outcome variable for test_df
test_df$status[test_df$status==1] <- "event"
test_df$status[test_df$status==0] <- "censored"

# Convert outcome variable to factor for each data frame
train_df$status <- as.factor(train_df$status)
test_df$status <- as.factor(test_df$status)

# Specify type of training method used and the number of folds
ctrlspecs <- trainControl(method="cv", 
                          number=3, 
                          savePredictions="all",
                          classProbs=TRUE)
# Set random seed for subsequent random selection and assignment operations
set.seed(1985)

# Specify logistic regression model to be estimated using training data
# and k-fold cross-validation process
model1 <- train(status ~ sexo + group + CTC + Mutant + edad, data=train_df, 
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)

# Print information about model
print(model1)

# Print results of final model estimated using training data
summary(model1)

# Estimate the importance of different predictors
varImp(model1)

# Predict outcome using model from training data based on testing data
predictions <- predict(model1, newdata=test_df)

# Create confusion matrix to assess model fit/performance on test data
confusionMatrix(data=predictions, test_df$status)

# Predict outcome using model from training data based on testing data
predictions <- predict(model1, newdata=train_df)

# Create confusion matrix to assess model fit/performance on test data
confusionMatrix(data=predictions, train_df$status)

#roc

#TESTING WITH LOGM test
class <- predict(model1, test_df)
probs <- predict(model1, test_df,'prob')

TEST.scored <- cbind(test_df,class,probs) %>% mutate(data = "TEST")

# TESTING TRAIN
class <- predict(model1, train_df)
probs <- predict(model1, train_df,'prob')

TRAIN.scored <- cbind(train_df,class,probs) %>% mutate(data = "TRAIN")

TRAIN_TEST_scored <- rbind(TRAIN.scored, TEST.scored)
TRAIN_TEST_scored$truelabel<-c(DAT4[index,4], DAT4[-index,4])
  
library('yardstick')
TRAIN_TEST_scored %>%
       group_by(data) %>%
       roc_curve(truth=truelabel,"censored") %>%
       autoplot()

TRAIN_TEST_scored %>%
  group_by(data) %>%
  roc_auc(truth=truelabel,"censored")











# Specify logistic regression model to be estimated using training data
# and k-fold cross-validation process
model2 <- train(status ~ group + CTC + Mutant , data=train_df, 
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)

# Print information about model
print(model2)

# Print results of final model estimated using training data
summary(model2)

# Estimate the importance of different predictors
varImp(model2)

# Predict outcome using model from training data based on testing data
predictions <- predict(model2, newdata=test_df)

# Create confusion matrix to assess model fit/performance on test data
confusionMatrix(data=predictions, test_df$status)

# Predict outcome using model from training data based on testing data
predictions <- predict(model2, newdata=train_df)

# Create confusion matrix to assess model fit/performance on test data
confusionMatrix(data=predictions, train_df$status)

#roc

#TESTING WITH LOGM test
class <- predict(model2, test_df)
probs <- predict(model2, test_df,'prob')

TEST.scored <- cbind(test_df,class,probs) %>% mutate(data = "TEST")

# TESTING TRAIN
class <- predict(model2, train_df)
probs <- predict(model2, train_df,'prob')

TRAIN.scored <- cbind(train_df,class,probs) %>% mutate(data = "TRAIN")

TRAIN_TEST_scored <- rbind(TRAIN.scored, TEST.scored)
TRAIN_TEST_scored$truelabel<-c(DAT4[index,4], DAT4[-index,4])

library('yardstick')
TRAIN_TEST_scored %>%
  group_by(data) %>%
  roc_curve(truth=truelabel,"censored") %>%
  autoplot()

TRAIN_TEST_scored %>%
  group_by(data) %>%
  roc_auc(truth=truelabel,"censored")

