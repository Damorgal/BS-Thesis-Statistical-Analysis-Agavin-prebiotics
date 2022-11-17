# Code in R for the statistical analysis and modeling of the exploratory randomized clinical trial: "Silva-Adame, M.B.; Martínez-Alvarado, A.; Martínez-Silva, V.A.; Samaniego-Méndez, V.; López, M.G. Agavins Impact on Gastrointestinal Tolerability-Related Symptoms during a Five-Week Dose-Escalation Intervention in Lean and Obese Mexican Adults: Exploratory Randomized Clinical Trial. Foods 2022, 11, 670. https://doi.org/10.3390/foods11050670"

# By: Diego A. Moreno-Galvan

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---


# Cross-validation and prediction error -----------------------------------------------------------
library(caret)

list = 0
list_good = 0
set.seed(1614)
for (i in 1:1000) {
    folds <- createFolds(features$GRUPO, k = 5)
    error = 0
    good = 0

    # Block 1 prediction
    test_set <- features[folds$Fold1,]
    train_set <- rbind(features[folds$Fold2,],
                       features[folds$Fold3,],
                       features[folds$Fold4,],
                       features[folds$Fold5,])
    # Training
    model <- glm(CLASIFICACION ~ ., family=binomial(link=logit), data=train_set)

    # Prediction
    probabilities = predict(model, test_set, type="response")
    y_pred <- ifelse(probabilities > 0.5, 1, 0)
    good = good + mean(y_pred == test_set$CLASIFICACION)*100
    error = error + sum((y_pred-test_set$CLASIFICACION)^2)/length(folds$Fold1)

    # Block 2 prediction
    test_set <- features[folds$Fold2,]
    train_set <- rbind(features[folds$Fold1,],
                       features[folds$Fold3,],
                       features[folds$Fold4,],
                       features[folds$Fold5,])
    model <- glm(CLASIFICACION ~ ., family=binomial(link=logit), data=train_set)
    probabilities = predict(model, test_set, type="response")
    y_pred <- ifelse(probabilities > 0.5, 1, 0)
    good = good + mean(y_pred == test_set$CLASIFICACION)*100
    error = error + sum((y_pred-test_set$CLASIFICACION)^2)/length(folds$Fold2)

    # Block 3 prediction
    test_set <- features[folds$Fold3,]
    train_set <- rbind(features[folds$Fold1,],
                       features[folds$Fold2,],
                       features[folds$Fold4,],
                       features[folds$Fold5,])
    model <- glm(CLASIFICACION ~ ., family=binomial(link=logit), data=train_set)
    probabilities = predict(model, test_set, type="response")
    y_pred <- ifelse(probabilities > 0.5, 1, 0)
    good = good + mean(y_pred == test_set$CLASIFICACION)*100
    error = error + sum((y_pred-test_set$CLASIFICACION)^2)/length(folds$Fold3)

    # Block 4 prediction
    test_set <- features[folds$Fold4,]
    train_set <- rbind(features[folds$Fold1,],
                       features[folds$Fold2,],
                       features[folds$Fold3,],
                       features[folds$Fold5,])
    model <- glm(CLASIFICACION ~ ., family=binomial(link=logit), data=train_set)
    probabilities = predict(model, test_set, type="response")
    y_pred <- ifelse(probabilities > 0.5, 1, 0)
    good = good + mean(y_pred == test_set$CLASIFICACION)*100
    error = error + sum((y_pred-test_set$CLASIFICACION)^2)/length(folds$Fold4)

    # Block 5 prediction
    test_set <- features[folds$Fold5,]
    train_set <- rbind(features[folds$Fold1,],
                       features[folds$Fold2,],
                       features[folds$Fold3,],
                       features[folds$Fold4,])
    model <- glm(CLASIFICACION ~ ., family=binomial(link=logit), data=train_set)
    probabilities = predict(model, test_set, type="response")
    y_pred <- ifelse(probabilities > 0.5, 1, 0)
    good = good + mean(y_pred == test_set$CLASIFICACION)*100
    error = error + sum((y_pred-test_set$CLASIFICACION)^2)/length(folds$Fold5)
    error = error/5
    good = good/5
    list = c(list, error)
    list_good = c(list_good, good)
}
mean(list[-1])
mean(list_good[-1])

# AUROC -------------------------------------------------------------------------------------------
probabilities = predict(fit, features, type="response")
sens = 0
spe = 0
x <- seq(0,1,by=0.01)
for(i in 1:length(x)) {
    y_pred <- ifelse(probabilities > x[i], 1, 0)
    sens = c(sens, sensitivity(factor(as.numeric(y_pred), levels=c(1,0)),
                               factor(features$CLASIFICACION, levels=c(1,0))))
    spe = c(spe, specificity(factor(as.numeric(y_pred), levels=c(1,0)),
                             factor(features$CLASIFICACION, levels=c(1,0))))
}
plot(1-spe[-1], sens[-1],
   pch=20, main='Curva ROC del modelo',
   xlab = '1-Specificity',
   ylab = 'Sensitivity', col='darkblue')
lines(1-spe[-1], sens[-1], col='blue')
abline(a=c(0,0),b=c(1,1),col='red',lty='dotted')
# AUC-ROC
aucroc <- 0.7857*0.1666+(0.2916-0.1666)*0.9285+(1-.2916)
aucroc


# Cross-validation prediction ---------------------------------------------------------------------
predicciones = 0
set.seed(1614)
rdata2 <- rdata
for (i in 1:3000) {
    # Lean-Agavins
    idx <- rdata2$GRUPO == 1 & rdata2$TRATAMIENTO == 1
    n_a <- rdata2[idx,]
    split <- sample.split(n_a$GRUPO, SplitRatio = 0.90)
    train.data  <- subset(n_a, split == TRUE)
    test.data <- subset(n_a, split == FALSE)
    # Lean-Placebo
    idx <- rdata2$GRUPO == 1 & rdata2$TRATAMIENTO == 0
    n_p <- rdata2[idx,]
    split <- sample.split(n_p$GRUPO, SplitRatio = 0.90)
    train.data  <- rbind(train.data, subset(n_p, split == TRUE))
    test.data <- rbind(test.data, subset(n_p, split == FALSE))
    # Obese-Agavins
    idx <- rdata2$GRUPO == 0 & rdata2$TRATAMIENTO == 1
    o_a <- rdata2[idx,]
    split <- sample.split(o_a$GRUPO, SplitRatio = 0.90)
    train.data  <- rbind(train.data, subset(o_a, split == TRUE))
    test.data <- rbind(test.data, subset(o_a, split == FALSE))
    # Obese-Placebo
    idx <- rdata2$GRUPO == 0 & rdata2$TRATAMIENTO == 0
    o_p <- rdata2[idx,]
    split <- sample.split(o_p$GRUPO, SplitRatio = 0.90)
    train.data  <- rbind(train.data, subset(o_p, split == TRUE))
    test.data <- rbind(test.data, subset(o_p, split == FALSE))

    # Training
    x <- model.matrix(CLASIFICACION ~ ., train.data)[,-1]
    y <- train.data$CLASIFICACION
    # Best lambda
    fit <- glmnet(x, y, alpha = 1, family = "binomial",
                  #lambda = cv.lasso$lambda.min)
                  #lambda = cv.lasso$lambda.1se)
                  lambda = exp(-4.7073854))
                  #lambda = 0.02342671)

    # Pred
    x.test <- model.matrix(CLASIFICACION ~ ., test.data)[,-1]
    probabilities <- predict(fit, newx = x.test, type="response")
    y_pred <- ifelse(probabilities > 0.5, 1, 0)
    predicciones <- c(predicciones, mean(y_pred == test.data$CLASIFICACION)*100)
}
predicciones = predicciones[-1]
hist(predicciones, main=c("Histograma de 3000 predicciones","modelo mínima desviación"),
   xlab="Porcentaje de Buena Clasificación", ylab="Frecuencia",
   col = "#33CC33", breaks = seq(0,100,10))



# Odds plot #######################################################################################
aux <- features[3:(dim(features)[2]-1)]
suma <- rep(0,dim(aux)[1])
for (i in 1:dim(aux)[1]) {
    suma[i] <- sum(aux[i,]*coefi[4:(length(coefi))])
}
inter = coefi[1]
trat = coefi[3]
grup = coefi[2]

pacientes = 0
for (i in 1:dim(aux)[1]) {
    if (features$GRUPO[i] == 1) {
          if (features$TRATAMIENTO[i] == 1) {
              pacientes <- c(pacientes, inter+trat+grup+suma[i])
          } else {
              pacientes <- c(pacientes, inter+grup+suma[i])
          }
    }
    else {
          if (features$TRATAMIENTO[i] == 1) {
              pacientes <- c(pacientes, inter+trat+suma[i])
          } else {
              pacientes <- c(pacientes, inter+suma[i])
          }
    }
}
pacientes = pacientes[-1]
odds <- log(probabilities)-log(1-probabilities)

# Graph -----------------------------------------
graph <- data.frame(pacientes, odds)

ggplot(graph, aes(x = pacientes, y = odds))  +
geom_abline(slope=1, color='orange')+
geom_point(color="blue", size=1)+
#scale_x_discrete(limits=c("0", "2.5","5","7", "10","12"))+
#scale_x_continuous(limits=c(-15,15)) +
#scale_y_continuous(limits=c(-15,15)) +
labs(x = expression(paste("Combinación lineal de pacientes ",p[i])),
     y = expression(paste(log,
                          bgroup("(",frac(paste("P(Y=1","|","X=",p[i],")"),
                                          1-paste("P(Y=1","|","X=",p[i],")")),")"))),
     title = "Ajuste del modelo de regresión logística")





# Logistic regression with LASSO regularization ###################################################
# Best lambda using cross-validation

set.seed(123)
x <- model.matrix(CLASIFICACION ~ ., features)[,-1]
y <- features$CLASIFICACION
cv.lasso <- cv.glmnet(x, y, alpha=1, family="binomial")
plot(cv.lasso)
cv.lasso$lambda.1se = exp(-4.7073854)
fit <- glmnet(x, y, alpha = 1, family = "binomial",
            lambda = cv.lasso$lambda.min)
            #lambda = cv.lasso$lambda.1se)
            #lambda = exp(-4.7073854))

# Coeficients
coefi <- coef(fit)
coefi
dif_residuos <- fit$nulldev - deviance(fit)
df = fit$df
p_value <- pchisq(q = dif_residuos,df = df, lower.tail = FALSE)
paste("Diferencia de residuos:", round(dif_residuos, 4))
paste("Grados de libertad:", df)
paste("p-value:", p_value)

# Predictions
test.data = features
x.test <- model.matrix(CLASIFICACION ~ ., test.data)[,-1]
probabilities <- predict(fit, newx = x.test, type="response")
y_pred <- ifelse(probabilities > 0.5, 1, 0)
mean(y_pred == test.data$CLASIFICACION)*100

# Agavins dose efect visualization ----------------------------------------------------------------
aux <- features[3:(dim(features)[2]-1)]
suma <- rep(0,dim(aux)[1])
for (i in 1:dim(aux)[1]) {
    suma[i] <- sum(aux[i,]*coefi[4:(length(coefi))])
}
inter = coefi[1]
trat = coefi[3]
grup = coefi[2]

N_A = N_P = O_A = O_P = 0
for (i in 1:dim(aux)[1]) {
    if (features$GRUPO[i] == 1) {
          if (features$TRATAMIENTO[i] == 1) {
            N_A <- c(N_A, 1/(1+exp(-(inter+trat+grup+suma[i]))))
          } else {
            N_P <- c(N_P, 1/(1+exp(-(inter+grup+suma[i]))))
          }
    }
    else {
          if (features$TRATAMIENTO[i] == 1) {
            O_A <- c(O_A, 1/(1+exp(-(inter+trat+suma[i]))))
          } else {
            O_P <- c(O_P, 1/(1+exp(-(inter+suma[i]))))
          }
    }
}
# Lean-Agavins
idx <- features$GRUPO==1 & features$TRATAMIENTO==1 & features$GENERO == 1
NAF = probabilities[idx]
idx <- features$GRUPO==1 & features$TRATAMIENTO==1 & features$GENERO == 0
NAH = probabilities[idx]
# Lean-Placebo
idx <- features$GRUPO==1 & features$TRATAMIENTO==0 & features$GENERO == 1
NPF = probabilities[idx]
idx <- features$GRUPO==1 & features$TRATAMIENTO==0 & features$GENERO == 0
NPH = probabilities[idx]
# Obese-Agavins
idx <- features$GRUPO==0 & features$TRATAMIENTO==1 & features$GENERO == 1
OAF = probabilities[idx]
idx <- features$GRUPO==0 & features$TRATAMIENTO==1 & features$GENERO == 0
OAH = probabilities[idx]
# Obese-Placebo
idx <- features$GRUPO==0 & features$TRATAMIENTO==0 & features$GENERO == 1
OPF = probabilities[idx]
idx <- features$GRUPO==0 & features$TRATAMIENTO==0 & features$GENERO == 0
OPH = probabilities[idx]

NAF = data.frame("Normopeso Mujer", "Agavinas", NAF)#[-1])
names(NAF) = c("GRUPO","TRATAMIENTO","fit")
NAH = data.frame("Normopeso Hombre", "Agavinas", NAH)#[-1])
names(NAH) = c("GRUPO","TRATAMIENTO","fit")
NPF = data.frame("Normopeso Mujer", "Placebo", NPF)#[-1])
names(NPF) = c("GRUPO","TRATAMIENTO","fit")
NPH = data.frame("Normopeso Hombre", "Placebo", NPH)#[-1])
names(NPH) = c("GRUPO","TRATAMIENTO","fit")
OAF = data.frame("Obeso Mujer", "Agavinas", OAF)#[-1])
names(OAF) = c("GRUPO","TRATAMIENTO","fit")
OAH = data.frame("Obeso Hombre", "Agavinas", OAH)#[-1])
names(OAH) = c("GRUPO","TRATAMIENTO","fit")
OPF = data.frame("Obeso Mujer", "Placebo", OPF)#[-1])
names(OPF) = c("GRUPO","TRATAMIENTO","fit")
OPH = data.frame("Obeso Hombre", "Placebo", OPH)#[-1])
names(OPH) = c("GRUPO","TRATAMIENTO","fit")

graph <- rbind(NAH, NPH, NAF, NPF, OAH, OPH, OAF, OPF)

ggplot(graph, aes(x = GRUPO, y = fit, fill=TRATAMIENTO))  +
stat_boxplot(geom ='errorbar') +
geom_boxplot() +
scale_y_continuous(limits=c(0,1)) +
labs(x = "Grupo", y = "Probabilidad de mejorar",
     title = "Predicciones del modelo de RL") +
scale_fill_manual(values=c("#33CC33", "red")) +
theme(legend.position="bottom")


