# Code in R for the statistical analysis and modeling of the exploratory randomized clinical trial: "Silva-Adame, M.B.; Martínez-Alvarado, A.; Martínez-Silva, V.A.; Samaniego-Méndez, V.; López, M.G. Agavins Impact on Gastrointestinal Tolerability-Related Symptoms during a Five-Week Dose-Escalation Intervention in Lean and Obese Mexican Adults: Exploratory Randomized Clinical Trial. Foods 2022, 11, 670. https://doi.org/10.3390/foods11050670"

# By: Diego A. Moreno-Galvan

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# (recommended): Method to select the best model ##################################################

# Correlations ------------------------------------------------------------------------------------
view(cor(rdata[-(16:27)]))          # Clinical variables
view(cor(rdata))                    # All data (clinical, symptoms, etc.)

# We select just the most correlated data (for a first reduction) ---------------------------------
rdata2 <- data.frame(rdata[1:2],
                     rdata$GLUCOSA,
                     rdata$AST_ALT,
                     rdata$VLDL,
                     rdata$CLASIFICACION,
                     #rdata$TRIGLICERIDOS,
                     #rdata$COLESTEROL_TOTAL,
                     #rdata$CREATININA,
                     #####rdata$HDL,
                     #rdata$LDL,
                     #rdata$UREA,
                     #rdata$ACIDO_URICO,
                     #rdata$AST,
                     #rdata$RANGO_INICIAL,
                     #rdata$RANGO_FINAL,
                     #rdata$DIF_RANGOS,
                     #rdata$FLATULENCIAS,
                     #rdata$DISTENSION,
                     #rdata$DIARREA,
                     #rdata$MOV_INT,
                     #rdata$EDAD,
                     #rdata$GENERO,
                     #rdata$DOLOR, rdata$APETITO,
                     #rdata$SACIEDAD)

names(rdata2) = c(  "GRUPO",
                    "TRATAMIENTO",
                    "GLUCOSA",
                    "AST_ALT",
                    "VLDL",
                    "CLASIFICACION")


# Logistic Regression -----------------------------------------------------------------------------
features <- rdata2
fit <- glm(CLASIFICACION ~ ., family=binomial(link=logit), data=features)
summary(fit)

# Confidence intervals ----------------------------------------------------------------------------
probability = 0.9
cint <- confint(fit, level=probability)
cint

# Confidence intervals plot ---------------------
num_var <- length(cint)

graph  <- data.frame(c("Intercepto", "Grupo","Tratamiento", "Glucosa", "AST/ALT", "VLDL"),
                     cint[1:(num_var/2)])
names(graph) = c("coeficiente", "intervalo")
graph2 <- data.frame(c("Intercepto", "Grupo","Tratamiento", "Glucosa", "AST/ALT", "VLDL"),
                     cint[(num_var/2+1):num_var])
names(graph2) = c("coeficiente", "intervalo")
graph3 <- data.frame(c("Intercepto", "Grupo","Tratamiento", "Glucosa", "AST/ALT", "VLDL"),
                     fit$coefficients)
names(graph3) = c("coeficiente", "intervalo")

graph <- rbind(graph, graph2, graph3)

ggplot(graph, aes(x = coeficiente, y=intervalo))  +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill='#33CC33', color="black") +
  geom_hline(yintercept=0, color='red')+
  scale_x_discrete(limits=c("Intercepto", "Grupo","Tratamiento", "Glucosa", "AST/ALT", "VLDL"))+
  labs(x = "Coeficientes correspondientes a:", y = "Valor del coeficiente",
       title = "Intervalos de Confianza para los Parámetros Estimados",
       subtitle = "Confianza de 0.9 de probabilidad")+
  theme_light()


# (optional): another option to select the best model #############################################
modelo <- step(fit, direction = "backward")
fit=modelo
# Deviance and null deviance
dif_residuos <- fit$null.deviance - fit$deviance
# Degrees of freedom
df <- fit$df.null - fit$df.residual
p_value <- pchisq(q = dif_residuos,df = df)#, lower.tail = FALSE)
paste("Diferencia de residuos:", round(dif_residuos, 4))
paste("Grados de libertad:", df)
paste("p-value:", 1-p_value)
anova(fit, test = "Chisq")
confint(object = fit, level = 0.9)

# Predictions -------------------------------------------------------------------------------------
test.data=features
test.data = cbind(rdata$GRUPO, features)
names(test.data)[1] = "GRUPO"
probabilities = predict(fit, test.data, type="response")
y_pred <- ifelse(probabilities > 0.5, 1, 0)
sum(y_pred == test.data$CLASIFICACION)/38*100
sum(y_pred)
sum(y_pred == test.data$CLASIFICACION & y_pred==0)
sum(y_pred == test.data$CLASIFICACION & y_pred==1)
sum(test.data$CLASIFICACION)
