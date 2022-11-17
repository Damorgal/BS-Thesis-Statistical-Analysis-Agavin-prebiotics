# Code in R for the statistical analysis and modeling of the exploratory randomized clinical trial: "Silva-Adame, M.B.; Martínez-Alvarado, A.; Martínez-Silva, V.A.; Samaniego-Méndez, V.; López, M.G. Agavins Impact on Gastrointestinal Tolerability-Related Symptoms during a Five-Week Dose-Escalation Intervention in Lean and Obese Mexican Adults: Exploratory Randomized Clinical Trial. Foods 2022, 11, 670. https://doi.org/10.3390/foods11050670"

# By: Diego A. Moreno-Galvan

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# Necesary libraries ------------------------------------------------------------------------------
library(ggplot2)
library(GGally)
library(tidyverse)
library(caTools)
library(broom)
library(glmnet)
library(matrixStats)

# Getting the data --------------------------------------------------------------------------------
data <- Proyecto_Agavinas[1:18]
limits <- Proyecto_Agavinas[19:31]
sintomas <- Proyecto_Agavinas[32:67]
caract <- Proyecto_Agavinas[68:69]
clasif <- rep(0, (dim(data)[1])/6)
for (i in 1:((dim(data)[1])/6)) {
  clasif[i] = data$CLASIFICACION[((i-1)*6)+1]
}

# Carry-forward imputation -------------------------------------------------------------------------
my_na <- is.na(data)
for (i in 1:dim(data)[1]) {
  for (j in 1:dim(data)[2]) {
    if (my_na[i,j]) {
      ant = i
      while (my_na[ant,j]) {
        ant = ant-1
      }
      if (data[i,]$DOSIS == 0.0) {
        sig = i
        while (my_na[sig,j]) {
          sig = sig+1
        }
        data[i,j] = data[sig,j]
      }
      else {
        data[i,j] = data[ant,j]
      }
      my_na[i,j] = FALSE
    }
  }
}

# Also you can choose a mean imputation ---------
my_na <- is.na(data)
for (i in 1:dim(data)[1]) {
  for (j in 1:dim(data)[2]) {
    if (my_na[i,j]) {
      ant = sig = i
      while (my_na[ant,j]) {
        ant = ant-1
      }
      while (my_na[sig,j]) {
        sig = sig+1
      }
      if (data[i,]$DOSIS == 0.0) {
        data[i,j] = data[sig,j]
      }
      else if(data[i,]$DOSIS == 12.0) {
        data[i,j] = data[ant,j]
      }
      else {
        data[i,j] = (data[ant,j] + data[sig,j])/2 
      }
      my_na[i,j] = FALSE
    }
  }
}

# Adding the column of classification (if a patient has better health) ----------------------------
for (i in 1:dim(data)[1]) {
  suma = 0
  for (j in 1:(dim(data)[2]-5)) {
    j = j+4
    liminf = as.numeric(limits[1,j-4])
    limsup = as.numeric(limits[2,j-4])
    
    if (j==5 | j==11 | j==16 ) {
      if ((liminf <= data[i,j]) & (data[i,j] <= limsup)) {
        suma = suma+1
      }
    }
  }
  data[i,dim(data)[2]] = suma
}

# Proposed linear combination to summarize the 5 weeks
# We propose to do a weighted average, where the weigh corresponds to the week dose ---------------
rdata <- data.frame(data[1:(dim(data)[1]/6), 2:3], 
                    data[1:(dim(data)[1]/6), 5:(dim(data)[2]-1)],
                    data[1:(dim(data)[1]/6), dim(data)[2]], data[1:(dim(data)[1]/6), dim(data)[2]],
                    data[1:(dim(data)[1]/6), dim(data)[2]], data[1:(dim(data)[1]/6), dim(data)[2]],
                    data[1:(dim(data)[1]/6), dim(data)[2]], data[1:(dim(data)[1]/6), dim(data)[2]],
                    data[1:(dim(data)[1]/6), dim(data)[2]], data[1:(dim(data)[1]/6), dim(data)[2]],
                    data[1:(dim(data)[1]/6), dim(data)[2]], data[1:(dim(data)[1]/6), dim(data)[2]],
                    data[1:(dim(data)[1]/6), dim(data)[2]], data[1:(dim(data)[1]/6), dim(data)[2]],
                    data[1:(dim(data)[1]/6), dim(data)[2]])

# Columns' names --------------------------------
names(rdata)[(dim(rdata)[2]-12):dim(rdata)[2]] = c("RANGO_INICIAL", 
                                                  "RANGO_FINAL",
                                                  "DIF_RANGOS",
                                                  "FLATULENCIAS",
                                                  "DISTENSION",
                                                  "MOV_INT",
                                                  "DOLOR",
                                                  "DIARREA",
                                                  "APETITO",
                                                  "SACIEDAD",
                                                  "GENERO",
                                                  "EDAD",
                                                  "CLASIFICACION")

# The linear combination ------------------------
for (i in 1:(dim(data)[1]/6)) { 
  s1 <- data[((i-1)*6)+1, 5:17]
  s2 <- data[((i-1)*6)+2, 5:17]
  s3 <- data[((i-1)*6)+3, 5:17]
  s4 <- data[((i-1)*6)+4, 5:17]
  s5 <- data[((i-1)*6)+5, 5:17]
  s6 <- data[((i-1)*6)+6, 5:17]
  pesos = 1 + 2.5 + 5 + 7 + 10 + 12
  com_lin <- s1 + s2*(2.5) + s3*(5) + s4*(7) + s5*(10) + s6*(12)
  com_lin = com_lin / pesos
  flat <- sintomas[i,2]*(2.5) + sintomas[i,3]*(5) +
    sintomas[i,4]*(7) + sintomas[i,5]*(10) + sintomas[i,6]*(12)
  dis <- sintomas[i,7]*(2.5) + sintomas[i,8]*(5) +
    sintomas[i,9]*(7) + sintomas[i,10]*(10) + sintomas[i,11]*(12)
  mov <- sintomas[i,12]*(2.5) + sintomas[i,13]*(5) +
    sintomas[i,14]*(7) + sintomas[i,15]*(10) + sintomas[i,16]*(12)
  dol <- sintomas[i,17]*(2.5) + sintomas[i,18]*(5) +
    sintomas[i,19]*(7) + sintomas[i,20]*(10) + sintomas[i,21]*(12)
  dia <- sintomas[i,22]*(2.5) + sintomas[i,23]*(5) +
    sintomas[i,24]*(7) + sintomas[i,25]*(10) + sintomas[i,26]*(12)
  ape <- sintomas[i,27]*(2.5) + sintomas[i,28]*(5) +
    sintomas[i,29]*(7) + sintomas[i,10]*(30) + sintomas[i,31]*(12)
  sac <- sintomas[i,32]*(2.5) + sintomas[i,33]*(5) +
    sintomas[i,34]*(7) + sintomas[i,35]*(10) + sintomas[i,36]*(12)
  
  rdata[i, 3:15] <- com_lin
  # Extra column (exploration): initial clinical variables in healthy rank
  rdata[i, 16] <- data[((i-1)*6)+1, 18]
  # Extra column (exploration): weighted average of the healthy ranks derivate
  rdata[i, 17] <- (data[((i-1)*6)+2, 18] - data[((i-1)*6)+1, 18])*(2.5) +
    (data[((i-1)*6)+3, 18] - data[((i-1)*6)+2, 18])*(5) + 
    (data[((i-1)*6)+4, 18] - data[((i-1)*6)+3, 18])*(7) +
    (data[((i-1)*6)+5, 18] - data[((i-1)*6)+4, 18])*(10) + 
    (data[((i-1)*6)+6, 18] - data[((i-1)*6)+5, 18])*(12) 
  rdata[i, 17] <- rdata[i, 17]/(pesos-1)
  # Extra column (exploration): healthy ranks difference (last and initial)
  rdata[i, 18] <- data[((i-1)*6)+6, 18] - data[((i-1)*6)+1, 18] 
  # Extra column (exploration): Symptoms
  rdata[i, 19] <- flat/(pesos-1)        # flatulence
  rdata[i, 20] <- dis/(pesos-1)         # abdominal distension
  rdata[i, 21] <- mov/(pesos-1)         # intestinal movement
  rdata[i, 22] <- dol/(pesos-1)         # pain
  rdata[i, 23] <- dia/(pesos-1)         # diarrhea
  rdata[i, 24] <- ape/(pesos-1)         # appetite
  rdata[i, 25] <- sac/(pesos-1)         # satiety
  # Characteristics
  rdata[i, 26] <- caract[i, 1]
  rdata[i, 27] <- caract[i, 2]
  rdata[i, 28] <- clasif[i]             # Clinical classification
  rdata[i, 1]  <- data[((i-1)*6)+1, 2]  # Group (by BMI: lean or obese)
  rdata[i, 2]  <- data[((i-1)*6)+1, 3]  # Treatment: Agavins or placebo
}


# We change classification to just 0 = No evidence of getting better, and 1 = Got better ----------
rdata$CLASIFICACION <- ifelse(rdata$CLASIFICACION > 0, 1, 0)
# We change lean group "N" to 1 and obese to 0 ----------------------------------------------------
rdata$GRUPO <- ifelse(rdata$GRUPO == "N", 1, 0)
# We change agavins treatment "A" to 1 and placebo to 0 -------------------------------------------
rdata$TRATAMIENTO <- ifelse(rdata$TRATAMIENTO == "A", 1, 0)
# We change femenine genre "F" to 1 and male to 0 -------------------------------------------------
rdata$GENERO <- ifelse(rdata$GENERO == "F", 1, 0)


# (recommended): We center the data to avoid bias -------------------------------------------------
aux <- rdata[3:27]
means <- colMeans(aux[sapply(aux, is.numeric)])

for (i in 1:(dim(rdata)[1])) {
  rdata[i, 3:15] <- rdata[i, 3:15] - means[1:13]     # Clinical var.
  rdata[i, 19:25] <- rdata[i, 19:25] - means[17:23]  # Symptoms
  rdata[i, 27] <- rdata[i, 27] - means[25]           # Age
}

# (not necessary): We normalize the data to avoid bias --------------------------------------------
rdata$GLUCOSA = rescale(rdata$GLUCOSA, to=c(-1,1))
rdata$UREA = rescale(rdata$UREA, to=c(-1,1))
rdata$CREATININA = rescale(rdata$CREATININA, to=c(-1,1))
rdata$ACIDO_URICO = rescale(rdata$ACIDO_URICO, to=c(-1,1))
rdata$AST = rescale(rdata$AST, to=c(-1,1))
rdata$ALT = rescale(rdata$ALT, to=c(-1,1))
rdata$AST_ALT = rescale(rdata$AST_ALT, to=c(-1,1))
rdata$COLESTEROL_TOTAL = rescale(rdata$COLESTEROL_TOTAL, to=c(-1,1))
rdata$TRIGLICERIDOS = rescale(rdata$TRIGLICERIDOS, to=c(-1,1))
rdata$HDL = rescale(rdata$HDL, to=c(-1,1))
rdata$LDL = rescale(rdata$LDL, to=c(-1,1))
rdata$VLDL = rescale(rdata$VLDL, to=c(-1,1))
rdata$INDICE_ATEROGENICO = rescale(rdata$INDICE_ATEROGENICO, to=c(-1,1))
