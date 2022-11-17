# Code in R for the statistical analysis and modeling of the exploratory randomized clinical trial: "Silva-Adame, M.B.; Martínez-Alvarado, A.; Martínez-Silva, V.A.; Samaniego-Méndez, V.; López, M.G. Agavins Impact on Gastrointestinal Tolerability-Related Symptoms during a Five-Week Dose-Escalation Intervention in Lean and Obese Mexican Adults: Exploratory Randomized Clinical Trial. Foods 2022, 11, 670. https://doi.org/10.3390/foods11050670"

# By: Diego A. Moreno-Galvan

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---

# MULTINOMIAL distribution accordance to evaluate the proposed model ------------------------------

# Groups-----------------------------------------
idx <- test.data$GRUPO==1 & test.data$TRATAMIENTO==1 # Lean-Agavins
NA_pi <- probabilities[idx]
NA_yi <- test.data$CLASIFICACION[idx]
idx <- test.data$GRUPO==1 & test.data$TRATAMIENTO==0 # Lean-Placebo
NP_pi <- probabilities[idx]
NP_yi <- test.data$CLASIFICACION[idx]
idx <- test.data$GRUPO==0 & test.data$TRATAMIENTO==1 # Obese-Agavins
OA_pi <- probabilities[idx]
OA_yi <- test.data$CLASIFICACION[idx]
idx <- test.data$GRUPO==0 & test.data$TRATAMIENTO==0 # Obese-Placebo
OP_pi <- probabilities[idx]
OP_yi <- test.data$CLASIFICACION[idx]

# Simulations of 1000 multinomial var. ----------
m <- 1000
tot <- 1000
NA_adj_tot <- 0
NP_adj_tot <- 0
OA_adj_tot <- 0
OP_adj_tot <- 0
#-----
NA_alp_tot <- 0
NP_alp_tot <- 0
OA_alp_tot <- 0
OP_alp_tot <- 0
#----
NA_bet_tot <- 0
NP_bet_tot <- 0
OA_bet_tot <- 0
OP_bet_tot <- 0

for(i_tot in 1:tot) {
    NA_sim <- NA_pi
    NP_sim <- NP_pi
    OA_sim <- OA_pi
    OP_sim <- OP_pi
    for(i in 1:m) {
      NA_sim <- cbind(NA_sim, rbinom(length(NA_pi), 1, NA_pi))
      NP_sim <- cbind(NP_sim, rbinom(length(NP_pi), 1, NP_pi))
      OA_sim <- cbind(OA_sim, rbinom(length(OA_pi), 1, OA_pi))
      OP_sim <- cbind(OP_sim, rbinom(length(OP_pi), 1, OP_pi))
    }

    # Parameter estimation --------------------------
    NA_q <- rep(0,4)
    NP_q <- rep(0,4)
    OA_q <- rep(0,4)
    OP_q <- rep(0,4)
    eps <- .001

    for(i in 2:(m+1)) {
      # Lean-Agav
      q <- sum(NA_sim[, i]==NA_yi & NA_yi==1)
      q <- c(q, sum(NA_sim[, i]!=NA_yi & NA_yi==0))
      q <- c(q, sum(NA_sim[, i]!=NA_yi & NA_yi==1))
      q <- c(q, sum(NA_sim[, i]==NA_yi & NA_yi==0))
      idx <- q==0
      q[idx] <- eps
      NA_q <- cbind(NA_q, q)
      
      # Lean-Plac
      q <- sum(NP_sim[, i]==NP_yi & NP_yi==1)
      q <- c(q, sum(NP_sim[, i]!=NP_yi & NP_yi==0))
      q <- c(q, sum(NP_sim[, i]!=NP_yi & NP_yi==1))
      q <- c(q, sum(NP_sim[, i]==NP_yi & NP_yi==0))
      idx <- q==0
      q[idx] <- eps
      NP_q <- cbind(NP_q, q)
      
      # Obese-Agav
      q <- sum(OA_sim[, i]==OA_yi & OA_yi==1)
      q <- c(q, sum(OA_sim[, i]!=OA_yi & OA_yi==0))
      q <- c(q, sum(OA_sim[, i]!=OA_yi & OA_yi==1))
      q <- c(q, sum(OA_sim[, i]==OA_yi & OA_yi==0))
      idx <- q==0
      q[idx] <- eps
      OA_q <- cbind(OA_q, q)
      
      # Obese-Plac
      q <- sum(OP_sim[, i]==OP_yi & OP_yi==1)
      q <- c(q, sum(OP_sim[, i]!=OP_yi & OP_yi==0))
      q <- c(q, sum(OP_sim[, i]!=OP_yi & OP_yi==1))
      q <- c(q, sum(OP_sim[, i]==OP_yi & OP_yi==0))
      idx <- q==0
      q[idx] <- eps
      OP_q <- cbind(OP_q, q)
    }

    # Delta parameter init --------------------------
    NA_delta <- 0
    NP_delta <- 0
    OA_delta <- 0
    OP_delta <- 0
    for(i in 2:(m+1)) {
      NA_delta <- c(NA_delta, log(NA_q[1,i]/NA_q[3,i]) + log(NA_q[4,i]/NA_q[2,i]))
      NP_delta <- c(NP_delta, log(NP_q[1,i]/NP_q[3,i]) + log(NP_q[4,i]/NP_q[2,i]))
      OA_delta <- c(OA_delta, log(OA_q[1,i]/OA_q[3,i]) + log(OA_q[4,i]/OA_q[2,i]))
      OP_delta <- c(OP_delta, log(OP_q[1,i]/OP_q[3,i]) + log(OP_q[4,i]/OP_q[2,i]))
    }

    NA_q <- NA_q[,-1]
    NP_q <- NP_q[,-1]
    OA_q <- OA_q[,-1]
    OP_q <- OP_q[,-1]

    # Alphas greater than 0 -------------------------
    NA_adj <- sum(log(NA_q[1,]/NA_q[3,]) > 0.)
    NP_adj <- sum(log(NP_q[1,]/NP_q[3,]) > 0.)
    OA_adj <- sum(log(OA_q[1,]/OA_q[3,]) > 0.)
    OP_adj <- sum(log(OP_q[1,]/OP_q[3,]) > 0.)

    NA_alp_tot <- c(NA_alp_tot, NA_adj)
    NP_alp_tot <- c(NP_alp_tot, NP_adj)
    OA_alp_tot <- c(OA_alp_tot, OA_adj)
    OP_alp_tot <- c(OP_alp_tot, OP_adj)

    # Betas greater than 0 --------------------------
    NA_adj <- sum(log(NA_q[4,]/NA_q[2,]) > 0.)
    NP_adj <- sum(log(NP_q[4,]/NP_q[2,]) > 0.)
    OA_adj <- sum(log(OA_q[4,]/OA_q[2,]) > 0.)
    OP_adj <- sum(log(OP_q[4,]/OP_q[2,]) > 0.)

    NA_bet_tot <- c(NA_bet_tot, NA_adj)
    NP_bet_tot <- c(NP_bet_tot, NP_adj)
    OA_bet_tot <- c(OA_bet_tot, OA_adj)
    OP_bet_tot <- c(OP_bet_tot, OP_adj)

    # Deltas greater than 0.4 (60% good - 40% bad) --
    NA_adj <- sum(NA_delta[-1] > 0.4)
    NP_adj <- sum(NP_delta[-1] > 0.4)
    OA_adj <- sum(OA_delta[-1] > 0.4)
    OP_adj <- sum(OP_delta[-1] > 0.4)

    NA_adj_tot <- c(NA_adj_tot, NA_adj)
    NP_adj_tot <- c(NP_adj_tot, NP_adj)
    OA_adj_tot <- c(OA_adj_tot, OA_adj)
    OP_adj_tot <- c(OP_adj_tot, OP_adj)
}

# Print results of multinomial experiment -------
mean(NA_adj_tot)
mean(NP_adj_tot)
mean(OA_adj_tot)
mean(OP_adj_tot)

mean(NA_alp_tot)
mean(NP_alp_tot)
mean(OA_alp_tot)
mean(OP_alp_tot)

mean(NA_bet_tot)
mean(NP_bet_tot)
mean(OA_bet_tot)
mean(OP_bet_tot)



# Validation Plots ################################################################################

# AST_ALT Plot ----------------------------------
x <- features$AST_ALT
y <- log(probabilities/(1-probabilities))
coef <- fit$coefficients
idx <- features$TRATAMIENTO == 1
colores <- rep('red',38)
colores[idx] = 'green'
idx <- features$GRUPO == 1
tipo <- rep(4,38)
tipo[idx] = 19
plot(x, y, pch=tipo, col=colores,
     main=expression(paste('AST/ALT centrado vs ln ', frac(p[i],1-p[i]))),
     xlab='AST/ALT centrado', ylab='Logits estimados')
lines(x, as.numeric(coef[1]) + as.numeric(coef['AST_ALT'])*x)
abline(h=0,col='blue',lty='dotted')
legend('topright', legend=c('Agavinas', 'Placebo','Normopeso', 'Obeso'),
       pch=c(15,15,4,19), col=c('green','red','black','black'),
       cex=0.7, y.intersp=0.3, text.width=0.15)

# Glucose Plot ----------------------------------
x <- features$GLUCOSA
plot(x,y,pch=tipo,col=colores,
     main=expression(paste('Glucosa centrada vs ln ', frac(p[i],1-p[i]))),
     xlab='Glucosa centrada', ylab='Logits estimados')
lines(x, as.numeric(coef[1])+ as.numeric(coef['GLUCOSA'])*x)
abline(h=0,col='blue',lty='dotted')
legend('topright', legend=c('Agavinas', 'Placebo','Normopeso', 'Obeso'),
       pch=c(15,15,4,19), col=c('green','red','black','black'),
       cex=0.7, y.intersp=0.3, text.width=2.8)

# VLDL Plot -------------------------------------
x <- features$VLDL
plot(x,y,pch=tipo,col=colores,
     main=expression(paste('VLDL centrado vs ln ', frac(p[i],1-p[i]))),
     xlab='VLDL centrado', ylab='Logits estimados')
lines(x, as.numeric(coef[1])+ as.numeric(coef['VLDL'])*x)
abline(h=0,col='blue',lty='dotted')
legend('topright', legend=c('Agavinas', 'Placebo','Normopeso', 'Obeso'),
       pch=c(15,15,4,19), col=c('green','red','black','black'),
       cex=0.7, y.intersp=0.3, text.width=5.6)

# (optional): rank diff. plot -------------------
x <- features$DIF_RANGOS
plot(x,y,pch=tipo,col=colores,main='Diferencia de rangos vs ln(pi/(1-pi)',
     xlab='Diferencia de rangos', ylab='Logits estimados')
lines(x, as.numeric(coef[1]) + as.numeric(coef['DIF_RANGOS'])*x)


# Getting-better estimated probability plot ####################################################

x <- features$CLASIFICACION
x <- x+.05*rdata$GRUPO -.025
y <- probabilities
idx <- rdata$GRUPO == 0
tipo <- rep(4,38)
tipo[idx] = 20
plot(x,y,pch=tipo,col=colores,main='Clasificación por Expertos vs Estimada',
     xlab='Clasificación de expertos', ylab='Probabilidad de mejora estimada',
     xlim=c(-0.3,1.3), ylim=c(0,1))
abline(h=0.5,col='blue',lty='dotted')
legend('topright', legend=c('Agavinas', 'Placebo','Normopeso', 'Obeso'),
       pch=c(15,15,4,19), col=c('green','red','black','black'),
       cex=0.7, y.intersp=0.3, text.width=.15)

       
       
# Getting-better estimated probability of IMAGINARY patients with AVERAGE clinical var. levels plot

coefi = fit$coefficients

# Lean-Agavins
idx <- features$GRUPO ==1 & features$TRATAMIENTO == 1 & rdata$GENERO==1#& features$GENERO==1
naf <- features[idx,]
aux <- naf[3:(dim(naf)[2]-1)]
# With de MEAN
means_naf <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_naf <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))
idx <- features$GRUPO == 1 & features$TRATAMIENTO == 1 & rdata$GENERO==0#features$GENERO==0
nah <- features[idx,]
aux <- nah[3:(dim(nah)[2]-1)]
# With de MEAN
means_nah <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_nah <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))

# Lean-Placebo
idx <- features$GRUPO == 1 & features$TRATAMIENTO == 0 &rdata$GENERO==1#& features$GENERO==1
npf <- features[idx,]
aux <- npf[3:(dim(npf)[2]-1)]
# With de MEAN
means_npf <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_npf <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))
idx <- features$GRUPO == 1 & features$TRATAMIENTO == 0 &rdata$GENERO==0#& features$GENERO==0
nph <- features[idx,]
aux <- nph[3:(dim(nph)[2]-1)]
# With de MEAN
means_nph <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_nph <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))

# Obese-Agavins
idx <- features$GRUPO == 0 & features$TRATAMIENTO == 1&rdata$GENERO==1#& features$GENERO==1
oaf <- features[idx,]
aux <- oaf[3:(dim(oaf)[2]-1)]
# With de MEAN
means_oaf <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_oaf <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))
idx <- features$GRUPO == 0 & features$TRATAMIENTO == 1&rdata$GENERO==0#& features$GENERO==0
oah <- features[idx,]
aux <- oah[3:(dim(oah)[2]-1)]
# With de MEAN
means_oah <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_oah <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))

# Obese-Placebo
idx <- features$GRUPO == 0 & features$TRATAMIENTO == 0&rdata$GENERO==1#& features$GENERO==1
opf <- features[idx,]
aux <- opf[3:(dim(opf)[2]-1)]
# With de MEAN
means_opf <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_opf <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))
idx <- features$GRUPO == 0 & features$TRATAMIENTO == 0&rdata$GENERO==0#& features$GENERO==0
oph <- features[idx,]
aux <- oph[3:(dim(oph)[2]-1)]
# With de MEAN
means_oph <- colMeans(aux[sapply(aux, is.numeric)])
# With de MEDIAN
#means_oph <- colMedians(data.matrix(aux[sapply(aux, is.numeric)]))

# Sum by genre ----------------------------------
sumaNAF = sum(means_naf[1:(length(means_naf))]*coefi[4:(length(coefi))])
sumaNAH = sum(means_nah[1:(length(means_nah))]*coefi[4:(length(coefi))])
sumaNPF = sum(means_npf[1:(length(means_npf))]*coefi[4:(length(coefi))])
sumaNPH = sum(means_nph[1:(length(means_nph))]*coefi[4:(length(coefi))])
sumaOAF = sum(means_oaf[1:(length(means_oaf))]*coefi[4:(length(coefi))])
sumaOAH = sum(means_oah[1:(length(means_oah))]*coefi[4:(length(coefi))])
sumaOPF = sum(means_opf[1:(length(means_opf))]*coefi[4:(length(coefi))])
sumaOPH = sum(means_oph[1:(length(means_oph))]*coefi[4:(length(coefi))])
inter = coefi[1]
trat = coefi[3]
grup = coefi[2]

# Women -----------------------------------------
NAF <- 1/(1+exp(-(inter+trat+grup+sumaNAF)))
NPF <- 1/(1+exp(-(inter+grup+sumaNPF)))
OAF <- 1/(1+exp(-(inter+trat+sumaOAF)))
OPF <- 1/(1+exp(-(inter+sumaOPF)))

# Men -------------------------------------------
NAH <- 1/(1+exp(-(inter+trat+grup+sumaNAH)))
NPH <- 1/(1+exp(-(inter+grup+sumaNPH)))
OAH <- 1/(1+exp(-(inter+trat+sumaOAH)))
OPH <- 1/(1+exp(-(inter+sumaOPH)))

# Graph -----------------------------------------
graph <- data.frame(c("Normopeso Mujer", "Normopeso Hombre",NA,NA),
                    c(NA,NA,"Obeso Mujer", "Obeso Hombre"),
                    c(NAF, NAH, OAF, OAH), c(NPF, NPH, OPF, OPH))
names(graph) = c("GRUPO1", "GRUPO2", "AGAVINAS", "PLACEBO")

ggplot(graph) +
  geom_point(aes(x = GRUPO2, y = AGAVINAS, color="Agavinas"), size=3, shape=8) +
  geom_point(aes(x = GRUPO1, y = AGAVINAS, color="Agavinas"), size=3) +
  geom_point(aes(x = GRUPO2, y = PLACEBO, color="Placebo"), size=3, , shape=8) +
  geom_point(aes(x = GRUPO1, y = PLACEBO, color="Placebo"), size=3) +
  #geom_line(aes(x = prob, y = curve, color="Agavinas")) +
  #scale_x_continuous(breaks = seq(0,12,1)) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1)) +
  labs(x = "Grupo", y = "Probabilidad de mejorar") +
  ggtitle("Probabilidades Estimadas del Modelo", "Promedio de los niveles por Grupo")+
  scale_colour_manual(name="Tratamiento",
                      values=c(Agavinas="darkblue", Placebo="red")) +
                      #values=c(Agavinas="#33CC33", Placebo="red")) +
  theme_bw() + theme(legend.position="bottom")
