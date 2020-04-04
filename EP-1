library(plyr)
library(reshape2)
library(ggplot2)

genomes <- read.csv(file = 'path')
genome <- character(68)
len_gen <- vector()

for(i in 1:68){
  genome[i] <- as.character(genomes[i,2])
  len_gen[i] <- nchar(genome[i])
}

AA_count <- vector()
AT_count <- vector()
AG_count <- vector()
AC_count <- vector()
TA_count <- vector()
TT_count <- vector()
TG_count <- vector()
TC_count <- vector()
GA_count <- vector()
GT_count <- vector()
GG_count <- vector()
GC_count <- vector()
CA_count <- vector()
CT_count <- vector()
CG_count <- vector()
CC_count <- vector()
A_count <- vector()
T_count <- vector()
G_count <- vector()
C_count <- vector()

#Este laço efetua a varredura de cada sequência, contando o número total de nucleotídeos
#e dinucleotídeos váĺidos de cada possibilidade (ignorando "N", "Y" e "W").
for(j in 1:68){
  AA_count[j] = as.integer(0)
  AT_count[j] = as.integer(0)
  AG_count[j] = as.integer(0)
  AC_count[j] = as.integer(0)
  TA_count[j] = as.integer(0)
  TT_count[j] = as.integer(0)
  TG_count[j] = as.integer(0)
  TC_count[j] = as.integer(0)
  GA_count[j] = as.integer(0)
  GT_count[j] = as.integer(0)
  GG_count[j] = as.integer(0)
  GC_count[j] = as.integer(0)
  CA_count[j] = as.integer(0)
  CT_count[j] = as.integer(0)
  CG_count[j] = as.integer(0)
  CC_count[j] = as.integer(0)
  A_count[j] = as.integer(0)
  T_count[j] = as.integer(0)
  G_count[j] = as.integer(0)
  C_count[j] = as.integer(0)
  for(k in 1:len_gen[j]){
    if(substr(genome[j], k, k)=="A") A_count[j] <- A_count[j] + 1
    else{
      if(substr(genome[j], k, k)=="T") T_count[j] <- T_count[j] + 1
      else{
        if(substr(genome[j], k, k)=="G") G_count[j] <- G_count[j] + 1
        else{
          if(substr(genome[j], k, k)=="C") C_count[j] <- C_count[j] + 1
        }
      }
    }
    if(substr(genome[j], k, k+1)=="AA") AA_count[j] <- AA_count[j] + 1
    else{
      if(substr(genome[j], k, k+1)=="AT") AT_count[j] <- AT_count[j] + 1
      else{
        if(substr(genome[j], k, k+1)=="AG") AG_count[j] <- AG_count[j] + 1
        else{
          if(substr(genome[j], k, k+1)=="AC") AC_count[j] <- AC_count[j] + 1
          else{
            if(substr(genome[j], k, k+1)=="TA") TA_count[j] <- TA_count[j] +1
            else{
              if(substr(genome[j], k, k+1)=="TT") TT_count[j] <- TT_count[j] +1
              else{
                if(substr(genome[j], k, k+1)=="TG") TG_count[j] <- TG_count[j] +1
                else{
                  if(substr(genome[j], k, k+1)=="TC") TC_count[j] <- TC_count[j] +1
                  else{
                    if(substr(genome[j], k, k+1)=="GA") GA_count[j] <- GA_count[j] +1
                    else{
                      if(substr(genome[j], k, k+1)=="GT") GT_count[j] <- GT_count[j] +1
                      else{
                        if(substr(genome[j], k, k+1)=="GC") GC_count[j] <- GC_count[j] +1
                        else{
                          if(substr(genome[j], k, k+1)=="GG") GG_count[j] <- GG_count[j] +1
                          else{
                            if(substr(genome[j], k, k+1)=="CA") CA_count[j] <- CA_count[j] +1
                            else{
                              if(substr(genome[j], k, k+1)=="CT") CT_count[j] <- CT_count[j] +1
                              else{
                                if(substr(genome[j], k, k+1)=="CG") CG_count[j] <- CG_count[j] +1
                                else{
                                  if(substr(genome[j], k, k+1)=="CC") CC_count[j] <- CC_count[j] +1
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        } 
      }
    }
  }
}

freqs <- data.frame("A_freq" = A_count[1:68]/len_gen[1:68], "T_freq" = T_count[1:68]/len_gen[1:68], "G_freq" = G_count[1:68]/len_gen[1:68], "C_freq" = C_count[1:68]/len_gen[1:68], "AA_freq" = (AA_count[1:68] * len_gen[1:68])/(A_count[1:68] * A_count[1:68]), "AT_freq" = (AT_count[1:68] * len_gen[1:68])/(A_count[1:68] * T_count[1:68]), "AG_freq" = (AG_count[1:68] * len_gen[1:68])/(A_count[1:68] * G_count[1:68]), "AC_freq" = (AC_count[1:68] * len_gen[1:68])/(A_count[1:68] * C_count[1:68]), "TA_freq" = (TA_count[1:68] * len_gen[1:68])/(T_count[1:68] * A_count[1:68]), "TT_freq" = (TT_count[1:68] * len_gen[1:68])/(T_count[1:68] * T_count[1:68]), "TG_freq" = (TG_count[1:68] * len_gen[1:68])/(T_count[1:68] * G_count[1:68]), "TC_freq" = (TC_count[1:68] * len_gen[1:68])/(T_count[1:68] * C_count[1:68]), "GA_freq" = (GA_count[1:68] * len_gen[1:68])/(G_count[1:68] * A_count[1:68]), "GT_freq" = (GT_count[1:68] * len_gen[1:68])/(G_count[1:68] * T_count[1:68]), "GG_freq" = (GG_count[1:68] * len_gen[1:68])/(G_count[1:68] * G_count[1:68]), "GC_freq" = (GC_count[1:68] * len_gen[1:68])/(G_count[1:68] * C_count[1:68]), "CA_freq" = (CA_count[1:68] * len_gen[1:68])/(C_count[1:68] * A_count[1:68]), "CT_freq" = (CT_count[1:68] * len_gen[1:68])/(C_count[1:68] * T_count[1:68]), "CG_freq" = (CG_count[1:68] * len_gen[1:68])/(C_count[1:68] * G_count[1:68]), "CC_freq" = (CC_count[1:68] * len_gen[1:68])/(C_count[1:68] * C_count[1:68]))
#freqs

standardization <- melt(freqs)
summary <- ddply(standardization, c("variable"), summarise, mean = mean(value), sd = sd(value))
#summary

for(i in 1:20){
  for(j in 1:68){
    freqs[j, i] <- (freqs[j, i] - summary[i, 2])/summary[i, 3]
  }
}

#Este loop calcula a distância entre cada par de sequências
distances <- matrix( , nrow = 68, ncol = 68)
for(i in 1:68){
  for(j in 1:68){
    distances[i, j] <- 0.0
    for(k in 2:20) distances[i, j] <- distances[i, j] + (freqs[i, k] - freqs[j, k])^2
    distances[i, j] <- sqrt(distances[i, j])
  }
}

#Este loop transforma a matriz de distâncias em uma matriz de adjacências
for(i in 1:68){
  for(j in 1:68){
    if(distances[i, j]<2.0) distances[i, j] <- 1
    else{
      distances[i, j] <- 0
    }
    if(i==j) distances[i,j] <- 0
  }
}

#Este loop calcula o grau de cada vértice  
degrees <- vector()
for(i in 1:68){
  degrees[i] <- 0
  for(j in 1:68){
    degrees[i] <- degrees[i] + distances[i, j]
  }
}

#Produz os gráficos de barras com as contagens dos graus dos vértices
barplot <- data.frame(seq = 1:68, deg = degrees[1:68])
graph <- ggplot(data=barplot, aes(x=seq, y=deg)) +
  geom_bar(stat="identity")
graph

#Imprime as ids dos genomas com maior grau na rede
for(i in 1:68){
  if(degrees[i]<51){
    print(genomes[i, 1])
  }
}
