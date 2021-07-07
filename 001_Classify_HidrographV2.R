#LIMPANDO MEMoRIA E CHAMANDO BIBLIOTECAS NECESSÁRIAS
rm(list = ls())
cat("\014")
library(lubridate)
library(RcmdrMisc)
library(hydroGOF)
library(plotly)
library(dplyr)


################## sub routines #######################

HidroClass <- function(disch.input, POT) {
  
  # disch.input is a dataframe with 1st column with dates and 3rd with streamflow
  
  eventos <- data.frame()
  
  #first day is always below POT
  eventos[1,1] <- 'abaixo'
  
  dia <- 2
  for (dia in 2:nrow(disch.input)) {
    
    # last day is always below POT
    if(dia == nrow(disch.input)){
      
      eventos[dia,1] <- "abaixo" 
      
      # end point match with start point of next event    
    }else if(disch.input[dia,3] < POT && disch.input[dia+1,3] >= POT && disch.input[dia-1,3] >= POT){
      
      eventos[dia,1] <- "inicio/final" 
      
      # start point of events  
    }else if(disch.input[dia,3] < POT && disch.input[dia+1,3] >= POT && disch.input[dia-1,3] < POT ){
      
      eventos[dia,1] <- "inicio"
      
      # event occur when discharge > POT  
    }else if (disch.input[dia,3] >= POT && (disch.input[dia+1,3] < POT)){ 
      
      eventos[dia,1] <- "evento"
      
      # event occur when discharge > POT  
    }else if (disch.input[dia,3] >= POT && (disch.input[dia+1,3] >= POT)){ 
      
      eventos[dia,1] <- "evento"
      
      # final is the first time step below POT  
    }else if (disch.input[dia,3] < POT && (disch.input[dia-1,3] >= POT )){
      
      eventos[dia,1] <- "final"
      
      # following time steps below POT are classified as below  
    }else if (disch.input[dia,3] < POT && (disch.input[dia+1,3] < POT )) {
      
      eventos[dia,1] <- "abaixo" 
    }
    
    
  }
  
  #correcting intersections
  for (dia in 2:nrow(disch.input)) {
    
    # day with start/end events is a event day
    if(eventos[dia,1] == "inicio/final" ){
      
      eventos[dia,1] <- "evento"
      
      # if the day after final is start, then we are in the same event (end is change to event)    
    }else if(eventos[dia,1] == "final" && eventos[dia+1,1] == "inicio"){
      
      eventos[dia,1] <- "evento" 
      
      # if the day is start and day before is event, then we are in the same event (start is change to event) 
    }else if(eventos[dia,1] == "inicio" && eventos[dia-1,1] == "evento" ){
      
      eventos[dia,1] <- "evento"
      
     
    }
    
    
  }
  # classifing hidrograph
  disch.input <- cbind(disch.input, eventos)
  

  
  return(disch.input)
  
}


EveCaracter <- function(disch.input) {
  
  # input is output from HidroClass
  
  salva.eve <- data.frame()
  dia <-1
  conta <- 1
  
  #loop to run all hidrograph classification
  while (dia <= nrow(disch.input)) {
    
    # go to next day if class is below
    if (disch.input[dia, 4] == "abaixo") {

      
      dia <- dia + 1
      
      # save day class as start 
    } else if (disch.input[dia, 4] == "inicio") {
      
      # first column at salva.eve is start day
      salva.eve[conta,1] <- disch.input[dia, 1]
      dia <- dia + 1
      #flag to sinalize is within an event
      flag <- 0
      
      # loop within an event
      while (flag == 0) {
        
        # conditional to account event days
        if (disch.input[dia, 4] == 'evento') {
          
          # runing days
          dia <- dia + 1
          # keep flag 0 inside event
          flag <- 0
          
          # conditional when arrive the final of event    
        } else if (disch.input[dia, 4] == 'final') {
          
          # second column of salva.eve is final day 
          salva.eve[conta,2] <- disch.input[dia, 1]
          
          # salve a event in a vector and identify max and length 
          ini <- which(disch.input[,1]==salva.eve[conta,1])
          fin <- which(disch.input[,1]==salva.eve[conta,2])
          vec.eve <- disch.input[ini:fin,]
          salva.eve[conta,3] <- vec.eve[which.max(vec.eve[,3]),1]
          salva.eve[conta,4] <- vec.eve[which.max(vec.eve[,3]),3]
          
          # go to next day
          dia <- dia + 1
          # flag 1 sinalize that final of event
          flag <- 1
          
          # conditional to identify a bug
        } else {
          print('deu pau')
          flag <- 1
        }
      }
      # counter to run rows in salva.eve
      conta <- conta+1 
      
      
    }
    
  }
  
  
  salva.eve <- cbind(salva.eve, salva.eve[,3] - salva.eve[,1])
  salva.eve[,5] <- as.numeric(salva.eve[,5])
  salva.eve[,5] <- salva.eve[,5]/max(salva.eve[,5])
  
  salva.eve <- cbind(salva.eve, salva.eve[,2] - salva.eve[,1])
  salva.eve[,6] <- as.numeric(salva.eve[,6])
  salva.eve[,6] <- salva.eve[,6]/max(salva.eve[,6])
  
  salva.eve <- cbind(salva.eve, salva.eve[,2] - salva.eve[,3])
  salva.eve[,7] <- as.numeric(salva.eve[,7])
  salva.eve[,7] <- salva.eve[,7]/max(salva.eve[,7])
  
  salva.eve <- cbind(salva.eve, salva.eve[,4] - POT)
  salva.eve <- cbind(salva.eve, salva.eve[,4]/POT)
  salva.eve <- cbind(salva.eve, salva.eve[,9] / salva.eve[,5])
  salva.eve <- cbind(salva.eve, salva.eve[,9] * salva.eve[,7])
  salva.eve <- cbind(salva.eve, salva.eve[,11] * salva.eve[,10] * (salva.eve[,9]**2) )
  salva.eve <- mutate(salva.eve, quantile_rank = percent_rank(salva.eve[,9]))

  
  colnames(salva.eve) <- c("ti", 
                           "tf", 
                           "tp",
                           "pico",
                           "tp - ti",
                           "tf - ti",
                           "tf - tp",
                           "dif exced",
                           "mult exced", 
                           "abrup",
                           "pers", 
                           "severit",
                           "rnk percent")
  
  return(salva.eve)
}



#set path to save files
setwd('D:/0005_R/rotinas/ClassificacaoCheias/01_classify_hydrograph/output')

caminho <- 'D:/0005_R/rotinas/ClassificacaoCheias/01_classify_hydrograph/input'


# fluviometric point names
nomes <- c("paraitinga", 
           "biritiba", 
           "jundiai", 
           "taiacupeba", 
           "billings", 
           "edgard", 
           "guarapiranga", 
           "ponte.nova", 
           "paiva.castro")


i=1
total <- data.frame()
row.ini <- 1
row.end <- 0
for (i in 1:length(nomes)) {
  
  # data source
  disch.input <-
    read.csv(file = paste0(caminho,"/chuva_vazao_", nomes[i],".csv"),
             header = TRUE, 
             sep = " ")

  
  # transform strings in date
  disch.input[,1] <-as.Date(disch.input[,1], format = "%Y-%m-%d")
  dates <- as.data.frame(disch.input[,1])
  
  # cutting interest period
  disch.input <- disch.input[123:nrow(disch.input),]
  
  # defining Point Of Threshold (POT)
  POT <- quantile (disch.input[,3], 0.85, na.rm = T)
  
  # converting NA in -99999
  disch.input[is.na(disch.input[,3])==TRUE,3] <- -99999
  
  # hidrograph classification 
  saida <- HidroClass(disch.input, POT)
  
  # characteristic of each flood events 
  eventos <- EveCaracter(saida)
  
  # correcting colnames in hidrograph output
  colnames(saida) <- c("time", "P", "Q", "Class")
  
  # write output files for each fluviometric point
  # hidrograph with classification
  write.csv(saida,paste0('hidrograph_class_',nomes[i],'.csv'), 
            row.names = F)
  # events 
  write.csv(eventos,paste0('events_', nomes[i],'.csv'),
            row.names = F)
  
  print(nrow(eventos))
  print(row.ini)
  print(row.end)

  
  total[row.ini:(nrow(eventos)+row.end), 1: (ncol(eventos)+1)] <- cbind(nomes[i], eventos)
  
  row.ini <- row.ini + nrow(eventos)
  row.end <- nrow(total)
}


# events 
write.csv(total,paste0('total_events.csv'),
          row.names = F)

