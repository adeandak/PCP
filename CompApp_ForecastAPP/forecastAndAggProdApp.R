#Carga de paquetes
library(readxl)
library(ggplot2)
source('forecastAndAggProdFunctions.R')

#---------Carga de datos
#Lee los datos del archivo de excel, en este ejemplo se saltan las primeras 10 filas pues no contienen datos relevantes,
#tambien se especifica la columna en la que se encuentran los datos

#Inputs para carga de Excel
pathArchivo="DAUPSA.xlsx"
saltos=10
columnaDatos=2
conTitulos=TRUE

#Carga de archivo e impresion de resultados
excelData <- read_excel(pathArchivo, skip = saltos, col_names = conTitulos)
print(excelData)


#-------Hacemos serie de tiempo, se debe especificar el inicio o el fin de la serie de tiempo, asi como su frecuencia,
#para esta solucion es necesario indicar el fin de la serie de tiempo.

#Inputs para creacion de serie de tiempo
frecuencia=12 #ie numero de periodos en un año
anoFin=2020
mesFin=7

#Creacion de serie de tiempo e impresion de resultados
tsData <- ts(excelData[,columnaDatos], end = c(anoFin,mesFin), frequency = frecuencia)
print(tsData)

#-------Realizacion del pronostico, se identifica automaticamente si es posible pronosticar con los datos cargados,
#en ese caso se selecciona el mejor metodo para realizarlo y se devuelve la serie de tiempo del pronostico

#Inputs para realizar el pronostico
horizonte=6 #ie periodos a futuro que se desean pronosticar

#Realizacion del pronostico y graficacion de resultados
fcData=forecastToh(ts=tsData,endTs = c(anoFin,mesFin),f=frecuencia,hF=horizonte)
autoplot(tsData) + xlab("Year") + ylab("Demand") +
  autolayer(fcData, PI=FALSE, series = paste(fcData$method,"forecast"))
print(fcData)


#-------Construccion del plan agregado de produccion, se deben especificar los datos de produccion del periodo anterior,
#al igual que las condiciones iniciales de trabajadores, inventario y backorders, y los datos de costos de cada variable

#Inputs para realizar el plan agregado de produccion
modo=1  #el tipo de plan a realizar, 
        #1=base model
        #2=level workforce
        #3=no backorders
        #4=no inventory/backorders
diasLaboralesPAnterior=260
produccionPAnterior=41383
trabajadoresPAnterior=40
vectorDiasLaborales=c(21,20,23,21,22,22)
trabajadoresIniciales=35
inventarioInicial=0
backordersInicial=0
salarios=120
costoContratar=450
costoDespedir=600
costoAlmacen=5
costoBackorder=15
costoProducir=0

#en el caso de usar el modo 2 o 3 se debe especificar la fuerza laboral
fuerzaLaboralConstante=40

#Solucion del plan de produccion e impresion
aggPlan=aggProdPlan(lambda=unlist(fcData[2], use.names=FALSE),
                 tPeriods=horizonte,
                 lDays=vectorDiasLaborales,
                 W0=trabajadoresIniciales,
                 I0=inventarioInicial,
                 B0=backordersInicial,
                 wages=salarios,
                 hire=costoContratar,
                 fire=costoDespedir,
                 hold=costoAlmacen,
                 backOrder=costoBackorder,
                 pCost=costoProducir,
                 pLastPeriod=produccionPAnterior,
                 wLastPeriod=trabajadoresPAnterior,
                 dLastPeriod=diasLaboralesPAnterior,
                 levelWf = fuerzaLaboralConstante,
                 mode = modo)