#Carga de paquetes
library(readxl)
library(fpp2)
library(ggplot2)

#---------Carga de datos
#Lee los datos del archivo de excel, saltando las primeras 10 filas pues estaban sucias:
domesticAutoProduction <- read_excel("DAUPSA.xlsx", skip = 10, sheet = "FRED Graph", col_names = TRUE)

#-------Hacemos serie de tiempo:
tsDomesticAutoProduction <- ts(domesticAutoProduction[,2], start = c(1993,1), end = c(2020,3), frequency = 12)

#-------------------------
#Análisis de la producción
#-------------------------
#-------Grafica la producción de autos domésticos:
autoplot(tsDomesticAutoProduction) + geom_point(color = "red", size =1 ) + ylab("Auto production") + 
  xlab("Year")  + theme(panel.background = element_rect(fill = "white", colour = "grey50"),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major = element_line(colour = "grey20", linetype = "dashed", size = 0.1)) +
  scale_x_continuous(breaks = scales::extended_breaks(20)) 

#------- Verificar condiciones---------
ggAcf(tsDomesticAutoProduction, lag.max = 342)
#Las autocorrelaciones son significativamente distintas de cero, por lo que podemos hacer un análisis de la prodcucción
Box.test(tsDomesticAutoProduction, lag = 1, type = "Ljung")
#Como p-value es menor a 0.5 sabemos que los datos en general no son aletaorios por lo que podemos hacer pronósticos

#------------------------------------------
#Encontrar el mejor método para pronosticar
#------------------------------------------

trainingData <- window(tsDomesticAutoProduction, end = c(2019, 9))
testData <- window(tsDomesticAutoProduction, start = c(2019,10))

#---------Hacer pronósticos-----
averageMethod <- meanf(trainingData, h=6) 
naiveMethod <- naive(trainingData, h= 6) 
seasonalNaiveMethod <- snaive(trainingData, h=6)
fcSimpleExpSmoo <- ses(trainingData, initial = c("optimal"), h= 6) 
holtForecast <- holt(trainingData, h=6)
hwAditive <- hw(trainingData, seasonal = "additive", h= 6)
hwMultiplicative <- hw(trainingData, seasonal = "multiplicative", h = 6)

#-------Checar la exactitud de los métodos
accuracy(averageMethod)
accuracy(naiveMethod)
accuracy(seasonalNaiveMethod)
accuracy(fcSimpleExpSmoo)
accuracy(holtForecast)
accuracy(hwAditive)
accuracy(hwMultiplicative)
#El método con menor RMSE es el de Holt

#-------Checar exactitud con respecto a los datos de prueba

accuracy(averageMethod, testData)
accuracy(naiveMethod, testData)
accuracy(seasonalNaiveMethod, testData)
accuracy(fcSimpleExpSmoo, testData)
accuracy(holtForecast, testData)
accuracy(hwAditive, testData)
accuracy(hwMultiplicative, testData)
#El método con mejor U de theil es el de Holt Winters Multiplicativo, sin embargo la de Holt es cercana 

#----------------------------------------
#Pronosticamos con el método seleccionado
#----------------------------------------

theSelectedMethod <- holt(tsDomesticAutoProduction, h = 6)
autoplot(tsDomesticAutoProduction) + xlab("Year") + ylab("Production") +
  autolayer(theSelectedMethod, PI=FALSE, series = "Holt")

#----Análisis de residuales

ggAcf(theSelectedMethod$residuals, lag.max = 342) 
#los residulaes parecen ser aleatorios pues sus correlaciones no son significactivamente distintas de cero.
Box.test(theSelectedMethod$residuals, lag = 1, type = "Ljung") 
#el p-value es de 0.6282 por lo que podemos inferir que los errores son aleatorios.


#-------------------
#----Pronosticar los siguientes 12 meses (hasta Enero 2021)
#------------------

theSelectedMethod2 <- holt(tsDomesticAutoProduction, h = 10)
autoplot(tsDomesticAutoProduction) + xlab("Year") + ylab("Production") +
  autolayer(theSelectedMethod2, PI=FALSE, series = "Holt")

#----Análisis de residuales

ggAcf(theSelectedMethod2$residuals, lag.max = 342) 
#los residulaes parecen ser aleatorios pues sus correlaciones no son significactivamente distintas de cero.
Box.test(theSelectedMethod2$residuals, lag = 1, type = "Ljung") 
#el p-value es de 0.6282 por lo que podemos inferir que los errores son aleatorios.
