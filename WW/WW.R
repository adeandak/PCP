#@authors: adeandak, sgonzalezab, jhernandezott
#
#' @title Wagner-Whitin algorithm
#'
#' @description  \code{WW} implements the Wagner-Whitin algorithm. Considering time-varying demand, costs, holding costs and setup cost, 
#'the algorithm builds production plans that minimizes the total setup and holding costs
#'in a finite horizon of time, assuming zero starting inventory
#' and no backlogging
#'
#' @param lambda  A numeric vector containing the demand per unit time
#' @param c  A numeric vector for the set-up cost per unit and period
#' @param h  A numeric vector for the holding cost per unit and period
#' @param K  A numeric vector for the costs per unit
#'
#' @return  cost The total cost of the optimal plan
#' @return  orders The order schedule of the optimal plan
#' @return  costMatrix The individual cost of all possible ordering combinations
#' @return  auxMatrix Breakdown of the total cost of the optimal plan by period and cost type
#' @return  analysis String indicating optimal ordering policy
#' @return  petersonSilver Value of the Peterson-Silver´s Rule.

WW=function(lambda = c(0,0,0,0,0,0),
            h = c(0,0,0,0,0,0),
            c = c(0,0,0,0,0,0),
            K = c(0,0,0,0,0,0)){
  #Peterson-Silver
  V=length(lambda)*sum(lambda^2)/(sum(lambda))^2-1
  if(V<0.25){
    print("Peterson-Silver es menor a 0.25, se recomienda usar EPQ.")
  }

  #inicio de la funcion
  n=length(lambda) #numero de periodos a considerar
  minAnt=0
  CM=matrix(data = 0, nrow = n, ncol = n)
  
  hc=matrix(data = 0, nrow = n, ncol = n)
  for(ren in 1:(n-1)){
    for(col in (ren+1):n){
      hc[ren,col]=hc[ren,(col-1)]+sum(h[ren:(col-1)])*lambda[col]
    }
  }
  
  for(i in 1:n){  #para recorrer los renglones
    for(j in i:n){ #para recorrer las columnas
      CM[i,j]=minAnt+K[i]+c[i]*sum(lambda[i:j])+hc[i,j]   #guardamos el costo en la matriz
    }
    minAnt=min(CM[1:i,i])  #guardamos el minimo de la columna anterior
  }
  
  #resultados
  CostoTotal=minAnt
  
  i=n
  Q=rep(0,n)
  auxString=rep("",n)
  while(i>0){
    minI=which.min(CM[1:i,i])
    auxString[i]=paste("Se pide en ",minI," para ",minI," a ",i)
    Q[minI]=sum(lambda[minI:i])
    i=minI-1
  }
  auxString=auxString[auxString!=""]
  
  l=Q-lambda
  for(k in 2:n){
    l[k]=l[k]+l[k-1]
  }
  cQ=c*Q
  hI=h*l
  KQ=K*(Q!=0)
  totalCost=cQ+KQ+hI
  matAnalysis=rbind(lambda,Q,l,c,cQ,KQ,hI,totalCost)
  
  res=list()
  res$cost=CostoTotal
  res$orders=Q
  res$costMatrix=CM
  res$auxMatrix=matAnalysis
  res$analysis=auxString
  res$petersonSilver=V
  return(res)
}

