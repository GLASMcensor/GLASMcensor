# GLASMcensor: Graph Laplacian Augmented Screening Method for high dimensional censored data
```
#install.packages("devtools")
library(devtools)
install_github("GLASMcensor/GLASMcensor")
```
# Example
```
# ddata=cbind(Ti,censor,X)
ddata=ddata[order(Ti),] #The dataset is well ordered.
X=ddata[,3:(dim(ddata)[2])]
Ti=ddata[,1]
censor=ddata[,2]  #censor_i=1 means Y_i>C_i and censor_i =0 means Y_i<C_i
#amatrix is the adjacency matrix of the network structure of X.

library(GLASMcensor)
weight_C = GLASMcensor::KMestimate_CTime(Ti,censor)
measure_ori_MBKR=GLASMcensor::cMBKRall(Ti,censor,X,weight_C) # The cBKR marginal screening measure
para_GL=1
measure_GL_MBKR=GLASMcensor::abrw(amatrix,abs(uti_ori_MBKR),para_GL)  # the graph laplacian augmented screening measure
```

# 


# Development
The packages is developed based on the openmp, Rcpp and RcppEigen.

We recommend the package *huge* to do Gaussian Graph Model estimation. 
