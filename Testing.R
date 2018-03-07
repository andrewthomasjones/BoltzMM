library(Matrix)

k<-4
num<-1000
mat <- matrix(rnorm(k*k),k,k)
mat<-forceSymmetric(mat)
diag(mat)<-0

mat<-as.matrix(mat)

xval <-rnorm(k)
bvec <-rnorm(k)

tic()
a<-BoltzMM::pfvbm(xval, bvec, mat)
b<-BoltzMM::allpfvbm(bvec, mat)
c<-colMeans(BoltzMM::rfvbm(num,bvec,mat))
toc()

tic()
a<-BoltzMM::allpfvbm_R(bvec, mat)
b<-BoltzMM::pfvbm_R(xval, bvec, mat)
c<-colMeans(BoltzMM::rfvbm_R(num,bvec,mat))
toc()

data<-as.matrix(iris[,1:4])
mat <- matrix(0,4,4)
bvec <-rep(0,k4)

fitfvbm_R(data,bvec,mat)
fitfvbm(data,bvec,mat)
