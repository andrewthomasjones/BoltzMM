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


mat <- matrix(0,20,20)
bvec <-rep(0,20)
data <- rfvbm(1000,bvec,mat)

tic()
x<-fitfvbm_R(data,bvec,mat)
toc()

tic()
y<-fitfvbm(data,bvec,mat)
toc()


x$Mmat
y$Mmat

t(x$bvec)
y$bvec

x$pll
y$pll

