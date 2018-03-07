k<-10
n<-100

xval <-rnorm(k)
bvec <-rep(0,k)
mat <- matrix(0,k,k)
data <- rfvbm(n,bvec,mat)


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

