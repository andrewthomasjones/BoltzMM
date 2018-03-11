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

### Test
DATA <- rfvbm(1000,c(0,1,0),matrix(c(0,1,1,1,0,1,1,1,0),3,3))
FIT <- fitfvbm(DATA,c(0,1,0),matrix(c(0,1,1,1,0,1,1,1,0),3,3))
tic()
fvbmpartiald_R(DATA,FIT)
toc()
tic()
fvbmpartiald(DATA,FIT)
toc()

# Test
DATA <- rfvbm(100,c(2,4,5),matrix(0.5,3,3)+diag(-0.5,3))
FIT <- fitfvbm(DATA,c(2,4,5),matrix(0.5,3,3)+diag(-0.5,3))
COV <- fvbmcov_R(DATA,FIT)
fvbmcov(DATA,FIT)
COV


