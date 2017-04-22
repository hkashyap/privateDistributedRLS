#reading input data from file into two 1D arrays
#relative path didn't work. Please use your Absolute path for datafile
phos2<-scan("D:/WSN/R_practice/data6_2a.dat",list(0,0,0))
maty<-phos2[[1]]
matx<-cbind(phos2[[2]],phos2[[3]])
maty<-matrix(maty)

#parameter initialization

p=2
delta=0.001
lambda=1
sample=nrow(matx)

#matrix initialization

P<-diag(1,p)
P<-P/delta
W<-matrix(0,p,1)
V<-matrix(0,1,p)
W1_plot<-matrix(0,sample,1)
W2_plot<-matrix(0,sample,1)
EST_OUT_PLOT<-matrix(0,sample,1)
IP_PLOT<-matrix(0,sample,2)
MSE<-matrix(0,sample,1)

#computation

for(n in 1:sample)
{
	V[1,]<-matx[n,]				
	sample_out=maty[n,]
	est_out=V %*% W
	error=sample_out-est_out
	z=V %*% P
	g1<-(lambda+(V %*% t(z)))
	g2<-g1[1,1]
	g<-z/g2		
	W=W + error[1,1] * t(g)
	P<-(P- t(g)%*%z)/lambda
	print(W)
	W1_plot[n,1]<-W[1,1]
	W2_plot[n,1]<-W[2,1]
	EST_OUT_PLOT[n,1]<-est_out[1,1]
	IP_PLOT[n,]<-matx[n,]
	MSE[n,]<-error[1,1]*error[1,1]
}

#plotting co-efficients and expected o/p curves

par(mfrow=c(2,2))
plot(W1_plot, col="red",xlab="Iterations",ylab="Co-efficient 1")
lines(lowess(W1_plot),col="blue")
plot(W2_plot,xlab="Iterations",ylab="Co-efficient 2")
lines(lowess(W2_plot),col="red")
plot(MSE, col="red",xlab="Iterations",ylab="Squared Error", type="b")
