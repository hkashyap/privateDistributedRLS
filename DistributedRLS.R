#Script Distributed RLS

#parameter initialization

p=10
p_d=2
r=0.75
N=10
f=0.001
lambda=1
sample=200
part<-matrix(rbind(1, 2,3,4,5,6,7,8,9,10,5,6,7,8,9,10,1,2,3,4),10,2)

#Initializing matrices average MSE calculations
avMSE_d<-matrix(0,sample,10)
avMSE<-matrix(0,sample,1)
avMSE_c<-matrix(0,sample,1)

#loop for running the RLS for 100 times

for(k in 1:100)
{

#y_k,X_k,u_k and v_k are initialised here using normal distribution generator

X<-matrix(0,sample,10)
y<-matrix(0,sample,1)
u<-matrix(0,sample,p)
for(i in 1:p)
{
	u[,i]<-rnorm(n=sample,mean=0,sd=0.5)
}
v<-rnorm(n=sample,mean=0,sd=0.5)

#defining C
co_eff<-runif(n=10,min=-5,max=5)
C<-matrix(co_eff,10,1)

#defining A
A<-diag(1,p)
A<-A*r

#initialising X_k(i.e. X_0)
X[1,]<-matrix(runif(n=10,min=-5,max=5),1,p)

#calculating o/p and state space model

for(i in 1:sample)
{
	y[i]<-(X[i,]%*%C)+v[i]
	if(i<sample)
	{
		X[i+1,]<-X[i,]%*%A
		for(j in 1:p)
		{
			X[i+1,j]<-X[i+1,j]+u[i,j]
		}
	}
}

#matrix initialization

P<-diag(1,p)
P<-P/f
W<-matrix(0,p,1)
V<-matrix(0,1,p)
MSE<-matrix(0,sample,1)

#computation using centralised RLS

for(n in 1:sample)
{
	V[1,]<-X[n,]				
	sample_out=y[n,]
	est_out=V %*% W
	error=sample_out-est_out
	z=V %*% P
	g1<-(lambda+(V %*% t(z)))
	g2<-g1[1,1]
	g<-z/g2		
	W=W + error[1,1] * t(g)
	P<-(P- t(g)%*%z)/lambda
	MSE[n,]<-error[1,1]*error[1,1]
	MSE[n,]<-10*log10(MSE[n,])
}

#Computation for DistributedRLS agentwise


MSE_d<-matrix(0,sample,N)
Y_est_d<-matrix(0,sample,N)

for(j in 1:N)
{
	#matrix initialization
	P<-diag(1,p_d)
	P<-P/f
	W<-matrix(0,p_d,1)
	V<-matrix(0,1,p_d)
	for(n in 1:sample)
	{
		for(i in 1:p_d)
		{
			V[1,i]<-X[n,part[j,i]]
		}
		sample_out=y[n,]
		est_out=V %*% W
		error=sample_out-est_out
		z=V %*% P
		g1<-(lambda+(V %*% t(z)))
		g2<-g1[1,1]
		g<-z/g2		
		W=W + error[1,1] * t(g)
		P<-(P- t(g)%*%z)/lambda
		MSE_d[n,j]<-error[1,1]*error[1,1]
		MSE_d[n,j]<-10*log10(MSE_d[n,j])
#	Y_est_d[n,j]=est_out
	}
	for(n in 1:sample)
	{
#		Y_est_d[n,j]<-Y_est_d[n,j]+(W[1,1]*X[n,part[j,1]])
#		Y_est_d[n,j]<-Y_est_d[n,j]+(W[2,1]*X[n,part[j,2]])
		Y_est_d[n,j]<-(W[1,1]*X[n,part[j,1]])+(W[2,1]*X[n,part[j,2]])	
	}

}

#Computation by central coordinator based on estimation by agents

#matrix initialization

P<-diag(1,N)
P<-P/f
W<-matrix(0,N,1)
V<-matrix(0,1,N)
MSE_c<-matrix(0,sample,1)

#computation using centralised RLS

for(n in 1:sample)
{
	V[1,]<-Y_est_d[n,]				
	sample_out=y[n,]
	est_out=V %*% W
	error=sample_out-est_out
	z=V %*% P
	g1<-(lambda+(V %*% t(z)))
	g2<-g1[1,1]
	g<-z/g2		
	W=W + error[1,1] * t(g)
	P<-(P- t(g)%*%z)/lambda
	MSE_c[n,]<-error[1,1]*error[1,1]
	MSE_c[n,]<-10*log10(MSE_c[n,])
}
avMSE_d<-avMSE_d+MSE_d
avMSE_c<-avMSE_c+MSE_c
avMSE<-avMSE+MSE
}

#calculating average matrices
avMSE_d<-avMSE_d/100
avMSE_c<-avMSE_c/100
avMSE<-avMSE/100

#plotting MSE curves for few agents
windows()
par(mfrow=c(2,2))

plot(avMSE_d[,1], col="red",xlab="Iterations",ylab="Average Squared Error in dB", type="b", main="Agent 1", pch=16)
plot(avMSE_d[,3], col="red",xlab="Iterations",ylab="Average Squared Error in dB", type="b", main="Agent 3", pch=16)
plot(avMSE_d[,7], col="red",xlab="Iterations",ylab="Average Squared Error in dB", type="b", main="Agent 7", pch=16)
plot(avMSE_d[,10], col="red",xlab="Iterations",ylab="Average Squared Error in dB", type="b", main="Agent 10", pch=16)


#creating a new window
windows()
plotMSE<-matrix(0,sample,2)
plotMSE[,1]<-avMSE
plotMSE[,2]<-avMSE_c
Iterations<-matrix(0,sample,2)
Iterations[,1]<-sequence(sample)
Iterations[,2]<-sequence(sample)
matplot(Iterations,plotMSE,type="b",pch=16,xlab="Iterations",ylab="Full information MSE(Black) and Coordinator's MSE(Red)in dB")

