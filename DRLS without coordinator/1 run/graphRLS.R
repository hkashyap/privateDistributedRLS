#Script Distributed RLS

#forming the random graph with 30 vertices and entries in Graph matrix representing edges
N=30
graph<-matrix(0,N,N)
neighbor<-matrix(0,N,1)
neighborList<-matrix(0,N,N)
maxneighnode=0
maxneigh=0
for(i in 1:N)
{
	prob<-runif(n=N,min=0,max=1)
	count=0
	for(j in 1:N)
	{
		if((prob[j]<=0.2)&(i!=j))
		{
			count=count+1
			neighborList[i,count]=j
			graph[i,j]=1
		}
	}
	if(count>maxneigh)
	{
		maxneigh=count
		maxneighnode=i
	}
	neighbor[i,1]<-count
}
print(neighbor)
print(maxneighnode)

#parameter initialization

p=10
p_d=2
r=0.75
f=0.001
lambda=1
sample=100
part<-matrix(0,N,p_d)
for(i in 1:N)
{
	for(j in 1:p_d)
	{
		part[i,j]<-sample(1:p,1)
	}
	if(part[i,1]==part[i,2])
	{
		if(part[i,1]>=5)
		{
			part[i,2]<-sample(1:(part[i,1]-1),1)
		}
		else
		{
			part[i,2]<-sample((part[i,1]+1):10,1)
		}
	}
}

#y_k,X_k,u_k and v_k are initialised here using normal distribution generator

X<-matrix(0,sample,p)
y<-matrix(0,sample,1)
u<-matrix(0,sample,p)
for(i in 1:p)
{
	u[,i]<-rnorm(n=sample,mean=0,sd=1)
}
v<-rnorm(n=sample,mean=0,sd=1)

#defining C
co_eff<-runif(n=p,min=-5,max=5)
C<-matrix(co_eff,p,1)

#defining A
A<-diag(1,p)
A<-A*r

#initialising X_k(i.e. X_0)
X[1,]<-matrix(runif(n=p,min=-5,max=5),1,p)

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

#matrix initialization
P_arr<-array(0,dim=c(N,N,N))
for(i in 1:N)
{
	P_arr[,,i]<-diag(1,N)
}
P_arr<-P_arr/f
W_arr<-array(0,dim=c(N,1,N))


for(n in 1:sample)
{
	for(j in 1:N)
	{
		p_d_x<-(neighbor[j,]+p_d)
		#matrix initialization
		P<-matrix(0,p_d_x,p_d_x)
		for(i in 1:p_d_x)
		{
			for(k in 1:p_d_x)
			{
				P[i,k]<-P_arr[i,k,j]
			}
		}

		W<-matrix(0,p_d_x,1)
		for(i in 1:p_d_x)
		{
			W[i,1]<-W_arr[i,1,j]
		}

		V<-matrix(0,1,p_d_x)
		for(i in 1:p_d_x)
		{
			if(i<=p_d)
			{
				V[1,i]<-X[n,part[j,i]]
			}
			else
			{
				if(n>1)
				{
					V[1,i]<-Y_est_d[n-1,neighborList[j,(i-p_d)]]
				}
			}
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
		Y_est_d[n,j]<-est_out
	
		#Storing estimation results
		for(i in 1:p_d_x)
		{
			for(k in 1:p_d_x)
			{
				P_arr[i,k,j]<-P[i,k]
			}
		}

		for(i in 1:p_d_x)
		{
			W_arr[i,1,j]<-W[i,1]
		}
	}
}

avMSE_d<-matrix(0,sample,1)
for(i in 1:sample)
{
	for(j in 1:N)
	{
		avMSE_d[i]<-avMSE_d[i]+MSE_d[i,j]
	}
	avMSE_d[i]<-avMSE_d[i]/N
}

#plotting MSE curves for a few agents
par(mfrow=c(2,2))

plot(MSE_d[,1], col="red",xlab="Iterations",ylab="Squared Error in dB", type="b", main="Agent 1",pch=16)
plot(MSE_d[,7], col="red",xlab="Iterations",ylab="Squared Error in dB", type="b", main="Agent 7",pch=16)
plot(MSE_d[,16], col="red",xlab="Iterations",ylab="Squared Error in dB", type="b", main="Agent 16",pch=16)
plot(MSE_d[,maxneighnode], col="red",xlab="Iterations",ylab="Squared Error in dB", type="b", main="Agent with max. neighbours",pch=16)

#Comparing Central Estimation and Average of Agentwise Distributed Estimation plots

#creating a new window
windows()
plotMSE<-matrix(0,sample,2)
plotMSE[,1]<-MSE
plotMSE[,2]<-avMSE_d
Iterations<-matrix(0,sample,2)
Iterations[,1]<-sequence(sample)
Iterations[,2]<-sequence(sample)
matplot(Iterations,plotMSE,type="b",pch=16,main="Comparison Between Central Estimation and Average of Agentwise Distributed Estimation",xlab="Iterations",ylab="Central Estimation(Black) and Average of Agentwise Distributed Estimation(Red)in dB")
