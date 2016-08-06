# Experiment MINLP 1 - Problem 1a-1d with V EVs  - z_rob(U) - U_1 max/min
# Version 1 - Lazy constraint

using JuMP, Gurobi
#, Ipopt, AmplNLWriter

N1=60
N = 24*N1
resolution=1/N1
println("Resolution ", resolution)

# Step 1: build the model.
#m = Model(solver = GurobiSolver(OutputFlag=1, LazyConstraints=1))
m = Model(solver = GurobiSolver(LazyConstraints=1))
#m = Model(solver=IpoptSolver(print_level=1))
#m = Model(solver=IpoptSolver())
#m = Model(solver=AmplNLSolver("/Users/Nikita/bonmin-osx/bonmin"))
#m = Model(solver=BonminNLSolver())

V=40 # EV fleet size
K=3 # Number of time intervals
d=zeros(2*N)
D=zeros(2*N, N)

# Data
ptx = [0,10,20,30,40,50,60,70,80,90,100]/100
pty = [51.5 50.1 48 46 43 39.5 37 32.5 27 17.5 0]

Sinl=0.2
Qr=29.07

Mmax=[ 37.9400   35.9800   34.1600   34.0600   32.8800   35.0500   39.4400   43.7000 46.8900   52.3200   56.2500   61.7800   63.3300   68.5300   72.3800   74.3300 71.4400   67.6700   61.8000   50.7600   52.0300   50.0100   47.5200   40.9700]

Mmin=[15.7300   15.4400   14.7800   14.2300   14.6400  15.1600   17.4700   22.8700   24.8300 28.3300   30.8900   35.1400   34.9100   35.8200   35.2200   34.5100   34.3200   34.2000 33.5200   32.7800  35.4100   32.5200   26.6600   24.6300]
 
for i=1:24
	for k=1:N1
		d[N1*(i-1)+k]=Mmax[i]
		d[N+N1*(i-1)+k]=-Mmin[i]
	end
end

for i=1:N
	D[i,i]=1
	D[N+i,i]=-1
end

numpts=length(ptx)

slope=zeros(numpts-1)
b=zeros(numpts-1)
coef=zeros(numpts-1)

# Step 2: define variables
@defVar(m,  p[v=1:V, t=1:N] >=0)
@defVar(m,  0<=Y[v=1:V, k=1:K] <=1, Int)
@defVar(m, y[1:2*N] >=0 )
@defVar(m,  z)
@defVar(m,  z)
@defExpr(yd, sum{d[i]*y[i], i=1:2*N})	
@defExpr(yD[j=1:N], sum{ D[i,j]*y[i], i=1:2*N})



# Step 3: add constraints

  @addConstraint(m,  -z+ yd <= 0)

for t=1:N
	@addConstraint(m, yD[t] == sum{p[v,t], v=1:V} )
end

for v=1:V
	@addConstraint(m, p[v,1] <= 48) # value of function f at S_inl = 48
end

for k=1:(numpts-1)
	slope[k]= (pty[k+1]-pty[k])/(ptx[k+1]-ptx[k])
	b[k]=pty[k] - (pty[k+1]-pty[k])/(ptx[k+1]-ptx[k])*ptx[k]
	coef[k]=slope[k]/Qr
end

#@defExpr(aux[k=1:(numpts-1), t=1:N+1], coef[k]*sum{resolution*p[r], r=1:t-1})

#@defNLExpr(aux2[t=1:N], -0.0001*(Sinl+ 1/Qr*sum{p[k], k=1:t-1} )^3 -0.0001*(Sinl+ 1/Qr*sum{p[k], k=1:t-1} )^2 -0.1378*(Sinl+ 1/Qr*sum{p[k], k=1:t-1} ) + 51.5074  )
for t=1:N
#		@addNLConstraint(m, aux2[t]>= p[t])
end

@defExpr(aux3[v=1:V], resolution*sum{p[v,t], t=1:N})
for v=1:V
	@addConstraint(m, aux3[v] == (1-Sinl)*Qr)
end

for k=1:K  # Boundaries on the number of cars for each time period
	@addConstraint(m,sum{Y[v,k], v=1:V}<= V/2)
	@addConstraint(m,sum{Y[v,k], v=1:V}>= V/4)
end

for v=1:V
	@addConstraint(m,sum{Y[v,k], k=1:K}==1 )  # Exactly one time window should realize
end

M=51*N/3 # Big-M
for v=1:V
	@addConstraint(m,sum{p[v,k], k=2*N/3+1:N}<=M*(1-Y[v,1]) )  
	@addConstraint(m,sum{p[v,k], k=1:N/3}<=M*(1-Y[v,2]) )  
	@addConstraint(m,sum{p[v,k], k=N/3+1:2*N/3}<=M*(1-Y[v,3]) )  
end


# Step 4: add the objective.
@setObjective(m, Min, z)

#Step 4.1: Set lazy function
tolerance=0.005
num=0
function lazy(cb)
	global num
	for v=1:V
#	    println("Lazy starts!")
		phat=zeros(N) # Projection point
		coef=zeros(N) # Coefficients of the gradient
		pVal = getValue(p)
		t0=N+1 # Default value for Minimal index that violates nonlinear contraint
		for t=1:N
			s1=0
			for k=1:t-1
				s1=s1+pVal[v,k]
			end
			s2=Sinl+1/Qr*s1
			if (pVal[v,t] >= (1+tolerance)*( -0.0001*s2^3 -0.0001*s2^2 - 0.1378*s2+51.5074) ) && (t0>t)
				t0=t 
			
				# Define a projection point phat
				for k=1:t0-1 
					phat[k]=pVal[v,k]
				end
				phat[t0]=( -0.0001*s2^3 -0.0001*s2^2 - 0.1378*s2+51.5074 ) 
				for k=t0+1:N
					phat[k]=0
				end
						
				# Calculate the gradient at the point phat
				for k=1:t0-1
					coef[k]=-(-3/10000*s2^2-2/10000*s2-0.1378)/Qr
				end
				coef[t0]=1
				for k=t0+1:N
					coef[k]=0
				end
				@addLazyConstraint(cb, sum{coef[k]*p[v,k], k=1:N}<= sum{coef[k]*phat[k], k=1:N})
#				println("=== Lazy added! === ")	
				num=num+1
			end
		end	# wrt to t
	end # wrt to v
end  # lazy function

# Step 5: solve the problem!
#setLazyCallback(m,lazy)
addLazyCallback(m, lazy)
solve(m)

# Step 6: print solutions
pVal=zeros(V,N)
yVal=zeros(2*N,1)
zVal=zeros(1)

pVal=getValue(p)
zVal=getValue(z)
YVal=getValue(Y)


# for j=1:N
# 	pVal[j] = getValue(p[j])
# end
# for j=1:2*N
# 	yVal[j]=getValue(y[j])
# end
# zVal=getValue(z)

println()
println("Robust cost: ")
println(resolution*zVal)
println("Schedule: ")
for v=1:V
	println(" ")
	println("v: ", v)
	for t=1:N
		if pVal[v,t]>=0.1	
			println("t:", t, "p[t]: ", pVal[v,t])
		end
	end
end
println("Num iterations: ", num)
println("  ")
println("Y:  ",YVal)

# println("Length:", length(pVal))
