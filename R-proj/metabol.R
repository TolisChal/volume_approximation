load("~/volume_approximation/R-proj/A.RData")
load("~/volume_approximation/R-proj/b.RData")

library(volesti)

d=dim(A)[2]
m = dim(A)[1]

P=gen_rand_hpoly(d,m)
P$A=A*0.1
P$b=b

q=inner_ball(P, method = "lpsolve")
print(q)
q=q[1:d]
retList=round_polytope(P, method = "max_ellipsoid")
