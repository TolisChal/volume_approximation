library(volesti)

q = R.matlab::readMat('polytope_iAB_RBC_283.mat')
q = q$polytope

inner_ball = get_max_inner_ball(q[[1]], q[[2]])

A = q[[1]]
b = q[[2]]

d = dim(A)[2]
P = gen_cube(d, 'H')
P$A = A
P$b = b

L = 5*d*inner_ball$radius

#print("sampling with billiard - barrier")
#time2 = system.time({samples2 = sample_points(P, n = 10000, 
#                                              random_walk = list("walk" = "SupOptSgBiW", "walk_length"=1, "burnin" = 0, 
#                                                                 "starting_point"=inner_ball$center,
#                                                                 "L" = L)) })

print("sampling with classic billiard")
time1 = system.time({samples1 = sample_points(P, n = 10000, 
                              random_walk = list("walk" = "aBiW", "walk_length"=1, "burnin" = 0, 
                                                 "starting_point"=inner_ball$center,
                                                 "L" = L)) })

ess1 = coda::effectiveSize(t(samples1))
#ess2 = coda::effectiveSize(t(samples2))

eff1 = min(ess1) / as.numeric(time1)[3]
#eff2 = min(ess2) / as.numeric(time2)[3]

save(ess1,file = "ess1.RData")
#save(ess2,file = "ess2.RData")

save(eff1,file = "eff1.RData")
#save(eff2,file = "eff2.RData")

save(time1,file = "time1.RData")
#save(time2,file = "time2.RData")

save(samples1,file = "samples1.RData")
#save(samples2,file = "samples2.RData")

#print("sampling with non optimized billiard - barrier")
#time3 = system.time({samples3 = sample_points(P, n = 10000, 
#                                              random_walk = list("walk" = "SgBiW", "walk_length"=1, "burnin" = 0, 
#                                                                 "starting_point"=inner_ball$center,
#                                                                 "L" = L)) })
#ess3 = coda::effectiveSize(t(samples3))
#eff3 = min(ess3) / as.numeric(time3)[3]

#save(ess3,file = "ess3.RData")
#save(eff3,file = "eff3.RData")
#save(time3,file = "time3.RData")
#save(samples3,file = "samples3.RData")


print("sampling with non optimized billiard - barrier")
time4 = system.time({samples4 = sample_points(P, n = 10000, 
                                              random_walk = list("walk" = "gBiWHes", "walk_length"=1, "burnin" = 0, 
                                                                 "starting_point"=inner_ball$center,
                                                                 "L" = L)) })
ess4 = coda::effectiveSize(t(samples4))
eff4 = min(ess4) / as.numeric(time4)[3]

save(ess4,file = "ess4.RData")
save(eff4,file = "eff4.RData")
save(time4,file = "time4.RData")
save(samples4,file = "samples4.RData")


