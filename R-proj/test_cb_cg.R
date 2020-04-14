library(volesti)

N=10

vol_vec=matrix(0,1,N)
steps_vec = matrix(0,1,N)
boc_vec = matrix(0,1,N)
nballs_vec = matrix(0,1,N)
time_vec=matrix(0,1,N)
d=500
  
  P = gen_cube(d, 'H')
  
  tims1=0
  tims2=0
  count = 1
  for (j in 1:N) {
    print(j)
    tim = system.time({ vol = test_volume(P, random_walk = "BiW", walk_length = 1, parameters = list("diameter" = 2*sqrt(d), "Window" = 250, "log_length"=TRUE, "verbose"=TRUE)) })
    #tim2 = system.time({ vol = volume(P, algo='CG') })
    
    vol_vec[count,j] = vol[1]
    nballs_vec[count,j] = vol[2]
    boc_vec[count,j] = vol[3]
    steps_vec[count,j] = vol[5]
    time_vec[count,j] = as.numeric(tim)[3]
    
    save(vol_vec, file = "vol_cb_cube_500_1000.RData")
    save(nballs_vec, file = "nballs_cb_cube_500_1000.RData")
    save(boc_vec, file = "boc_cb_cube_500_1000.RData")
    save(time_vec, file = "times_cb_cube_500_1000.RData")
    save(steps_vec, file = "steps_cb_cube_500_1000.RData")
  
  }
  
 
  