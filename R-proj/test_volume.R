library(volesti)

d=100
P = gen_cube(d,'H')
#path = system.file('extdata', package = 'volesti')
#P = file_to_polytope(paste0(path,'/birk10.ine'))
vol_vec=c()
time_vec=c()
for (i in 1:1) {
  tim = system.time({ vol = test_volume(P, random_walk = "BiW", walk_length = 1, parameters = list("diameter" = 2*sqrt(d))) })
  vol_vec = c(vol_vec, vol)
  time_vec = c(time_vec, as.numeric(tim)[3])
}

print(vol)
print(tim)
