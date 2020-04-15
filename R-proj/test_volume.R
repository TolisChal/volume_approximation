library(volesti)

d=30
P = gen_rand_zonotope(d,2*d)
tim = system.time({ vol = volume(P, parameters = list("Window"=200, "verbose"=TRUE, "hpoly"=TRUE)) })
#path = system.file('extdata', package = 'volesti')

