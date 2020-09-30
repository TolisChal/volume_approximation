library(volesti)

P = gen_cube(100,'H')
p = sample_points(P, n =2000, random_walk = list("walk"="adaptBiW", "walk_length"=1))