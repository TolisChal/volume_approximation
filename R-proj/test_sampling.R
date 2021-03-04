library(volesti)

P = gen_cube(2, 'H')


p = sample_points(P, n=5, random_walk = list("walk" = "OptSgBiW", "walk_length"=1))
