library(volesti)

P = gen_cube(2, 'H')


p = sample_points(P, n=5000, random_walk = list("walk" = "SupOptSgBiW", "walk_length"=1))
