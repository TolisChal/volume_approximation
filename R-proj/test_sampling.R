library(volesti)

P = gen_skinny_cube(2)


p = sample_points(P, n=1000, random_walk = list("walk" = "SgBiW", "walk_length"=5))