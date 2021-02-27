library(volesti)

P = gen_skinny_cube(4)


p = sample_points(P, n=5000, random_walk = list("walk" = "SgBiW", "walk_length"=1))