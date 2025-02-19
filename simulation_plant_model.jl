include("functions.jl")

import Pkg
using Random
using LinearAlgebra
using Plots
using Distributions
using FFTW
using Images
using NPZ


# Paramètres de simulation :
N_max = 20000  # nombre d'itérations de la simulation
T = 100.0  # Temps de simulation
dom = 2.0  # Domaine de simulation

NX_0 = 100  # Nombre initial de particules X
NY_0 = 0  # Nombre initial de particules Y
NZ_0 = 0  # Nombre initial de particules Z

diff = [0.001, 1, 5] # Coefficients diffusion pour chaque type de particules

lamb = 0.55  # Paramètre de mort naturelle des plantes (particules X)
K = 0.5  # Paramètre de saturation

P = 7.  # Intensité du processus de Poisson simulant la pluie sur l'environnement
raining_intensity = 60  # Nombre de particules ajoutées par un événement de pluie

L = 1.4  # Paramètre d'évaporation
C = [lamb + 10., 15 , L + 10, alpha * 10]  # Majorants de l'intensité du processus de poisson
rG = 2.  # Rayon de voisinage pour la densité de particules Y
rc = [0.5, 0.1, 10]  # Rayon de voisinage pour la densité de particules X Y Z
I_parameters = [150, 2]  # pente et ordonnée à l'origine pour la fonction infiltration de l'eau

# Initialisation des positions
X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0 = init(NX_0, NY_0, NZ_0, dom)
# X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0 = circular_patch_initialization(NX_0, NY_0, NZ_0, (-1, -1), 0.5)
# X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0 = two_circular_patch_initialization(NX_0, NY_0, NZ_0, (-1, 1), (1, -1), 0.5, 0.5)

X0 = [X_x0 ,Y_x0, Z_x0]
Y0 = [X_y0 ,Y_y0, Z_y0]

sim_plant = simulate_plant(X0, Y0, lamb, P, raining_intensity, K, L, I_parameters, C, rG, rc, diff, dom, N_max, T)
positions_Xx, positions_Xy, positions_Yx, positions_Yy, positions_Zx, positions_Zy, niter, last_time, naissances_X, morts_X, naissances_Y, morts_Y, naissances_Z, morts_Z = sim_plant

println("naissances plantes: ", naissances_X)
println("morts plantes: ", morts_X)
println("naissances eau sous sol: ", naissances_Y)
println("morts eau sous sol: ", morts_Y,)
println("morts eau surface: ", morts_Z)
println("Temps final : ", last_time,)

# On trace l'évolution du nombre de particules de chaque type en fonction du temps
times = [i for i in 1:niter+1]
labels = ["Plante", "Eau de sous-sol", "Eau de surface"]  # Légendes pour chaque type
base_colors = [:green, :blue, :red]  # Couleurs associées à chaque type

nb_plants = [length(plant_x_pos) for plant_x_pos in positions_Xx]
nb_ground_water = [length(ground_water_x_pos) for ground_water_x_pos in positions_Yx]
nb_surface_water = [length(surface_water_x_pos) for surface_water_x_pos in positions_Zx]
nb_particles = [nb_plants, nb_ground_water, nb_surface_water]

println("Nombre de points temporels : ", length(times))


particles_over_time = plot_particles_over_time(nb_particles, times, labels, base_colors)
display(particles_over_time)
savefig(particles_over_time, "figures/graphiques/particles_over_time.png")

plant_over_time = plot_plants_over_time(nb_plants, times)
display(plant_over_time)


# Définition des positions et des légendes
positions_x = [positions_Xx positions_Yx positions_Zx]  # Concaténation des coordonnées X
positions_y = [positions_Xy positions_Yy positions_Zy]  # Concaténation des coordonnées Y

# Création de l'animation
#anim = create_animation(positions_x, positions_y, dom, labels, base_colors, 3)

# Enregistrement en GIF
gif(anim, "figures/animations/tests.gif", fps=35)

# Ne pas sauvegarder si jamais la simulation se termine par la mort de toutes les plantes
if length(positions_Xx[end]) > 0
    # Sauvegarde des dernières positions de l'eau de surface au format .npz
    npzwrite("generated/surface_water_last_x_positions.npz", positions_Zx[end])
    npzwrite("generated/surface_water_last_y_positions.npz", positions_Zy[end])

    # Enregistrement de la dernière frame avec juste les plantes en noir
    plant_last_pos = get_image_from_positions(positions_Xx, positions_Xy)
    # display(plant_last_pos)
    savefig(plant_last_pos, "figures/images/plants.png")

    plants_gray_img = Gray.(load("figures/images/plants.png"))  # transforme l'image en niveau de gris
    plants_matrix = channelview(plants_gray_img)  # transforme l'image en matrice

    save("figures/images/plants.png", plants_matrix)  # sauvegarde la matrice en image
end




