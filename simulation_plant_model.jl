include("functions.jl")

import Pkg
using Random
using LinearAlgebra
using Plots
using Distributions
using FFTW
using Images
using NPZ
using JLD2


# Fixe la graine aléatoire pour la reproductibilité
Random.seed!(42)

# Paramètres de simulation :
N_max = 30000  # nombre d'itérations de la simulation
T = 100.0  # Temps de simulation

NX_0 = 100  # Nombre initial de particules X
NY_0 = 0  # Nombre initial de particules Y
NZ_0 = 0  # Nombre initial de particules Z

diff = [0.001, 0.8, 3] # Coefficients diffusion pour chaque type de particules

lamb = 0.14  # Paramètre de mort naturelle des plantes
K = 15  # Capacité de charge du milieu pour les plantes

P = 1.  # Intensité du processus de Poisson simulant les événements de pluie
raining_intensity = 100  # Nombre de particules d'eau de surface ajoutées par un événement de pluie

evaporation = 1.  # Taux d'évaporation de l'eau de sous-sol
evaporation_surface = 10.
evaporations = [evaporation, evaporation_surface]

upperbound_plant_birth = 3  # plus la valeur est grande, moins l'événement "naissance d'une plante" sera favorisé 
upperbound_plant_death = 6.  # plus la valeur est grande, moins l'événement "mort d'une plante" sera favorisé
upperbound_gwater = 3.  # plus la valeur est grande, moins l'événement "mort eau de sous-sol" l'eau" sera favorisé
upperbound_swater = 130.  # plus la valeur est grande, moins l'événement "infiltration de l'eau" sera favorisé
C = [upperbound_plant_birth, upperbound_plant_death, upperbound_gwater, upperbound_swater]  # Majorants de l'intensité du processus de poisson

domain = 2.
influence_scale = 0.25
r_plant_influence = domain * influence_scale  # Rayon d'influence des plantes autour d'elle --> consommation, compétition, infiltration
println("Rayon d'influence des plantes = ", r_plant_influence)

birth_scale = 0.15
r_plant_birth = domain * birth_scale  # Rayon du cercle dans lequel une plante fille est créée autour d'une plante mère 
println("Rayon de naissance des plantes = ", r_plant_birth)

I_parameters = [50, 2]  # pente et ordonnée à l'origine pour la fonction infiltration de l'eau

# Initialisation des positions
X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0 = init(NX_0, NY_0, NZ_0, domain)
# X_x0, Y_x0 = multiple_circular_patch_initialization(800, 40, 2, 50)

s = scatter(X_x0, Y_x0, label="Plantes", color=:green, markersize=4, title="Répartition initiale des plantes")
display(s)

X0 = [X_x0 ,Y_x0, Z_x0]
Y0 = [X_y0 ,Y_y0, Z_y0]

sim_plant = simulate_plant(X0, Y0, lamb, P, raining_intensity, K, evaporations, I_parameters, C, r_plant_influence, r_plant_birth, diff, domain, N_max, T)
plants_x_ts, plants_y_ts, g_water_x_ts, g_water_y_ts, s_water_x_ts, s_water_y_ts, niter, last_time, naissances_X, morts_X, naissances_Y, morts_Y, naissances_Z, morts_Z = sim_plant


# Sauvegarde des séries temporelles des positions des particules
@save "generated/plants_ts.jld2" plants_x_ts plants_y_ts
@save "generated/ground_water_ts.jld2" g_water_x_ts g_water_y_ts
@save "generated/surface_water_ts.jld2" s_water_x_ts s_water_y_ts

println("Nombre final de plantes = ", length(plants_x_ts[end]))
println("Temps final = ", last_time,)
println("Nombre de points temporels = ", niter + 1)

# On trace l'évolution du nombre de particules de chaque type en fonction du temps
times = [i for i in 1:niter+1]
labels = ["Plante", "Eau de sous-sol", "Eau de surface"]  # Légendes pour chaque type
base_colors = [:green, :blue, :red]  # Couleurs associées à chaque type

nb_plants = [length(plant_x_pos) for plant_x_pos in plants_x_ts]
nb_ground_water = [length(g_water_x) for g_water_x in g_water_x_ts]
nb_surface_water = [length(s_water_x) for s_water_x in s_water_x_ts]
nb_particles = [nb_plants, nb_ground_water, nb_surface_water]


particles_over_time = plot_particles_over_time(nb_particles, times, labels, base_colors)
display(particles_over_time)
# savefig(particles_over_time, "figures/graphiques/Scénarios/particles_over_time.png")

plant_over_time = plot_plants_over_time(nb_plants, times)
display(plant_over_time)
# savefig(plant_over_time, "figures/graphiques/Scénarios/particles_over_time.png")

# Sauvegarde de la série temporelle du nombre de plantes au cours du temps
npzwrite("generated/plants_timeserie.npz", nb_plants)


# Ne pas sauvegarder si jamais la simulation se termine par la mort de toutes les plantes
if length(plants_x_ts[end]) > 0
    # Sauvegarde des dernières positions des plantes au format .npz
    npzwrite("generated/plants_last_x_positions.npz", plants_x_ts[end])
    npzwrite("generated/plants_last_y_positions.npz", plants_y_ts[end])

    # Affichage état écosystème
    s = scatter(plants_x_ts[end], legend=:topright, plants_y_ts[end], label="Plantes", color=:green, markersize=4, title="Répartition des plantes à la fin de la simulation")
    display_epsilon = 0.1

    for (x, y) in zip(plants_x_ts[end], plants_y_ts[end])
        plot!(s, 
            [x .+ r_plant_influence * cos.(LinRange(0, 2π, 50))], 
            [y .+ r_plant_influence * sin.(LinRange(0, 2π, 50))],
            fill=(true, 0.1), 
            lw=0, 
            color=:green, 
            label=false  # Pas de légende pour les cercles
        )
    end
    xlims!(minimum(plants_x_ts[end]) - display_epsilon, maximum(plants_x_ts[end]) + display_epsilon)
    ylims!(minimum(plants_y_ts[end]) - display_epsilon, maximum(plants_y_ts[end]) + display_epsilon)
    display(s)
end
