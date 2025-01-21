
include("functions.jl")

import Pkg
using Random
using LinearAlgebra
using Plots
using Distributions
using Pkg
using StatsBase



# Paramètres de simulation :

# Paramètres de simulation
N0 = 12  # Nombre initial de particules

br = 0.1  # Paramètre br
dr = 0.5 # Paramètre drC
dt = 0.001 # Pas de temps
max_steps = 1000  # Nombre de pas de temps
T=10.0
C= 10
# Initialisation avec des positions entre -2*L et 2*L 
L = 1.0
chi = 80.0
X0 = rand(N0) .* 4*L .- 2*L
Y0= rand(N0) .* 4*L .- 2*L
L=1


sim=simulate_particles_bd(T,X0,Y0, chi, L, br, dr,C)


print("naissances: ", sim[3], "\n")
print("morts: ", sim[4], "\n")
print("niter: ", sim[5], "\n")



# Paramètres de simulation
positions_x, positions_y = sim

# Création de l'animation
anim = @animate for frame_number in 1:size(positions_x, 1)
    scatter(positions_x[frame_number, :], positions_y[frame_number, :],
            xlim=(-2*L, 2*L), ylim=(-2*L, 2*L), legend=false,
            xlabel="x", ylabel="y", title="Simulation de particules")
end

# Enregistrement en GIF
gif(anim, "particles_simulation.gif", fps=600)