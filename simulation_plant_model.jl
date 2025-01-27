
include("functions.jl")

import Pkg
using Random
using LinearAlgebra
using Plots
using Distributions


# Paramètres de simulation : 
T = 10.0
NX_0 = 10  # Nombre initial de particules X
NY_0 = 0  # Nombre initial de particules Y
NZ_0 = 0  # Nombre initial de particules Z
dom = 2.0  # Domaine de simulation
lamb = 10  # Paramètre de mort de la particule X
P = 1000.  # Paramètre de précipitation
K = 0.001  # Paramètre de saturation
L = 0.01  # Paramètre d'évaporation
alpha = 100.  # Paramètre d'infiltration dans le sous sol
C = [lamb + 10., 10 / K , L + 10, alpha * 10]  # Majorants de l'intensité du processus de poisson
rG = 1.  # Rayon de voisinage pour la densité de particules Y
rc = [0.5, 0.1, 10]  # Rayon de voisinage pour la densité de particules X Y Z
diff = [0.001, 1, 5] #paramètres de diffusion pour chaque type de particules

# Initialisation des positions
# X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0 = init(NX_0, NY_0, NZ_0, L)
# X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0 = circular_patch_initialization(NX_0, NY_0, NZ_0, (-1, -1), 0.5)
X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0 = two_circular_patch_initialization(NX_0, NY_0, NZ_0, (-1, 1), (1, -1), 0.5, 0.5)

X0 = [X_x0 ,Y_x0, Z_x0]
Y0 = [X_y0 ,Y_y0, Z_y0]

sim_plant=simulate_plant(T, X0, Y0, lamb, P, K, L,I, C,rG,rc,diff; dom= 1., eps=1e-4, N_max=1000)

#print nombre de naissances et de morts pour chaque type de particules
println("naissances plantes: ", sim_plant[8], "\n")
println("morts plantes: ", sim_plant[9], "\n")
println("naissances eau sous sol: ", sim_plant[10], "\n")
print("morts eau sous sol: ", sim_plant[11], "\n")
print("morts eau surface: ", sim_plant[13], "\n")


#On trace l'évolution du nombre de particules de chaque type en fonction du temps
niter = sim_plant[7]
times = [i for i in 1:niter]
plot(times, [length(sim_plant[1][i]) for i in 1:niter], label="plantes", xlabel="temps", ylabel="nombre de particules", title="Evolution du nombre de particules en fonction du temps")
plot!(times, [length(sim_plant[3][i]) for i in 1:niter], label="eau sous sol")
plot!(times, [length(sim_plant[5][i]) for i in 1:niter], label="eau surface")


positions_Xx, positions_Xy, positions_Yx, positions_Yy, positions_Zx, positions_Zy, niter = sim_plant  

L = dom  # Taille de la boîte
# Définition des positions et des légendes
positions_x = [positions_Xx positions_Yx positions_Zx]  # Concaténation des coordonnées X
positions_y = [positions_Xy positions_Yy positions_Zy]  # Concaténation des coordonnées Y
labels = ["Plante", "Eau de sous-sol", "Eau de surface"]  # Légendes pour chaque type
base_colors = [:red, :green, :blue]  # Couleurs associées à chaque type

# Création de l'animation
anim = @animate for frame_number in 1:5:size(positions_x, 1)
    p = plot(legend=:topright, xlim=(-2*L, 2*L), ylim=(-2*L, 2*L), xlabel="x", ylabel="y", title="Simulation de particules")
    for i in 1:3  # Boucle sur les trois types de particules
        scatter!(
            positions_x[frame_number, i], 
            positions_y[frame_number, i],
            label=labels[i], 
            color=base_colors[i], 
            ms=4,  # Taille des marqueurs
            legend=:topright
        )
    end
    p
end

# Enregistrement en GIF
gif(anim, "diffusion_differentes.gif", fps=60)