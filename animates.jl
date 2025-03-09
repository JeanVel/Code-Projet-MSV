include("functions.jl")

using Pkg
using JLD2
using KernelDensity


function create_animation(positions_x, positions_y, dom, labels, base_colors, frame_step=5)
    plants_x = positions_x[1]
    plants_y = positions_y[1]

    g_water_x = positions_x[2]
    g_water_y = positions_y[2]

    s_water_x = positions_x[3]
    s_water_y = positions_y[3]

    plant_label, g_water_label, s_water_label = labels
    plant_color, g_water_color, s_water_color = base_colors

    anim = @animate for frame_number in 1:frame_step:axes(positions_x, 1)
        p = plot(legend=:topright, xlim=(-2*dom, 2*dom), ylim=(-2*dom, 2*dom), xlabel="x", ylabel="y", title="Simulation de particules")
        scatter!(p, plants_x[frame_number], plants_y[frame_number], label=plant_label, color=plant_color, ms=4)

        

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
    return anim
end

@load "generated/plants_ts.jld2" plants_x_ts plants_y_ts
@load "generated/ground_water_ts.jld2" g_water_x_ts g_water_y_ts
@load "generated/surface_water_ts.jld2" s_water_x_ts s_water_y_ts

labels = ["Plante", "Eau de sous-sol", "Eau de surface"]  # Légendes pour chaque type
base_colors = [:green, :blue, :red]  # Couleurs associées à chaque type


# Grille pour affichage
min_x, max_x = minimum(plants_x_ts[1]), maximum(plants_x_ts[1])
min_y, max_y = minimum(plants_y_ts[1]), maximum(plants_y_ts[1])

density_estimation = kde((plants_x_ts[1], plants_y_ts[1]))

p = heatmap(density_estimation.x, density_estimation.y, density_estimation.density', color=:plasma, xlabel="X", ylabel="Y", title="Densité spatiale (KDE)",
                xlims=(min_x, max_x), ylims=(min_y, max_y))
scatter!(p, plants_x_ts[1], plants_y_ts[1], xlabel="X", ylabel="Y", color=:green)

display(p)