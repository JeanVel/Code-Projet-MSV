include("functions.jl")

using Pkg
using JLD2
using KernelDensity
using Plots


function animates_plants(plants_positions, dom, frame_step)
    plants_x, plants_y = plants_positions

    p = plot(legend=:topright, xlim=(-2 * dom, 2 * dom), ylim=(-2 * dom, 2 * dom), xlabel="x", ylabel="y", title="Simulation de plantes")
    frame_number = length(plants_x)

    anim = @animate for frame_number in 1:frame_step:frame_number
        scatter(p, plants_x[frame_number], plants_y[frame_number], color=:green, label="Plantes", ms=4)
    end

    return anim
end


function animates_g_water_density(g_water_positions, dom, frame_step)
    g_water_x, g_water_y = g_water_positions

    p = plot(legend=:topright, xlim=(-2 * dom, 2 * dom), ylim=(-2 * dom, 2 * dom), xlabel="x", ylabel="y", title="Simulation de l'eau de sous-sol")
    frame_number = length(g_water_x)

    anim = @animate for frame_number in 1:frame_step:frame_number
        if length(g_water_x[frame_number]) > 0
            g_water_density_estimation = kde((g_water_x[frame_number], g_water_y[frame_number]))
            heatmap(
                p, g_water_density_estimation.x, g_water_density_estimation.y, g_water_density_estimation.density',
                color=:plasma, xlabel="X", ylabel="Y", label="Densité d'eau dans le sol",
                xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        else
            heatmap(
                p, [-2 * dom, 2 * dom], [-2 * dom, 2 * dom], fill(0.0, 2, 2)',
                color=:plasma, xlabel="X", ylabel="Y", label="Densité d'eau dans le sol",
                xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        end
    end

    return anim
end


function animates_plants_and_g_water_density(plants_positions, g_water_positions, dom, frame_step)
    plants_x, plants_y = plants_positions
    g_water_x, g_water_y = g_water_positions

    frame_number = length(g_water_x)

    anim = @animate for frame_number in 1:frame_step:frame_number

        if length(g_water_x[frame_number]) > 0
            g_water_density_estimation = kde((g_water_x[frame_number], g_water_y[frame_number]))
            p2 = heatmap(
                g_water_density_estimation.x, g_water_density_estimation.y, g_water_density_estimation.density',
                color=:plasma, xlabel="X", ylabel="Y", xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom),
                legend=:topright, title="Simulation de l'eau de sous-sol"
            )
        else
            p2 = heatmap(
                [-2 * dom, 2 * dom], [-2 * dom, 2 * dom], fill(0.0, 2, 2)',
                color=:plasma, xlabel="X", ylabel="Y", xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom),
                legend=:topright, title="Simulation de l'eau de sous-sol"
            )
        end
        scatter!(p2, plants_x[frame_number], plants_y[frame_number], label="Plantes", color=:green, ms=4)
        
    end

    return anim
end


function animates_ecosystem(positions_x, position_y, dom, frame_step)
    plants_x, plants_y = positions_x[1], position_y[1]
    g_water_x, g_water_y = positions_x[2], position_y[2]
    s_water_x, s_water_y = positions_x[3], position_y[3]

    frame_number = length(g_water_x)

    anim = @animate for frame_number in 1:frame_step:frame_number
        p1 = plot(legend=:topright, xlim=(-2 * dom, 2 * dom), ylim=(-2 * dom, 2 * dom), xlabel="x", ylabel="y", title="Simulation de plantes et de l'eau de sous-sol")
        p2 = plot(legend=:topright, xlim=(-2 * dom, 2 * dom), ylim=(-2 * dom, 2 * dom), xlabel="x", ylabel="y", title="Simulation de l'eau de surface")

        if length(g_water_x[frame_number]) > 0
            g_water_density_estimation = kde((g_water_x[frame_number], g_water_y[frame_number]))
            heatmap(
                p1, g_water_density_estimation.x, g_water_density_estimation.y, g_water_density_estimation.density',
                color=:plasma, xlabel="X", ylabel="Y", label="Densité d'eau dans le sol",
                xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        else
            heatmap(
                p1, [-2 * dom, 2 * dom], [-2 * dom, 2 * dom], fill(0.0, 2, 2)',
                color=:plasma, xlabel="X", ylabel="Y", label="Densité d'eau dans le sol",
                xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        end

        scatter!(p1, plants_x[frame_number], plants_y[frame_number], label="Plantes", color=:green, ms=4)

        if length(s_water_x[frame_number]) > 0
            s_water_density_estimation = kde((s_water_x[frame_number], s_water_y[frame_number]))
            heatmap(
                p2, s_water_density_estimation.x, s_water_density_estimation.y, s_water_density_estimation.density',
                color=:viridis, xlabel="X", ylabel="Y", label="Densité d'eau de surface",
                xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        else
            heatmap(
                p2, [-2 * dom, 2 * dom], [-2 * dom, 2 * dom], fill(0.0, 2, 2)',
                color=:viridis, xlabel="X", ylabel="Y", label="Densité d'eau de surface",
                xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        end

        plot(p1, p2, layout=(1, 2), size=(1200, 400))
    end

    return anim
end


function create_animation(positions_x, positions_y, dom, frame_step=5)
    plants_x = positions_x[1]
    plants_y = positions_y[1]

    g_water_x = positions_x[2]
    g_water_y = positions_y[2]

    s_water_x = positions_x[3]
    s_water_y = positions_y[3]

    p1 = plot(legend=:topright, xlim=(-2 * dom, 2 * dom), ylim=(-2 * dom, 2 * dom), xlabel="x", ylabel="y", title="Simulation de plantes et de l'eau de sous-sol")
    p2 = plot(legend=:topright, xlim=(-2 * dom, 2 * dom), ylim=(-2 * dom, 2 * dom), xlabel="x", ylabel="y", title="Simulation de l'eau de surface")

    frame_number = size(plants_x, 1)
    println("FRAME NUMBER: ", frame_number)
    anim = @animate for frame_number in 1:frame_step:frame_number
        if length(g_water_x[frame_number]) > 0
            g_water_density_estimation = kde((g_water_x[frame_number], g_water_y[frame_number]))
            heatmap(
                p1, g_water_density_estimation.x, g_water_density_estimation.y, g_water_density_estimation.density',
                color=:plasma, xlabel="X", ylabel="Y", label="Densité d'eau dans le sol",
                xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        else
            heatmap(
            p1, [-2 * dom, 2 * dom], [-2 * dom, 2 * dom], fill(0.0, 2, 2)',
            color=:plasma, xlabel="X", ylabel="Y", label="Densité d'eau dans le sol",
            xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        end

        scatter!(p1, plants_x[frame_number], plants_y[frame_number], label="Plantes", color=:green, ms=4)

        if length(s_water_x[frame_number]) > 0
        s_water_density_estimation = kde((s_water_x[frame_number], s_water_y[frame_number]))
        heatmap(
            p2, s_water_density_estimation.x, s_water_density_estimation.y, s_water_density_estimation.density',
            color=:viridis, xlabel="X", ylabel="Y", label="Densité d'eau de surface",
            xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
        )
        else
            heatmap(
            p2, [-2 * dom, 2 * dom], [-2 * dom, 2 * dom], fill(0.0, 2, 2)',
            color=:viridis, xlabel="X", ylabel="Y", label="Densité d'eau dans le sol",
            xlims=(-2 * dom, 2 * dom), ylims=(-2 * dom, 2 * dom)
            )
        end

        plot(p1, p2, layout=(1, 2), size=(1200, 400))
    end
    return anim
end


@load "generated/plants_ts.jld2" plants_x_ts plants_y_ts
@load "generated/ground_water_ts.jld2" g_water_x_ts g_water_y_ts
@load "generated/surface_water_ts.jld2" s_water_x_ts s_water_y_ts


# Création de l'animation
# anim = animates_plants([plants_x_ts, plants_y_ts], 2, 10)
# anim = animates_g_water_density([g_water_x_ts[400:600], g_water_y_ts[400:600]], 2, 10)
anim = animates_plants_and_g_water_density([plants_x_ts[400:600], plants_y_ts[400:600]], [g_water_x_ts[400:600], g_water_y_ts[400:600]], 2, 10)
# anim = create_animation([plants_x_ts[400:600], g_water_x_ts[400:600], s_water_x_ts[400:600]], [plants_y_ts[400:600], g_water_y_ts[400:600], s_water_y_ts[400:600]], 2, 5)
println("ANIMATION CREATED")
gif(anim, "figures/animations/plants_gw_density.gif", fps=10)
println("ANIMATION SAVED")
