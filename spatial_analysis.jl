using JLD2
using Plots
using LinearAlgebra
using StatsBase


function ripley_k(points, distances, area)
    n = size(points, 1)
    K = zeros(length(distances))
    
    for (i, d) in enumerate(distances)
        count = 0
        for j in 1:n
            for k in 1:n
                if j != k
                    dist = LinearAlgebra.norm(points[j, :] .- points[k, :])
                    if dist <= d
                        count += 1
                    end
                end
            end
        end
        K[i] = (area / n) * count
        
    end
    
    return K
end


function k_confidence_interval(n_points, domain, distances, n_rep)
    K_values = zeros(n_rep, length(distances))
    dom_x, dom_y = domain

    for i in 1:n_rep
        points_x = [rand(1) * dom_x for _ in 1:n_points]
        points_y = [rand(1) * dom_y for _ in 1:n_points]
        points = hcat(points_x, points_y)

        K_values[i, :] = ripley_k(points, distances, 1.0)
    end
    
    K_mean = StatsBase.mean(K_values, dims=1)
    K_std = StatsBase.std(K_values, dims=1)
    
    K_upper = K_mean .+ 1.96 .* K_std
    K_lower = K_mean .- 1.96 .* K_std
    
    return K_mean[1, :], K_upper[1, :], K_lower[1, :]
end


@load "generated/plants_ts.jld2" plants_x_ts plants_y_ts

plants_x, plants_y = plants_x_ts[end], plants_y_ts[end]

plants_points = zeros(size(plants_x)[1], 2)
plants_points = hcat(plants_x, plants_y)

n_plants = size(plants_x)[1]
println("Nombre de plantes dans l'écosystème = ", n_plants)

domain = [maximum(plants_x) - minimum(plants_x), maximum(plants_y) - minimum(plants_y)]
println("Domaine spatiale = ", domain[1], " x ", domain[2])
max_dist = LinearAlgebra.norm(domain)
println("Distance maximale envisageable dans ce domaine = ", max_dist)

step = 0.1
distances = 0.1:step:max_dist

k_values = ripley_k(plants_points, distances, 1.0)

n_rep = 100
spatial_random_k, upper_k, lower_k = k_confidence_interval(n_plants, domain, distances, n_rep)

plot(distances, k_values, xlabel="s", ylabel="K(s)", label="K écosystème", lw=2)
plot!(
    distances, spatial_random_k,
    ribbon=(upper_k .- spatial_random_k, spatial_random_k .- lower_k), lw=2,
    label="Intervalle de confiance"
    )
plot!(distances, pi * distances.^2, label="K sous hypothèse CSR", lw=2)
