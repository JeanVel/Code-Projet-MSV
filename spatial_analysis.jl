using JLD2, Random, Plots, Statistics, LinearAlgebra, StatsBase, ProgressMeter

# === Génération de n points aléatoires sur le domaine ===
function generate_points(n, domain_size)
    return [[domain_size * rand(), domain_size * rand()] for _ in 1:n]
end

# === Calcul de la fonction K ===
function ripley_k(points, distances, domain_size)
    n = length(points)
    area = domain_size^2
    K = zeros(length(distances))
    
    for (i, d) in enumerate(distances)
        count = 0.0
        for j in 1:n
            for k in 1:n
                if j != k
                    dist = norm(points[j] - points[k])
                    if dist <= d
                        count += 1
                    end
                end
            end
        end
        K[i] = (area / (n * (n - 1))) * count
    end

    return K
end

# === Calcul de la fonction K et de l'intervalle de confiance associé pour des points placés aléatoirement dans l'espace ===
function k_confidence_interval(n_points, distances, domain, n_rep)
    K_values = zeros(n_rep, length(distances))

    
    @showprogress for i in 1:n_rep
        points = generate_points(n_points, domain)

        K_values[i, :] = ripley_k(points, distances, domain)
    end

    K_mean = StatsBase.mean(K_values, dims=1)
    K_std = StatsBase.std(K_values, dims=1)

    K_upper = K_mean .+ 1.96 .* K_std
    K_lower = K_mean .- 1.96 .* K_std

    return K_mean[1, :], K_upper[1, :], K_lower[1, :]
end


@load "generated/initial_conditions.jld2" X_x0 X_y0

plants_initial = hcat(X_x0, X_y0)

scatter(X_x0, X_y0, label="Plantes", color=:green, markersize=4, title="Répartition initiale des plantes")

domain = 1
d = 0.1:0.1:7.07
k_values_initial = ripley_k(plants_initial, d, domain)

# p = plot(d, k_values_initial, xlabel="s", ylabel="K(s)", label="Écosystème initial", lw=2, color=:blue)
# plot!(p, d, pi .* d .^ 2, label="K(s) = π s²", lw=2, linestyle=:dash, color=:red)
# savefig("figures/graphiques/Scénarios/4/ripley_k_plot_init.png")

@load "generated/plants_ts.jld2" plants_x_ts plants_y_ts

plants_x, plants_y = plants_x_ts[end], plants_y_ts[end]

plants_points = hcat(plants_x, plants_y)

n_plants = size(plants_x)[1]
println("Nombre de plantes dans l'écosystème = ", n_plants)

domain = 5
println("Domaine spatial = [", -domain, ", ", domain, "]²")
max_dist = domain * sqrt(2)
println("Distance maximale envisageable dans ce domaine = ", max_dist)

step = 0.1
distances = 0.1:step:max_dist

k_values = ripley_k(plants_points, distances, domain)

n_rep = 10
# spatial_random_k, upper_k, lower_k = k_confidence_interval(n_plants, distances, domain, n_rep)

p = plot(distances, k_values, xlabel="s", ylabel="K(s)", label="Écosystème", lw=2, color=:blue)
# p = plot!(d, k_values_initial, label="Écosystème initial", lw=2, color=:orange)
# plot!(p,
#     distances, spatial_random_k, label="Zone neutre", color=:green, alpha=1, lw=2,
#     ribbon=(upper_k .- spatial_random_k, spatial_random_k .- lower_k), fillalpha=0.2
# )
plot!(p, distances, pi * distances.^2, label="K(s) = π s²", lw=2, linestyle=:dash, color=:red)
display(p)
savefig("figures/graphiques/Scénarios/4/ripley_k_plot.png")
