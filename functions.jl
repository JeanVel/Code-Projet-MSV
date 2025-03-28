import Pkg
using Random
using LinearAlgebra
using Plots
using Distributions
using Random
using StatsBase
using KernelDensity


function density(x,y,pos_x,pos_y,rc) #fonction qui calcule la densite de particules (peu importe le type) autour de du point (x,y) rc
    dens = 0.0
    for j in eachindex(pos_x)
        if x != pos_x[j] && y != pos_y[j]
            dx = x - pos_x[j]
            dy = y - pos_y[j]
            if dx^2 + dy^2 <= rc^2
                dens += 1
            end
        end
    end

    if length(pos_x) != 0
        dens /= (length(pos_x))  # Densité
    else 
        dens = 0
    end
    return dens

end


function simulate_poisson_process(P::Float64, T::Float64)
    times = Float64[]  # Liste pour stocker les instants
    t = 0.0  # Temps initial

    while t < T
        dt = -log(rand()) / P  # Tirage d'un intervalle exponentiel
        t += dt
        if t < T
            push!(times, t)
        end
    end

    return times
end


function init(NX_0, NY_0, NZ_0, dom)
    # NX_0 particules X aléatoirement sur [dom , dom]
    X_x0 = rand(NX_0) .* 2 * dom .- dom
    X_y0 = rand(NX_0) .* 2 * dom .- dom

    # NY_0 particules Y aléatoirement sur [dom , dom]
    Y_x0 = rand(NY_0) .* 2 * dom .- dom
    Y_y0 = rand(NY_0) .* 2 * dom .- dom

    # NZ_0 particules Z aléatoirement sur [dom , dom]
    Z_x0 = rand(NZ_0) .* 2 * dom .- dom
    Z_y0 = rand(NZ_0) .* 2 * dom .- dom

    return X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0
end


function circular_patch_initialization(nb_plant_particles, nb_ground_water_particles, nb_surface_water_particles,
    position, radius)
    random_plants_angle = rand(Uniform(0, 2 * pi), nb_plant_particles)
    random_plants_radius = rand(Uniform(0, radius), nb_plant_particles)

    plants_x0 = position[1] .+ random_plants_radius .* map(x -> cos(x), random_plants_angle)
    plants_y0 = position[2] .+ random_plants_radius .* map(x -> sin(x), random_plants_angle)

    random_gw_angle = rand(Uniform(0, 2 * pi), nb_ground_water_particles)
    random_gw_radius = rand(Uniform(0, radius), nb_ground_water_particles)
    ground_water_x0 = position[1] .+ random_gw_radius .* map(x -> cos(x), random_gw_angle)
    ground_water_y0 = position[2] .+ random_gw_radius .* map(x -> sin(x), random_gw_angle)

    random_sw_angle = rand(Uniform(0, 2 * pi), nb_surface_water_particles)
    random_sw_radius = rand(Uniform(0, radius), nb_surface_water_particles)
    surface_water_x0 = position[1] .+ random_sw_radius .* map(x -> cos(x), random_sw_angle)
    surface_water_y0 = position[2] .+ random_sw_radius .* map(x -> sin(x), random_sw_angle)

    return plants_x0, ground_water_x0, surface_water_x0, plants_y0, ground_water_y0, surface_water_y0
end


function two_circular_patch_initialization(nb_plant_particles, nb_gw_particles, nb_sw_particles,
    position_1, position_2, radius_1, radius_2)
    plant_group = nb_plant_particles ÷ 2
    remaining_plants = nb_plant_particles - plant_group

    gw_group = nb_gw_particles ÷ 2
    remaining_gw = nb_gw_particles - gw_group

    sw_group = nb_sw_particles ÷ 2
    remaining_sw = nb_gw_particles - sw_group

    p_x0, gw_x0, sw_x0, p_y0, gw_y0, sw_y0 = circular_patch_initialization(plant_group, gw_group, sw_group, position_1, radius_1)
    p_x0_2, gw_x0_2, sw_x0_2, p_y0_2, gw_y0_2, sw_y0_2 = circular_patch_initialization(remaining_plants, remaining_gw, remaining_sw, position_2, radius_2)

    return vcat(p_x0, p_x0_2), vcat(gw_x0, gw_x0_2), vcat(sw_x0, sw_x0_2), vcat(p_y0, p_y0_2), vcat(gw_y0, gw_y0_2), vcat(sw_y0, sw_y0_2)
end


function place_particle_in_radius_around_center(center, radius)
    angle = rand(Uniform(0, 2 * pi))
    r = radius * rand()
    x = center[1] + r * cos(angle)
    y = center[2] + r * sin(angle)

    return x, y
end


function multiple_circular_patch_initialization(n_particles, n_patches, r_influence, L)
    # Compute the number of patches per row and column
    patches_per_row = ceil(Int, sqrt(n_patches))
    println("Patches par ligne = ", patches_per_row)
    patches_per_col = ceil(Int, n_patches / patches_per_row)
    println("Patches par colonne = ", patches_per_col)
    particles_per_patch = n_particles ÷ n_patches

    # Compute the spacing between patches
    x_spacing = L / (patches_per_row - 1)
    y_spacing = L / (patches_per_col - 1)

    epsilon = domain * 0.1  # Small value to not place the patches on the border

    # Initialize arrays to store particle positions
    all_x = Float64[]
    all_y = Float64[]
    
    for i in 0:patches_per_row-1
        for j in 0:patches_per_col-1
            center_x = - L + i * 2 * x_spacing
            center_y = - L + j * 2 * y_spacing

            # Distribute particles among patch
            for _ in 1:particles_per_patch
                x, y = place_particle_in_radius_around_center([center_x, center_y], r_influence)
                push!(all_x, x)
                push!(all_y, y)
            end
        end
    end

    return all_x, all_y
end


function infiltration(plant_density, slope, intercept)
    return slope .* plant_density .+ intercept
end


function competition(plant_density, K)
    return plant_density ./ K
end

### Z = eau surface , Y = eau sous sol , X = plante
function simulate_plant(X0, Y0, lamb, P, raining_intensity, K, evaporation, I_param, C, r_plant_influence, r_plant_birth, diff, dom, N_max, T)
    # Initialisation des positions des particules
    
    #particules X :
    positions_Xx = [X0[1]]
    positions_Xy = [Y0[1]]

    #particules Y : 
    positions_Yx = [X0[2]]
    positions_Yy = [Y0[2]]

    #particules Z :
    positions_Zx = [X0[3]]
    positions_Zy = [Y0[3]]


    t = 0.0
    niter = 0
    naissances_X = 0
    naissances_Y = 0
    naissances_Z = 0
    morts_X = 0
    morts_Y = 0
    morts_Z = 0

    # Simule un processus de poisson d'intensité P entre 0 et T, on stocke les instants dans une liste : 
    raining_times = simulate_poisson_process(P, T) # Particules de pluie qui tombent
    
    while niter < N_max && t < T
        
        NX_t = length(positions_Xx[end])  # Nombre de particules à l'instant t
        NY_t = length(positions_Yx[end])  # Nombre de particules à l'instant t
        NZ_t = length(positions_Zx[end])  # Nombre de particules à l'instant t


        tau_xb = rand(Exponential(1 / ( C[1] * NX_t)))
        tau_xd= rand(Exponential(1 / ( C[2] * NX_t)))
        tau_y = rand(Exponential(1 / ( C[3] * NY_t )))
        tau_z = rand(Exponential(1 / ( C[4] * NZ_t )))
        
        tau=min(tau_xb,tau_xd, tau_y, tau_z) # pas de temps
        Zxb, Zxd, Zy, Zz = rand(), rand(), rand(), rand()

        if NX_t == 0 
            println("Plus de plantes")
            println("NX_t :",NX_t)
            println("niter :",niter)            
            break
        end


        dB_Xx, dB_Xy  = diff[1]*sqrt(2 * tau) * randn(NX_t), diff[1]*sqrt(2 * tau) * randn(NX_t) # Incréments browniens de taille Nt
        dB_Yx, dB_Yy  = diff[2]*sqrt(2 * tau) * randn(NY_t), diff[2]*sqrt(2 * tau) * randn(NY_t) # Incréments browniens de taille Nt
        dB_Zx, dB_Zy  = diff[3]*sqrt(2 * tau) * randn(NZ_t), diff[3]*sqrt(2 * tau) * randn(NZ_t) # Incréments browniens de taille Nt

        
        # Mise à jour des positions avec le schéma d'Euler à l'instant t + tau
        new_Xx = positions_Xx[end] .+ dB_Xx
        new_Xy = positions_Xy[end] .+ dB_Xy 
        #on ajoute les nouvelles positions à la liste
        push!(positions_Xx, new_Xx)
        push!(positions_Xy, new_Xy)

        new_Yx = positions_Yx[end] .+ dB_Yx  
        new_Yy = positions_Yy[end] .+ dB_Yy 
        #on ajoute les nouvelles positions à la liste
        push!(positions_Yx, new_Yx)
        push!(positions_Yy, new_Yy)        
        
        # Mise à jour des positions avec le schéma d'Euler à l'instant t + tau
        new_Zx = positions_Zx[end] .+ dB_Zx 
        new_Zy = positions_Zy[end] .+ dB_Zy
        #on ajoute les nouvelles positions à la liste
        push!(positions_Zx, new_Zx)
        push!(positions_Zy, new_Zy)

        # Mort et naissance des particules eau de surface :
        is_raining = any(r_t -> t < r_t < t + tau, raining_times)
        if is_raining
            for _ in 1:raining_intensity
                push!(positions_Zx[end], rand() * 2 * dom - dom)
                push!(positions_Zy[end], rand() * 2 * dom - dom)
            end
        end

        NZ_t = length(positions_Zx[end])  # mise à jour du nombre de particules d'eau de surface (Z)

        # Densité de plantes autour des particules d'eau de surface -> infiltration
        density_plants_around_sw = [density(x, y, positions_Xx[end], positions_Xy[end], r_plant_influence) for (x,y) in zip(positions_Zx[end], positions_Zy[end])]
        
        # Choix de la particule Z : 
        if NZ_t != 0  # s'il reste de l'eau de surface
            rate_vect_zd = evaporation[2] .+ infiltration(density_plants_around_sw, I_param[1], I_param[2])
            I_z = sample(1:NZ_t, Weights(rate_vect_zd))

            rate_I_zd = rate_vect_zd[I_z]
            b_z = rate_I_zd / C[4]
        
            if b_z >= 1
                error("Erreur : b > 1, choisissez un meilleur majorant C4.")
            end

            if Zz <= b_z  # Mort
            
                #ajouter une particule de type Y à la place de la particule Z
                push!(positions_Yx[end], positions_Zx[end][I_z])
                push!(positions_Yy[end], positions_Zy[end][I_z])
                naissances_Y += 1

                # Supprimer la particule I_z
                deleteat!(positions_Zx[end], I_z)
                deleteat!(positions_Zy[end], I_z)
                morts_Z += 1

            end
        end 


        # Mort et naissance eau de sous-sol :
        
        # Densité de plantes autour des particule d'eau dans le sol -> consommation
        rate_vect_yd = [evaporation[1] .+ density(x,y,positions_Xx[end], positions_Xy[end], r_plant_influence) for (x,y) in zip(positions_Yx[end],positions_Yy[end])]
        NY_t=length(positions_Yx[end])

        if NY_t!=0 #si il y a de l'eau en surface
            I_y = sample(1:NY_t, Weights(rate_vect_yd))
            rate_I_yd = rate_vect_yd[I_y]
            b_y= rate_I_yd/C[3]
            if b_y >= 1
                error("Erreur : b > 1, choisissez un meilleur majorant C3.")
            end
            if Zy <= b_y  # Mort
                # Supprimer la particule I_y
                deleteat!(positions_Yx[end], I_y)
                deleteat!(positions_Yy[end], I_y)
                morts_Y += 1
            end
        end 

        
        # Mort et naissance des plantes
        # Naissance 
        # Densité de plantes autour des particules d'eau de surface -> consommation
        if NX_t != 0 
            gw_density_around_plants=[density(x,y,positions_Yx[end], positions_Yy[end], r_plant_influence) for (x,y) in zip(positions_Xx[end],positions_Xy[end])]
            if sum(gw_density_around_plants) != 0 
                I_xb = sample(1:NX_t, Weights(gw_density_around_plants))
                b_xb = gw_density_around_plants[I_xb] / C[1]
                if b_xb >= 1
                    error("Erreur : b > 1, choisissez un meilleur majorant C1.")
                end
                if Zxb <= b_xb  # Naissance
                    # naissance d'une plante dans un rayon r_plant_birth autour de la plante mère
                    r = r_plant_birth * sqrt(rand())
                    theta = 2 * π * rand()
                    x_new = positions_Xx[end][I_xb] + r * cos(theta)
                    y_new = positions_Xy[end][I_xb] + r * sin(theta)
                    push!(positions_Xx[end], x_new)
                    push!(positions_Xy[end], y_new)
                    naissances_X += 1
                end
            end 

            NX_t=length(positions_Xx[end])

            # Choix de celle qui va mourir : 
            
            # Densité de plantes autour des plantes  -> compétition
            plant_density_around_plants = [density(x,y,positions_Xx[end], positions_Xy[end], r_plant_influence) for (x,y) in zip(positions_Xx[end], positions_Xy[end])]
            
            rate_vect_xd = lamb .+ competition(plant_density_around_plants, K)
            if sum(rate_vect_xd) !=0 
                I_xd=sample(1:NX_t, Weights(rate_vect_xd))
                b_xd = rate_vect_xd[I_xd] / C[2]
                if b_xd >= 1
                    error("Erreur : b > 1, choisissez un meilleur majorant C2.")
                end 
                if Zxd <= b_xd  # Mort
                    # Supprimer la particule I_xd
                    deleteat!(positions_Xx[end], I_xd)
                    deleteat!(positions_Xy[end], I_xd)
                    morts_X += 1
                end

            end
        end 
        t += tau  # Mise à jour du temps
        niter += 1
        
    end

    # Retour des résultats
    return (positions_Xx, positions_Xy, positions_Yx, positions_Yy, positions_Zx, positions_Zy, niter, t,
            naissances_X, morts_X, naissances_Y, morts_Y, naissances_Z, morts_Z)   
end


function plot_particles_over_time(nb_particles, times, labels, base_colors)
    particles_over_time = plot(xlabel="Temps", ylabel="Nombre de particules", legend=:topright)

    for (N, l, c) in zip(nb_particles, labels, base_colors)
        plot!(particles_over_time, times, N, label=l, color=c, lw=2)
    end

    return particles_over_time
end


function plot_plants_over_time(nb_plants, times, color=:green)
    plant_over_time = plot(xlabel="Temps", ylabel="Nombre de particules de plantes",legend=:none)
    plot!(plant_over_time, times, nb_plants, color=color, lw=2)

    return plant_over_time
end
