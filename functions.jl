
import Pkg
using Random
using LinearAlgebra
using Plots
using Distributions
using Pkg
using Random



function simulate_particles(n_particles, T, dt, X0, Y0, chi, L, eps=1e-4)
    n_steps = Int(T / dt)
    sqrt_2_dt = sqrt(2 * dt)
    N = n_particles
    
    # Initialisation des positions
    positions_x = zeros(Float64, n_steps, n_particles)
    positions_y = zeros(Float64, n_steps, n_particles)
    positions_x[1, :] .= X0
    positions_y[1, :] .= Y0
    
    for t in 2:n_steps
        # Incréments brownien pour chaque particule en x et y
        dB_x = sqrt_2_dt * randn(n_particles)
        dB_y = sqrt_2_dt * randn(n_particles)
        
        # Calcul des interactions entre les particules
        interactions_x = zeros(Float64, n_particles)
        interactions_y = zeros(Float64, n_particles)
        for i in 1:n_particles
            for j in 1:n_particles
                if i != j
                    # Calcul de la distance sur le tore
                    dx = positions_x[t-1, i] - positions_x[t-1, j]
                    dy = positions_y[t-1, i] - positions_y[t-1, j]
                    
                    dist_squared = dx^2 + dy^2
                    interactions_x[i] -= (chi / (2 * π * N)) * dx / (dist_squared + eps)
                    interactions_y[i] -= (chi / (2 * π * N)) * dy / (dist_squared + eps)
                end
            end
        end
        
        # Mise à jour des positions avec le schéma d'Euler
        positions_x[t, :] .= (positions_x[t-1, :] .+ dB_x .+ dt .* interactions_x) 
        positions_y[t, :] .= (positions_y[t-1, :] .+ dB_y .+ dt .* interactions_y)
    end
    
    return positions_x, positions_y
end


function simulate_particles_bd(T, X0, Y0, chi, L, br, dr, C; eps=1e-4, rc=0.1, N_max=2000)
    # Initialisation des positions des particules
    positions_x = [X0]
    positions_y = [Y0]
    t = 0.0
    niter = 0
    naissances = 0
    morts = 0
    
    while t < T && niter < N_max
        niter += 1
        Nt = length(positions_x[end])  # Nombre de particules à l'instant t

        tau = rand(Exponential(1 / ((br + dr * C) * Nt)))
        Z = rand()
        I = rand(1:Nt)  # Tirage d'un indice uniforme entre 1 et Nt
        t += tau  # Mise à jour du temps


        # Calcul des interactions entre les particules
        interactions_x = zeros(Float64, Nt)
        interactions_y = zeros(Float64, Nt)
        for i in 1:Nt
            for j in 1:Nt
                if i != j
                    dx = positions_x[end][i] - positions_x[end][j]
                    dy = positions_y[end][i] - positions_y[end][j]
                    dist_squared = dx^2 + dy^2
                    interactions_x[i] -= (chi / (2 * π * Nt)) * dx / (dist_squared + eps)
                    interactions_y[i] -= (chi / (2 * π * Nt)) * dy / (dist_squared + eps)
                end
            end
        end

        dB_x = sqrt(2 * tau) * randn(Nt) # Incréments browniens de taille Nt
        dB_y = sqrt(2 * tau) * randn(Nt)
        
        # Mise à jour des positions avec le schéma d'Euler à l'instant t + tau
        new_x = positions_x[end] .+ dB_x .+ tau .* interactions_x  #de taille Nt
        new_y = positions_y[end] .+ dB_y .+ tau .* interactions_y


        #on ajoute les nouvelles positions à la liste
        push!(positions_x, new_x)
        push!(positions_y, new_y)

    
        # Calcul de la densité autour de la particule I dans un rayon rc
        dens = 0.0
        for j in 1:length(new_x)
            if j != I
                dx = new_x[I] - new_x[j]
                dy = new_y[I] - new_y[j]
                if dx^2 + dy^2 <= rc^2
                    dens += 1
                end
            end
        end
        dens /= (Nt)  # Densité 

        a = br / (br + dr * C)
        b = a + (dr * dens) / (br + dr * C)

        if b >= 1
            println("Erreur : b > 1, choisissez un meilleur majorant C.")
            break
        end

        if Z <= a  # Naissance
            # naissance d'une particule à une distance rc coté de la particule I 
            r = rc * sqrt(rand())
            theta = 2 * π * rand()
            x_new = new_x[I] + r * cos(theta)
            y_new = new_y[I] + r * sin(theta)
            push!(positions_x[end], x_new)
            push!(positions_y[end], y_new)


            naissances += 1

        elseif Z > a && Z <= b  # Mort
            # Supprimer la particule I
            deleteat!(positions_x[end], I)
            deleteat!(positions_y[end], I)
            morts += 1
        end


    end

    # Retour des résultats
    return positions_x, positions_y, naissances, morts, niter
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


function density(x,y,pos_x,pos_y,rc) #fonction qui calcule la densite de particules (peu importe le type) autour de du point (x,y) rc
    dens = 0.0
    for j in 1:length(pos_x)
        if x != pos_x[j] && y != pos_y[j]
            dx = x - pos_x[j]
            dy = y - pos_y[j]
            if dx^2 + dy^2 <= rc^2
                dens += 1
            end
        end
    end
    dens /= (length(pos_x))  # Densité
    return dens

end

function init(NX_0, NY_0, NZ_0, dom)
    # NX_0 particules X aléatoirement sur [dom , dom]
    X_x0 = rand(NX_0) .* 2*dom .- dom
    X_y0 = rand(NX_0) .* 2*dom .- dom

    # NY_0 particules Y aléatoirement sur [dom , dom]
    Y_x0 = rand(NY_0) .* 2*dom .- dom
    Y_y0 = rand(NY_0) .* 2*dom .- dom

    # NZ_0 particules Z aléatoirement sur [dom , dom]
    Z_x0 = rand(NZ_0) .* 2*dom .- dom
    Z_y0 = rand(NZ_0) .* 2*dom .- dom

    return X_x0, Y_x0, Z_x0, X_y0, Y_y0, Z_y0
end


function two_circular_patch_initialization(nb_plant_particles, nb_ground_water_particles, nb_surface_water_particles,
    position_1, radius_1, position_2, radius_2)
    

end


function elliptic_patch_initialization(nb_plant_particles, nb_ground_water_particles, nb_surface_water_particles,
    position, a, b)

end


function simulate_plant(T, X0, Y0, lamb, P, K,  L ,alpha, C,rG,rc,diff; dom= 1., eps=1e-4, N_max=2000)
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

    #simule un processus de poisson d'intensité P entre 0 et T, on stocke les instants dans une liste : 
    times = simulate_poisson_process(P, T) #Particules de pluie qui tombent 
    
    while t < T && niter < N_max
        niter += 1
        NX_t = length(positions_Xx[end])  # Nombre de particules à l'instant t
        NY_t = length(positions_Yx[end])  # Nombre de particules à l'instant t
        NZ_t = length(positions_Zx[end])  # Nombre de particules à l'instant t


        tau_x = rand(Exponential(1 / ( C[1] * NX_t)))
        tau_y = rand(Exponential(1 / ( C[2] * NY_t ) ))
        tau_z = rand(Exponential(1 / ( C[3] * NZ_t ) ))
        

        tau=min(tau_x,tau_y,tau_z) # pas de temps
        Zx , Zy, Zz = rand(),rand(),rand()


        if NX_t != 0 && NY_t != 0 && NZ_t != 0
            I_x ,I_y, I_z = rand(1:NX_t),rand(1:NY_t) , rand(1:NZ_t) 
        else 
            println("Erreur : pas de particules")
            println("NX_t :",NX_t)
            println("NY_t :",NY_t)
            println("NZ_t :",NZ_t)
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

        ### mort et naissance de particule Z : 

        # Nombre de d'élements de times qui tombent entre t et t+tau :
        n_rain = 0
        for i in 1:length(times)
            if times[i] > t && times[i] < t + tau
                n_rain += 1
            end
        end

        #on ajoute n_rain particules de pluie à la liste des particules Z aléatoirement sur le [dom, dom]
        for i in 1:n_rain
            push!(positions_Zx[end], rand() .* 2*dom .- dom)
            push!(positions_Zy[end], rand() .* 2*dom .- dom)
        end

        #densité de plantes autour de particule I_z dans un rayon rc[3]
        dens_xz = density(positions_Zx[end][I_z],positions_Zx[end][I_z],positions_Xx[end],positions_Xy[end],rc[3])
        b_z=dens_xz/C[3]
        
        if b_z >= 1
            println("Erreur : b > 1, choisissez un meilleur majorant C.")
            break
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
        ### mort et naissance de particule Y :
        
        #densite de plantes autour de particule I_y dans un rayon rG
        dens_xy = density(positions_Yx[end][I_y],positions_Yx[end][I_y],positions_Xx[end],positions_Xy[end],rG)
        b_y= (L+dens_xy)/C[2]
        if b_y >= 1
            println("Erreur : b > 1, choisissez un meilleur majorant C.")
            break
        end
        if Zy <= b_y  # Mort
            # Supprimer la particule I_y

            deleteat!(positions_Yx[end], I_y)
            deleteat!(positions_Yy[end], I_y)
            morts_Y += 1
        end

        ### mort et naissance de particule X :

        #densite d'eau de sous sol (Y) autour de particule I_x dans un rayon rc[1]
        dens_yx = density(positions_Xx[end][I_x],positions_Xx[end][I_x],positions_Yx[end],positions_Yy[end],rG)
        #densite de plantes autour de particule I_x dans un rayon rc[1]
        dens_xx = density(positions_Xx[end][I_x],positions_Xx[end][I_x],positions_Xx[end],positions_Xy[end],rc[1])
        
        a_x=dens_yx/C[1]
        b_x=a_x+(lamb+(dens_yx*(dens_xx/K)))/C[1]
        if b_x >= 1
            println("Erreur : b > 1, choisissez un meilleur majorant C.")
            break
        end

        if Zx <= a_x  # Naissance
            # naissance d'une particule à coté de la particule I_x
            r = rc[1] * sqrt(rand())
            theta = 2 * π * rand()
            x_new = positions_Xx[end][I_x] + r * cos(theta)
            y_new = positions_Xy[end][I_x] + r * sin(theta)
            push!(positions_Xx[end], x_new)
            push!(positions_Xy[end], y_new)
            naissances_X += 1


        elseif Zx > a_x && Zx <= b_x  # Mort
            # Supprimer la particule I_x
            deleteat!(positions_Xx[end], I_x)
            deleteat!(positions_Xy[end], I_x)
            morts_X += 1
        end


        t += tau  # Mise à jour du temps

    end

    # Retour des résultats
    return positions_Xx, positions_Xy, positions_Yx, positions_Yy, positions_Zx, positions_Zy, niter, naissances_X, morts_X, naissances_Y, morts_Y, naissances_Z, morts_Z   
end