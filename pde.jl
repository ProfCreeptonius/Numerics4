using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using Plots
using JLD2
using Distances

function generate_K(grid_points)
    N, d = size(grid_points)
    @assert d == 2
    k = zeros(N, N, 2)
    w = zeros(N, N)

    @inbounds for j in 1:N, i in 1:N
        x = grid_points[j,1] - grid_points[i,1] #first component of x' -x
        y = grid_points[j,2] - grid_points[i,2] #2nd component of x' -x
        w[i,j] = ww = exp(-sqrt(x^2 + y^2))
        k[i,j, 1] = ww * x
        k[i,j, 2] = ww * y
    end
    return k, w
end

function g(x)
    if x>0
        return 0.1+x
    else
        return 0.1
    end
end

function gamma(grid_points, rho, y, eta,dV)
    Nx, Ny, I, J = size(rho)
    rate = zeros((Nx*Ny, I, J, J))
    m = dV*sum(rho, dims=(1,2))[1,1,:,:]
    dists = pairwise(Euclidean(), grid_points', y)

    for i in 1:I, j in 1:J, j2 in 1:J
        j == j2 && continue
        @. rate[:,i,j,j2] = eta * exp.(-dists[:, j2]) * g((m[i,j2]-m[mod1(i+1,Int(I)),j2])/(m[i,j2]+m[mod1(i+1,Int(I)),j2]))
    end

    return reshape(rate,(Nx, Ny, I, J, J))
end

function agent_force(rho, K_matrix, W_matrix,  dV)
    force = zeros(size(rho)..., 2)
    @views for d in 1:2
        f = vec(force[:,:,d])
        f .= dV * K_matrix[:,:,d] * vec(rho)
    end

    norms = dV * W_matrix * vec(rho)

    return force./reshape(norms, size(rho))
end

function follower_force(z, grid_points, N_x, N_y)
    pointwise_int =   z' .- grid_points
    force = reshape(pointwise_int,N_x,N_y, 2)
    return force
end


function second_derivative((N_x,N_y), (dx, dy))
    # here matrix should have different shape for Ny not equal Nx
    M = Tridiagonal(ones(N_x-1), fill(-2., N_x), ones(N_x-1))
    # boundary conditions change in time to ensure flux balance, they are added when solving pde
    M[1,1] = -1.
    M[end,end] = -1.
    M .= (1/dx^2)*M

    return M
end

function centered_difference((N_x,N_y), (dx, dy))
    #centered first difference for force, doesnt work for different x, y grids
    C = 1/(2*dx)*Tridiagonal(-ones(N_x-1), zeros(N_x), ones(N_x-1))
    # at the boundary a one-sided difference scheme is used
    C[1,1:2] = 1/(dx)* [-1,1]
    C[end,end-1:end] = 1/(dx)* [-1,1]

    return C
end


function parameters(;J=4, b=2.5, eta=15.0, controlspeed = 0.25, frictionI = 2) #a=3 makes interesting case too
    a = 1. #a=1 in paper, interaction strength between agents
    #b = 2. # interaction strength between agents and influencers
    c = 1. # interaction strength between agents and media
    #eta = 15.0 #rate constant for changing influencer
    n = 250 #number of agents, important for random initial conditions
    #J = 4 #number of influencers
    n_media = 2
    sigma = 0.5 # noise on individual agents
    sigmahat = 0 # noise on influencers
    sigmatilde = 0 # noise on media
     # friction for influencers
    frictionM = 100  #friction for medi

    q = (; n, J, n_media, frictionM, frictionI, a, b, c, eta, sigma, sigmahat, sigmatilde, controlspeed)
    return q
end

function PDEconstruct()
    # Define the constants for the PDE
    dx = 0.1
    dy = dx

    domain = [-2.5 2.5; -2.5 2.5]
    N_x = Int((domain[1,2]-domain[1,1])/dx+1)
    N_y = Int((domain[2,2]-domain[2,1])/dy+1) #so far only works if N_y = N_x
    N = N_x*N_y
    dV = dx*dy # to integrate the agent distribution
    X = [x for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    Y = [y for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    grid_points = [vec(X) vec(Y)]
    # matrix of K evaluated at gridpoints
    K_matrix, W_matrix  = generate_K(grid_points)

    M = second_derivative((N_x,N_y), (dx, dy))
    C = centered_difference((N_x,N_y), (dx, dy))

    p = (; grid_points, dx, dy, dV, X, Y, N, N_x, N_y, domain, K_matrix, W_matrix, C,  M)

    return p
end

function initialconditions(P)
    (; N_x , N_y,  J, n, dV, domain, dx) = P
    rho_0 = zeros(N_x, N_y, 2, 4)
    mid_y =Int(round(N_y/2))
    mid_x =Int(round(N_x/2))
    start_x = Int((domain[1,2] - 2)/dx + 1)
    end_x = N_x - start_x +1
    rho_0[start_x:mid_x, start_x:mid_y,:,4] .= 1
    rho_0[start_x:mid_x, mid_y+2:end_x,:,2] .= 1
    rho_0[mid_x+2:end_x, start_x:mid_y,:,3] .= 1
    rho_0[mid_x+2:end_x, mid_y+2:end_x,:,1] .= 1
    rho_0[mid_x+1, mid_x+1,:,:] .= 0.5

    u0 = rho_0/(sum(rho_0)*dV)
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y4_0 = [-1.,-1.]
    y2_0 = [-1.,1.]
    y3_0 = [1.,-1.]
    y1_0 = [1.,1. ]

    y0 = [y1_0  y2_0 y3_0 y4_0]
    counts = n/(2*J)*ones(J,2) # proportion of agents that follow each influencer
    controlled = zeros(J)
    return ArrayPartition(u0,z0,y0), counts, controlled
end

#gaussian that integrates to 1 and centered at center
gaussian(x, center, sigma=0.2) = 1/(2*pi*sigma^2) * exp(-1/(2*sigma^2)*norm(x-center)^2)

function inf_initialconditions(P)

    (; grid_points, N_x , N_y,  dV, J,n) = P
    random_pos = rand(n,2).*4 .-2 #n uniform samples in domain
    rho_0 = zeros(N_x, N_y, 2, J)
    counts = zeros(J,2)
    y0 = zeros(2,J)

    function add_agent(agenti, infi, state,  rho_0, counts, y0)
        rho_0[:,:,state,infi]+= reshape([gaussian(grid_points[j,:], random_pos[agenti,:]) for j in 1:N_x*N_y], N_x, N_y)
        counts[infi, state]+=1
        y0[:,infi]+=random_pos[agenti,:]
        return rho_0, counts, y0
    end

    for i in 1:n
        state = rand([1,2])
        if random_pos[i,2]>0
            if random_pos[i,1]>0
                rho_0, counts, y0 = add_agent(i, 1, state,  rho_0, counts, y0)
            else
                rho_0, counts, y0 = add_agent(i, 2, state,  rho_0, counts, y0)
            end
        else
            if random_pos[i,1]>0
                rho_0, counts, y0 = add_agent(i, 3, state,  rho_0, counts, y0)
            else
                rho_0, counts, y0 = add_agent(i, 4, state,  rho_0, counts, y0)
            end
        end
    end

    rho_0 = rho_0/(sum(rho_0)*dV)
    u0 = rho_0
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y0= y0./dropdims(sum(counts, dims=2), dims=2)'
    controlled = zeros(J)
    return ArrayPartition(u0,z0,y0), counts, controlled
end

function noinf_initialconditions(P)

    (; grid_points, N_x , N_y,  dV, J ,n) = P
    random_pos = rand(n,2).*4 .-2 #n uniform samples in domain
    rho_0 = zeros(N_x, N_y, 2, 1)
    counts = zeros(1,2)
    y0 = zeros(2,1)

    function add_agent(agenti, infi, state,  rho_0, counts, y0)
        rho_0[:,:,state,infi]+= reshape([gaussian(grid_points[j,:], random_pos[agenti,:]) for j in 1:N_x*N_y], N_x, N_y)
        counts[infi, state]+=1
        y0[:,infi]+=random_pos[agenti,:]
        return rho_0, counts, y0
    end

    for i in 1:n
        state = rand([1,2])

        rho_0, counts, y0 = add_agent(i, 1, state,  rho_0, counts, y0)

    end
    rho_0 = rho_0/(sum(rho_0)*dV)
    u0 = rho_0
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    controlled = zeros(1)
    return ArrayPartition(u0,z0,y0), counts, controlled
end

function constructinitial(scenario,P)
    if scenario=="4inf"
        uzy0, counts, controlled = inf_initialconditions(P)
    elseif scenario=="noinf"
        uzy0, counts, controlled = noinf_initialconditions(P)
    elseif scenario =="uniform"
        uzy0, counts, controlled = initialconditions(P)
    elseif scenario=="controlled"
            uzy0, counts, controlled = inf_initialconditions(P)
    end
    return uzy0, counts, controlled
end


function f(duzy,uzy,P,t)
    yield()
    (; grid_points, N_x, N_y,  a, b, c, sigma, eta, K_matrix, W_matrix, dx,dy, dV, C,  M , N, J, frictionM, frictionI, controlled, controlspeed) = P
    D = sigma^2 * 0.5
    u, z, y2 = uzy.x
    du, dz, dy2 = duzy.x


    rhosum = sum(u, dims=(3,4))[:,:,1,1]
    rhosum_j = sum(u, dims=3)
    rhosum_i = sum(u, dims=4)
    m_i = dV*sum(u,dims = (1,2,4))
    m_j = dV*sum(u,dims=(1,2,3))
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    rate_matrix = gamma(grid_points, u, y2, eta,dV)
    reac = zeros(N_x, N_y)
    for i in 1:2
        for j in 1:J
            rho = @view  u[:,:,i,j]
            drho = @view du[:,:,i,j]
            zi = @view z[:,i]
            dzi = @view dz[:,i]
            yj = @view y2[:,j]
            dyj = @view dy2[:,j]

            force = c * follower_force(zi, grid_points, N_x, N_y) + a * Fagent + b * follower_force(yj, grid_points, N_x, N_y)

            dive =  C * (rho .* force[:,:,1]) + (rho .* force[:,:,2]) * C'

            reac .= 0
            for j2=1:J
                if j2!= j
                    @. @views reac += -rate_matrix[:,:,i,j,j2] .* rho + rate_matrix[:,:,i,j2,j] .* u[:,:,i,j2]
                end
            end

            dif = D*(M*rho + rho*M')
            #balance fluxes at boundary (part of boundady conditions)
            dif[1,:]+= -D/dx * (force[1,:,1].*rho[1,:])
            dif[end,:]+= D/dx * (force[end,:,1].*rho[end,:])
            dif[:,1]+= -D/dy * (force[:,1,2].*rho[:,1])
            dif[:,end]+= D/dy * (force[:,end,2].*rho[:,end])

            drho .=  dif - dive + reac

            mean_rhoi = 1/m_i[i] * dV*reshape(rhosum_i[:,:,i,:],1,N)*grid_points
            dzi .= 1/(frictionM) * (mean_rhoi' - zi)

            if controlled[j] == 0
                mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
                dyj .= 1/(frictionI) * (mean_rhoj' - yj)
            else #controll movement
                dyj .= controlspeed* ([1.5 1.5]' - yj)
            end

        end
    end


end



function sol2uyz(sol, t)
    u = sol(t).x[1]
    z = sol(t).x[2]
    y = sol(t).x[3]
    return u,z,y
end


function solve(tmax=0.1; alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())

    P = (; scenario, p..., q...)
    uzy0, counts, controlled = constructinitial(scenario,P)
    P = (; P..., controlled)

    # Solve the ODE
    prob = ODEProblem(f,uzy0,(0.0,tmax),P)
    @time sol = DifferentialEquations.solve(prob, alg, save_start=true)
    
    return sol, P, counts
end

function solvecontrolled(tcontrol = 0.05, tmax=0.1; alg=nothing, scenario="controlled", p = PDEconstruct(), q= parameters(), savedt=0.05, atol = 1e-6, rtol = 1e-3)
    
    P1 = (; scenario, p..., q...)
    uzy0, counts, controlled = constructinitial(scenario,P1)
    P1 = (; P1..., controlled)

    # Solve the ODE
    prob1 = ODEProblem(f,uzy0,(0.0,tcontrol),P1)

    @time sol1 = DifferentialEquations.solve(prob1, alg, saveat = 0:savedt:tcontrol,save_start=true, abstol = atol, reltol = rtol)
    

    #add new influencer
    u,z,y = sol2uyz(sol1,tcontrol)
    P2 = merge(P1, (;J=P1.J+1,controlled = [P1.controlled..., 1] ))
    (;N_x, N_y,J, grid_points,N) = P2
    u2 = zeros(N_x, N_y, 2, J)
    u2[:,:,:,1:J-1] = u
    y2 = zeros(2, J)
    y2[:,1:J-1] = y
    y2[:,J] =  1/sum(u2,dims=(1,2,4))[1,1,2,1] * reshape(sum(u2, dims=4)[:,:,2,:],1,N)*grid_points
    uzy0 = ArrayPartition(u2,z,y2)

    # solve ODE with added influencer
    prob2 = ODEProblem(f,uzy0,(0.0,tmax-tcontrol),P2)
    @time sol2 = DifferentialEquations.solve(prob2, alg,  saveat = 0:savedt:(tmax-tcontrol),save_start=true, abstol = atol, reltol = rtol)

    return [sol1, sol2], [P1, P2], counts
end

function solveplot(tmax=0.1; alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())
    sol,P, counts = solve(tmax; alg=alg, scenario=scenario, p=p, q=q)

    u,z,y = sol2uyz(sol, tmax)

    plotarray(u,z,y, P, tmax)

    return sol, P, counts
end


function solveensemble(tmax=0.1, N=10; savepoints = 4, alg=nothing, scenario="4inf", p = PDEconstruct(), q= parameters())
    P = (; p..., q...)
    (; N_x, N_y,J) = P

    zs = zeros(2, 2, savepoints, N)
    ys = zeros(2, J, savepoints, N)
    us = zeros(N_x, N_y, 2, J, savepoints,  N)
    savetimes = LinRange(0, tmax, savepoints)
    av_counts = zeros(J,2)
    Threads.@threads for i=1:N
        sol, _ , counts= solve(tmax; alg=alg, scenario=scenario, p=p, q=q)
        av_counts = av_counts +  counts*(1/N)
        for j in 1:savepoints
            u,z,y = sol2uyz(sol, savetimes[j])
            us[:,:,:,:,j,i] = u
            zs[:,:,j,i] = z
            ys[:,:,j,i] = y
        end

    end

    @save string("data/pde_ensemble_",scenario,".jld2") us zs ys
    return us, zs, ys, P, av_counts
end

function plotensemble(us, zs, ys, P, tmax; title1 = "img/pde_ensemble", title2 = "img/pde_ensemble_influencer_", clmax = 0.5, scenario="4inf")
    N = size(ys,4)
    savepoints = size(us,5)
    av_u = sum(us,dims=6)*(1/N)
    av_z = sum(zs,dims=4 )*(1/N)
    av_y = sum(ys,dims=4 )*(1/N)
    savetimes = LinRange(0,tmax,savepoints)

    for k in 1:savepoints
        plotarray(av_u[:,:,:,:,k], av_z[:,:,k], av_y[:,:,k], P, savetimes[k]; save=false, clmax = clmax, scenario = scenario)
        savefig(string(title1,string(k),scenario,".png"))

        if scenario=="4inf"
            (;X, Y, domain, dx, dy, J) = P
            x_arr = domain[1,1]:dx:domain[1,2]
            y_arr = domain[2,1]:dy:domain[2,2]

            yall = reshape(ys[:,:,k,:], (2,N*J))'
            evalkde = [sumgaussian([X[i,j], Y[i,j]], yall) for i in 1:size(X,1), j in 1:size(X,2)]
            heatmap(x_arr, y_arr, evalkde', c=:berlin, title=string("Distribution of influencers at time ", string(round(savetimes[k], digits=2)))) |> display
            savefig(string(title2,string(k),".png"))
        end
    end


end

function psensemble(tmax=0.1, N=10; alg=nothing, scenario="4inf")
    us, zs, ys, P, av_counts  = solveensemble(tmax, N; alg=alg, scenario=scenario)
    plotensemble(us, zs, ys, P, tmax)
    return us, zs, ys, P, av_counts
end

function plot_solution(rho, z, y, x_arr, y_arr; title="", labelz="", labely="", clim=(-Inf, Inf), scenario="4inf")
    subp = heatmap(x_arr,y_arr, rho', title = title, c=:berlin, clims=clim)
    scatter!(subp, z[1,:], z[2,:], markercolor=:yellow,markersize=6, lab=labelz)
    if scenario=="4inf"
        scatter!(subp, y[1,:], y[2,:], markercolor=:red,markersize=6, lab=labely)
    end
    return subp
end


function plotarray(u,z,y, P, t; save=true, clmax = maximum(u), scenario="4inf")
    (; domain, dx, dy, dV, J) = P

    #u,z,y = sol2uyz(sol, t)

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]

    array = Any[]
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄", "y₅"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄" "ρ₋₁₅";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄" "ρ₁₅"]
    for j in 1:J
        for i in 1:2
            # make a plot and add it to the array
            cl = (0, clmax) #limits colorbar
            title = string(dens_labels[i,j],"(", string(round(t, digits=2)), "), prop = ", string(round(sum(u[:,:,i,j]*dV), digits = 3)))
            push!(array, plot_solution(u[:,:,i,j], z[:,i], y[:,j], x_arr, y_arr; title = title,labely = y_labels[j], labelz = z_labels[i], clim=cl, scenario=scenario))
        end
    end
    plt = plot(array..., layout=(J,2),size=(1000,min(J*250,1000)))
    plt |> display

    if save==true
        savefig(string("img/pde_array_",scenario,".png"))
    end

    return plt
end

#=
function gifarray

    tmax=sol.t[end]
    #cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar

    pdegif = @animate for t = 0:dt:tmax
        u,z,y = sol2uyz(sol, t)
        plotarray(u,z,y, P, t; save=false)
    end
    Plots.gif(pdegif, string("img/pde_array_",scenario,".gif"), fps = 10)
end
=#

gifarray(sol, P, args...; kwargs...) = gifarray([sol], [P], args...; kwargs...)

function gifarray(sols::Vector, Ps::Vector, dt=0.01; scenario = "4inf")
    T = 0
    anim = Animation()
    for (sol, P) in zip(sols, Ps)
        for t in 0:dt:sol.t[end]
            u,z,y = sol2uyz(sol, t)
            plt = plotarray(u,z,y, P, t+T; save=false)
            frame(anim, plt)
        end
        T += sol.t[end]
    end
    Plots.gif(anim, string("img/pde_array_",scenario,".gif"), fps = 10)
end

function plotsingle(u,z,y,P,t; save=true, scenario="4inf")
    (; domain, dx, dy, J) = P
    #u,z,y = sol2uyz(sol, t)

    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]

    dens = dropdims(sum(u, dims=(3,4)), dims=(3,4))

    subp = heatmap(x_arr,y_arr, dens', title = string("t=", string(round(t, digits=2))), c=:berlin) 

    scatter!(subp, z[1,:], z[2,:], markercolor=:yellow,markersize=4, lab = "media")

    if scenario!="noinf"
        scatter!(subp, y[1,:], y[2,:], markercolor=:red,markersize=4, lab="influencers")
    end

    subp |> display

    if save==true
        savefig(string("img/pde_single_",scenario,".png"))
    end
    return subp
end

gifsingle(sol, P, args...; kwargs...) = gifsingle([sol], [P], args...; kwargs...)

function gifsingle(sols::Vector, Ps::Vector, dt=0.01; scenario = "4inf")
    T = 0
    anim = Animation()
    for (sol, P) in zip(sols, Ps)
        for t in 0:dt:sol.t[end]
            u,z,y = sol2uyz(sol, t)
            plt = plotsingle(u,z,y,P,t+T,save=false, scenario=scenario)
            frame(anim, plt)
        end
        T += sol.t[end]
    end
    Plots.gif(anim, string("img/pde_single_",scenario,".gif"), fps = 10)
end

function test_f()
    p = PDEconstruct()
    q = parameters()
    P = (; p..., q...)
    uzy0 = initialconditions(P)
    duzy = copy(uzy0)
    @time f(duzy, uzy0, P, 0)
    return duzy
end

function solvenoinf(tmax; alg=nothing)
    sol, P, _ = solve(tmax; alg=alg,  scenario="noinf", p = PDEconstruct(), q= parameters(J=1, b=0, eta=0))
    return sol, P
end

function solveplotnoinf(tmax; alg=nothing)
    sol, P = solvenoinf(tmax; alg=alg)
    (;scenario) =P
    u,z,y = sol2uyz(sol, tmax)
    plotsingle(u,z,y,P,tmax; scenario=scenario)
    plotarray(u,z,y, P, tmax;  scenario=scenario)
    return sol, P
end
