using DifferentialEquations, LinearAlgebra
using RecursiveArrayTools
using Plots

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
    @inbounds for i=1:I #current attitude
        for j in 1:J #current influencer
            for j2 in 1:J #future influencer
                if j2!=j
                    for x in 1:Nx*Ny
                        rate[x,i,j,j2] = exp(-norm(grid_points[x,:]-y[:,j2]))*g((m[i,j2]-m[mod1(i+1,Int(I)),j2])/(m[i,j2]+m[mod1(i+1,Int(I)),j2]))
                    end
                end
            end
        end
    end
    return eta*reshape(rate,(Nx, Ny, I, J, J))
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

function construct()
    # Define the constants for the PDE
    sigma = 0.5
    D = sigma^2 * 0.5
    a = 1.0
    b = 2.0
    c = 2.0
    eta = 50.0
    Gamma_0 = 100
    gamma_0 = 25
    J=4 # number of influencers
    dx = 0.05
    dy = dx
    dV = dx*dy
    domain = [-2 2; -2 2]
    N_x = Int((domain[1,2]-domain[1,1])/dx+1)
    N_y = N_x #so far only works if N_y = N_x
    N = N_x*N_y
    X = [x for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    Y = [y for x in domain[1,1]:dx:domain[1,2], y in domain[2,1]:dy:domain[2,2]]
    grid_points = [vec(X) vec(Y)] 
    # matrix of K evaluated at gridpoints
    K_matrix, W_matrix  = generate_K(grid_points)

    M = second_derivative((N_x,N_y), (dx, dy))
    C = centered_difference((N_x,N_y), (dx, dy))
 
    p = (; grid_points, N_x, N_y, domain,  a, b, c, eta, K_matrix, W_matrix, dx,dy, dV, C,  D, M , N, J, Gamma_0, gamma_0)

    return p
end

function initialconditions(p)
    (; N_x , N_y,  J) = p
    rho_0 = zeros(N_x, N_y, 2, 4) 
    mid_y =Int(round(N_y/2))
    mid_x =Int(round(N_x/2))
    rho_0[1:mid_x, 1:mid_y,:,4] .= 1/(2*J*4)
    rho_0[1:mid_x, mid_y+2:end,:,2] .= 1/(2*J*4)
    rho_0[mid_x+2:end, 1:mid_y,:,3] .= 1/(2*J*4)
    rho_0[mid_x+2:end, mid_y+2:end,:,1] .= 1/(2*J*4)
    #cat(1/32*ones(N_x,N_y), 1/32*ones(N_x,N_y), dims=3)
    #rho_0 = fill(1/(2*J*4*4), N_x, N_y, 2, J)
    u0 = rho_0
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y4_0 = [-1.,-1.]
    y2_0 = [-1.,1.]
    y3_0 = [1.,-1.]
    y1_0 = [1.,1. ]

    y0 = [y1_0  y2_0 y3_0 y4_0]
    return ArrayPartition(u0,z0,y0)
end

#gaussian that integrates to 1 and centered at center
gaussian(x,center, sigma=0.15) = 1/(2*pi*sigma^2) * exp(-1/(2*sigma^2)*norm(x-center)^2)

function random_initialconditions(p)
    
    (; grid_points, N_x , N_y,  dV, J) = p
    random_pos = rand(128,2).*4 .-2 #128 uniform samples in domain
    rho_0 = zeros(N_x, N_y, 2, 4) 
    counts = [0 0 0 0]
    y0 = zeros(2,4)
    for i in 1:128
        if random_pos[i,2]>0
            if random_pos[i,1]>0
                rho_0[:,:,rand([1,2]),1]+= reshape([gaussian(grid_points[j,:], random_pos[i,:]) for j in 1:N_x*N_y], N_x, N_y)
                counts[1]+=1
                y0[:,1]+=random_pos[i,:]
            else
                rho_0[:,:,rand([1,2]),2]+= reshape([gaussian(grid_points[j,:], random_pos[i,:]) for j in 1:N_x*N_y], N_x, N_y)
                counts[2]+=1
                y0[:,2]+=random_pos[i,:]
            end
        else
            if random_pos[i,1]>0
                rho_0[:,:, rand([1,2]),3]+= reshape([gaussian(grid_points[j,:], random_pos[i,:]) for j in 1:N_x*N_y], N_x, N_y)
                counts[3]+=1
                y0[:,3]+=random_pos[i,:]
            else
                rho_0[:,:, rand([1,2]), 4]+= reshape([gaussian(grid_points[j,:], random_pos[i,:]) for j in 1:N_x*N_y], N_x, N_y)
                counts[4]+=1
                y0[:,4]+=random_pos[i,:]
            end
        end
    end
    rho_0 = rho_0/(sum(rho_0)*dV)
    u0 = rho_0
    z2_0 = [1.,1.]
    z1_0 = [-1.,-1.]
    z0 = [z1_0  z2_0]
    y0= y0./counts
    return ArrayPartition(u0,z0,y0)
end

function f(duzy,uzy,p,t)
    yield()
    (; grid_points, N_x, N_y, domain,  a, b, c, eta, K_matrix, W_matrix, dx,dy, dV, C,  D, M , N, J, Gamma_0, gamma_0) = p
    u, z, y2 = uzy.x
    du, dz, dy2 = duzy.x

    rhosum = sum(u, dims=(3,4))[:,:,1,1]
    rhosum_j = sum(u, dims=3)
    rhosum_i = sum(u, dims=4)
    m_i = dV*sum(u,dims = (1,2,4))
    m_j = dV*sum(u,dims=(1,2,3))
    Fagent = agent_force(rhosum, K_matrix, W_matrix, dV)
    rate_matrix = gamma(grid_points, u, y2, eta,dV)
    for i in 1:2
        for j in 1:J
            rho = @view  u[:,:,i,j]
            drho = @view du[:,:,i,j]
            zi = @view z[:,i]
            dzi = @view dz[:,i]
            yj = @view y2[:,j]
            dyj = @view dy2[:,j]

            force = c * follower_force(zi, grid_points, N_x, N_y) + a * Fagent + b * follower_force(yj, grid_points, N_x, N_y)            
            
            div =  C * (rho .* force[:,:,1]) + (rho .* force[:,:,2]) * C'
            reac = zeros(N_x, N_y)
            for j2=1:J
                if j2!= j
                    reac += -rate_matrix[:,:,i,j,j2] .* rho + rate_matrix[:,:,i,j2,j] .* u[:,:,i,j2]
                end
            end
            
            dif = D*(M*rho + rho*M')  
            #balance fluxes at boundary (part of boundady conditions)
            dif[1,:]+= -D/dx * (force[1,:,1].*rho[1,:]) 
            dif[end,:]+= D/dx * (force[end,:,1].*rho[end,:])
            dif[:,1]+= -D/dy * (force[:,1,2].*rho[:,1])
            dif[:,end]+= D/dy * (force[:,end,2].*rho[:,end])

            drho .=  dif - div + reac

            mean_rhoi = 1/m_i[i] * dV*reshape(rhosum_i[:,:,i,:],1,N)*grid_points
            dzi .= 1/(Gamma_0*m_i[i]) * (mean_rhoi' - zi)

            mean_rhoj = 1/m_j[j] * dV*reshape(rhosum_j[:,:,:,j],1,N)*grid_points
            dyj .= 1/(gamma_0*m_j[j]) * (mean_rhoj' - yj)
        end
    end

end

function solve(tmax=0.1; alg=nothing, init="random", p = construct())
    
    if init=="random"
        uzy0 = random_initialconditions(p)
    else
        uzy0 = initialconditions(p)
    end
    
    # Solve the ODE
    prob = ODEProblem(f,uzy0,(0.0,tmax),p)
    @time sol = DifferentialEquations.solve(prob, alg, progress=true,save_everystep=true,save_start=true)
    
    return sol, p
end

function ensemble(tmax=0.1, N=10; alg=nothing, init="random")
    zs = Any[]
    ys = Any[]
    us = Any[]
    p = construct()
    for i=1:N
        sol,p = solve(tmax; alg=alg, init=init, p=p)
        push!(us, sol(tmax).x[1])
        push!(zs, sol(tmax).x[2])
        push!(ys, sol(tmax).x[3])
    end
    return us, zs, ys
end

function plot_ensemble(us, zs, ys, p)
    ys = cat(ys..., dims=3)
    N=size(ys)[3]
    zs = cat(zs..., dims=3)
    us = cat(us...,dims=5)
    av_u = sum(us,dims=5)*(1/N)

    (; domain, dx, dy,J) = p

    #PLOTTING
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    
    plot_array = Any[]  
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"]
    for j in 1:J    
        for i in 1:2
            # make a plot and add it to the plot_array
            push!(plot_array, plot_solution2(av_u[:,:,i,j], zs[:,i,:], ys[:,j,:], x_arr, y_arr; title = string(dens_labels[i,j]) ,labely = y_labels[j], labelz = z_labels[i]))
 
        end
    end
    plot(plot_array..., layout=(4,2),size=(1000,1000)) |> display
    savefig("img/ensemble_average_agents.png")
    
    y_all = reshape(ys, (2, N*J))
    plot(histogram2d(y_all[1,:], y_all[2,:], bins=(-2:0.1:2, -2:0.1:2)), title="Final distribution of influencers") |> display
    savefig("img/influencer_dist.png")
end


function plot_solution(rho, z, y, x_arr, y_arr; title="", labelz="", labely="", clim=(-Inf, Inf))
    subp = heatmap(x_arr,y_arr, rho', title = title, c=:berlin, clims=clim)
    scatter!(subp, [z[1]], [z[2]], markercolor=[:yellow],markersize=6, lab=labelz)
    scatter!(subp, [y[1]], [y[2]], markercolor=[:red],markersize=6, lab=labely)
    return subp
end

function plot_solution2(rho, z, y, x_arr, y_arr; title="", labelz="", labely="", clim=(-Inf, Inf))
    subp = heatmap(x_arr,y_arr, rho', title = title, c=:berlin, clims=clim)
    scatter!(subp, z[1,:], z[2,:], markercolor=[:yellow],markersize=6, lab=labelz)
    scatter!(subp, y[1,:], y[2,:], markercolor=[:red],markersize=6, lab=labely)
    return subp
end

function solveplot(tmax=0.1; alg=nothing)
    sol, p = solve(tmax; alg=nothing)


    (; domain, dx, dy,J) = p
    #PLOTTING
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    
    plot_array = Any[]  
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"]
    for j in 1:J    
        for i in 1:2
            # make a plot and add it to the plot_array
            push!(plot_array, plot_solution(sol(tmax).x[1][:,:,i,j], sol(tmax).x[2][:,i], sol(tmax).x[3][:,j], x_arr, y_arr; title = string(dens_labels[i,j],"(", string(tmax), ")") ,labely = y_labels[j], labelz = z_labels[i]))
        end
    end
    plot(plot_array..., layout=(4,2),size=(1000,1000)) |> display

    #p2 = plot_solution(sol(tmax).x[1][:,:,2], sol(tmax).x[2][:,2], x_arr, y_arr; title=string("ρ₋₁(",string(tmax),")"), label="z₋₁")
    #plot(p1, p2, layout=[4 4], size=(1000,4*400)) 
    savefig("img/finaltime_pde.png")
    return sol, p
end
 
function creategif(sol,p, dt=0.01)
    #rho1(t)=sol(t).x[1][:,:,1]
    #'rho2(t)=sol(t).x[1][:,:,2]
    #z1(t)=sol(t).x[2][:,1]
    #z2(t)=sol(t).x[2][:,2]
    
    tmax=sol.t[end]
    (; domain, dx, dy, J) = p
    x_arr = domain[1,1]:dx:domain[1,2]
    y_arr = domain[2,1]:dy:domain[2,2]
    cl = (0, maximum(maximum(sol(t).x[1]) for t in 0:dt:tmax)) #limits colorbar
    z_labels = ["z₋₁","z₁" ]
    y_labels = ["y₁", "y₂", "y₃", "y₄"]
    dens_labels = [  "ρ₋₁₁" "ρ₋₁₂" "ρ₋₁₃" "ρ₋₁₄";"ρ₁₁" "ρ₁₂" "ρ₁₃" "ρ₁₄"]

    pdegif = @animate for t = 0:dt:tmax
        plot_array = Any[]  
        for j in 1:J    
            for i in 1:2
                # make a plot and add it to the plot_array
                push!(plot_array, plot_solution(sol(t).x[1][:,:,i,j], sol(t).x[2][:,i], sol(t).x[3][:,j], x_arr, y_arr; title = string(dens_labels[i,j],"(", string(t), ")") ,labely = y_labels[j], labelz = z_labels[i], clim = cl))
            end
        end
        plot(plot_array..., layout=(4,2),size=(1000,1000)) |> display
    end
    Plots.gif(pdegif, "img/evolution.gif", fps = 30)
end

function test_f()
    p = construct()
    uzy0 = initialconditions(p)
    duzy = copy(uzy0)
    @time f(duzy, uzy0, p, 0)
    return duzy
end
