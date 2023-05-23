using LinearAlgebra
using StatsBase
using Random
using JLD2
 

""" create ABM initial conditions """
function ABMinit((p,q))
    (; n, L, M) = q
    (; domain) = p

    @assert M == 2 
    @assert L == 4

    # individuals' opinions
    x =  rand(n,2).*(domain[1,2]-domain[1,1]) .+domain[1,1]  
    # media opinions
    media = [-1. -1.; 1. 1.]  
    # individuals following different influencers in the different quadrants
    xI = Any[]
    push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x>0,x[:,2])))
    push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x>0,x[:,2])))
    push!(xI,  intersect(findall(x-> x>0,x[:,1]), findall(x-> x<=0,x[:,2])) )
    push!(xI,  intersect(findall(x-> x<=0,x[:,1]), findall(x-> x<=0,x[:,2]))) 

    # network between individuals and influencers
    FolInfNet=zeros(n, L)
    for i in 1:L
        FolInfNet[xI[i],i] .=1
    end

    # initial opinions of influencers
    inf = zeros(L,2)
    for i in 1:L
        inf[i,:] = sum(x[xI[i],:],dims=1)/size(xI[i],1) 
    end

    # initial medium per individual, either medium -1 or 1
    state = (rand(n).>0.5)*2 .-1 

    # initialization of interaction network between individuals
    IndNet = ones(n,n)  # everyone connected to everyone including self-connections

    return x, media, inf, FolInfNet, state, IndNet
end

"""
    attraction(x, IndNet)

Function that calculates attraction for the vector `x` with parameters...
"""
function attraction(x, IndNet)
    n = size(IndNet, 1)
    force = zeros(n,2)
    for j in 1:n
        Neighb=findall(x->x==1, IndNet[j,:])
        if isempty(Neighb)
            force[j,:]=[0 0]
        else
            fi = [0 0]
            wsum=0
            for i in eachindex(Neighb)
                d = x[Neighb[i],:]-x[j,:]
                w = exp(-sqrt(d[1]^2 +d[2]^2))
                fi = fi + w*d'
                wsum = wsum+w
            end
            force[j,:] = fi/wsum
        end
    end
    return force
end

function influence(x,media,inf,FolInfNet,state,(p,q))
    (; n, b, c, L) = q
    force1 =zeros(size(x))
    force2 =zeros(size(x))
    for j in 1:n
        if state[j]==1
            force1[j,:] = media[2,:]-x[j,:]
        else  
            force1[j,:] = media[1,:]-x[j,:]
        end

        for k in 1:L
            if FolInfNet[j,k]==1
                force2[j,:]=inf[k,:] -x[j,:]
            end
        end
    end
    force = c*force1 + b*force2
    return force
end

function changeinfluencer(state,x,FolInfNet,inf,(p,q))
    (; eta, n, L) =q
    dt = p.dt
    
    theta =0.1 # threshold for r-function
    
    # compute structural similiarity = fraction of followers with same state
    fraction = zeros(L)
    for i=1:L
        fraction[i]= sum(FolInfNet[:,i].*state)/sum(FolInfNet[:,i])
    end
    
    # compute distance of followers to influencers
    dist = zeros(n, L)
    for i=1:L
        for j=1:n
            d = x[j,:]-inf[i,:]
            dist[j,i]= exp(-sqrt(d[1]^2+d[2]^2))
        end
    end
    
    # compute attractiveness of influencer for followers
    attractive = zeros(n, L)
    for j=1:n
        for i=1:L
            g2 = state[j]*fraction[i]
            if g2<theta 
                g2=theta
            end
            attractive[j,i]= eta * dist[j,i]*g2 
        end

        r=rand()
        lambda = sum(attractive[j,:])
        if r<1-exp(-lambda*dt) 
            p = attractive[j,:]/lambda
            r2=rand()
            k=1
            while sum(p[1:k])<r2
                k=k+1
            end
            FolInfNet[j,:]=zeros(L)      
            FolInfNet[j,k]=1
        end
    end

    return FolInfNet
end

function boundaryconditions(x,domain)
    # reflective boundary conditions
    indx1 = findall(x->x>domain[1,2],x[:,1])
    indy1 = findall(x->x>domain[2,2],x[:,2])
    indx2 = findall(x->x<domain[1,1],x[:,1])
    indy2 = findall(x->x<domain[2,1],x[:,2])
    x[indx1,1] = - x[indx1,1] .+ 2* domain[1,2] 
    x[indy1,2] = - x[indy1,2] .+ 2* domain[2,2] 
    x[indx2,1] = - x[indx2,1] .+ 2* domain[1,1] 
    x[indy2,2] = - x[indy2,2] .+ 2* domain[2,1] 
    return x
end


"""Simulate the ABM """
function ABMsolve(NT = 100;  p = ABMconstruct(), q=parameters(), init="4inf",chosenseed=0)
    Random.seed!(chosenseed)
    (;dt, domain) = p
    (;n, M, L, sigma, sigmahat, sigmatilde, a,  frictionI, frictionM) = q

    x, media, inf, FolInfNet, state, IndNet  = ABMinit((p,q))

    xs = [copy(x)] # list of opinions of individuals in time
    infs = [copy(inf)] # list of opinions of influencers in time
    meds = [copy(media)] # list of opinions of media in time
    stateinfs = [copy(FolInfNet) * collect(1:L)] # list of which influencer an individual follows in time

    for k in 2:NT+1
        xold = x

        # opinions change due to interaction with opinions of friends, influencers and media
        individualforce = attraction(xold,IndNet)
        leaderforce = influence(xold,media,inf,FolInfNet,state,(p,q))
        totalforce = a * individualforce + leaderforce
        x = xold + dt*totalforce + sqrt(dt)*sigma*randn(n,2); 

        infold = inf
        # influencer opinions adapt slowly to opinions of followers with friction
        masscenter=zeros(L,2)
        for i in 1:L
            if sum(FolInfNet[:,i])>0 
                masscenter[i,:] =sum(FolInfNet[:,i] .* xold, dims = 1) /sum(FolInfNet[:,i])    
                inf[i,:] =  inf[i,:] + dt/frictionI * (masscenter[i,:]-inf[i,:]) + 1/frictionI*sqrt(dt)*sigmahat*randn(2,1) 
            else
                inf[i,:] =  inf[i,:] + 1/frictionI*sqrt(dt)*sigmahat*randn(2,1) 
            end
        end
        
        # media opinions change very slowly based on opinions of followers with friction
        masscenter=zeros(M,2)
        states = [-1, 1]
        for i in 1:M
            xM = findall(x->x==states[i], state)
            if size(xM,1)>0
                masscenter[i,:] = sum(xold[xM,:], dims=1)/size(xM,1)
                media[i,:] = media[i,:]  + dt/frictionM * (masscenter[i,:] -media[i,:]) + 1/frictionM * sqrt(dt)*sigmatilde*randn(2,1)
            else
                media[i,:] = media[i,:]  + 1/frictionM * sqrt(dt)*sigmatilde*randn(2,1)
            end
        end
        
        # apply reflective boundary conditions
        x = boundaryconditions(x,domain)
        inf = boundaryconditions(inf,domain)
        media = boundaryconditions(media,domain)

        # individual may jump from one influencer to another
        # jumps according to rate model
        FolInfNet = changeinfluencer(state,xold,FolInfNet,infold,(p,q))

        xs = push!(xs,copy(x))
        infs = push!(infs, copy(inf))
        meds = push!(meds, copy(media))
        stateinfs = push!(stateinfs, copy(FolInfNet * collect(1:L)))
    end

    return xs, stateinfs, infs, meds, state, (p,q)

end

function ABMsolveplot(;NT = 200, ts = [1 11 51 101 151 201],  p = ABMconstruct(), q=parameters(), init="4inf", save=true,seed=0)
    @time xs, stateinfs, infs, meds, state, (p,q) = ABMsolve(NT;  p=p, q=q, init=init,chosenseed=seed)
    ABMplotsnapshots(xs, stateinfs, infs, meds, state, (p,q), ts; name=init,save=save)
    ABMplotfollowernumbers(stateinfs,state,(p,q), save=save)
    return xs, stateinfs, infs, meds, state, (p,q)
end
