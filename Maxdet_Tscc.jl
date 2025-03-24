module Maxdet_Tscc


########### Imports ###########

using LinearAlgebra
using Convex
using SCS

include("./start_subdivision.jl")

using .start_subdivision




########### Functions ###########

"""
    solve()

    Input a set of secondary cones and solve the determinant maximization problem [see. Val03 p80] for every secondary cone.
"""
function solve(cones::Set{}, d::Int64)
    c = collect(cones)
    t = Set{Tuple{BigFloat,Matrix}}()
    dim_cone = Int(d*(d+1)/2)

    for o in 1:length(cones)

        #Setup the programm
        constraint = []
        x = Variable(dim_cone)

        #PSD constraint
        Gx = sum(x[i]*vec2mat(Diagonal(ones(dim_cone))[i,:]) for i = 1:dim_cone)
        objective = -logdet(Gx)
        push!(constraint,Gx ⪰ 0)

        #Containment in T-secondary cone constraint
        halfspaces = collect(keys(c[o].Inequalities))
        Mat_vec_half = Vector{Matrix{Rational{BigInt}}}(undef,dim_cone)
        for j = 1:dim_cone
            m = zeros(length(halfspaces),length(halfspaces))
            for (k,H) in enumerate(halfspaces)
                m[k,k] = -(H.a)[j]
            end
            Mat_vec_half[j] = deepcopy(m)
        end
        push!(constraint,(sum(x[j]*Mat_vec_half[j] for j=1:dim_cone)) ⪰ 0)
        hyperplanes = collect(c[o].Equalities)
        Mat_vec_half = Vector{Matrix{Rational{BigInt}}}(undef,dim_cone)
        for (j,H) in enumerate(hyperplanes)
            push!(constraint,sum(x[i]*H.a[i] for i=1:dim_cone) == 0)
        end

        #Constraint for every delone cell      NOTE: Could be faster by only taking traslational-inequivalent delone cells
        D = collect(c[o].DeloneSub.Cells)
        ml = zeros(d+1,d+1)
        ml[1,1] = 4
        for (i,Cell) in enumerate(D)
            V = start_subdivision.affine_independent_set(Cell.Vertices,d,d+1)
            V = filter!( !=(zeros(d)) , V)
            storevec = zeros(dim_cone,length(V),length(V))
            for (k,vert1) in enumerate(V)
                for (l,vert2) in enumerate(V)
                    storevec[:,k,l] = start_subdivision.mat2vecdouble(vert1*vert2')
                end
            end
            Mat_vec = Vector{Matrix{Rational{BigInt}}}(undef,dim_cone)
            for j = 1:dim_cone
                m1 = zeros(d+1,d+1)
                for k = 1:d+1
                    for l = 1:d+1
                        if k ==1 && l ==1
                            m1[1,1] = 0
                        elseif k == 1 && l != 1
                            m1[1,l] = storevec[j,l-1,l-1]
                        elseif k!= 1 && l ==1
                            m1[k,1] = storevec[j,k-1,k-1]
                        else
                            m1[k,l] = storevec[j,k-1,l-1]
                        end
                    end
                end
                Mat_vec[j] = deepcopy(m1)
            end
            push!(constraint,(ml+sum(x[j]*Mat_vec[j] for j=1:dim_cone)) ⪰ 0)
        end

        #Solve
        p = minimize(objective, constraint)
        solve!(p,SCS.Optimizer; silent = true)
        Q = vec2mat(Convex.evaluate(x))
        density = sqrt((1^d)/det(Q))
        push!(t,(density,copy(Q)))
    end
    t1 = collect(t)
    
    i = argmin(first.(t1))
    return t1[i]
end

"""
    vec2mat(v)

return the symmetric Matrix `M` which is encoded in the vector `v`.

See also [`mat2vec`](@ref)
"""
function vec2mat(v)
    n = Int((sqrt(8*size(v)[1] +1)-1)/2)
    M = zeros(typeof(v[1]),(n,n))
    k = 1
    for i in 1:n
        for j in 1:i-1
            M[i,j] = v[k]
            #M[i,j] = v[k]/2
            M[j,i] = M[i,j]
            k = k+1
        end
        M[i,i] = v[k]
        k = k+1
    end
    return M
end

#End of Module
end