module Tscc


########### Imports ###########

using Oscar
using LinearAlgebra
using Polyhedra
using DataStructures
using CDDLib
using Permutations

include("./start_subdivision.jl")

using .start_subdivision




########### Redefined functions ###########

"""
    deepcopy(S::DeloneCell)

Create a deep copy of `S`.
"""
Base.:deepcopy(S::start_subdivision.DeloneCell) = start_subdivision.DeloneCell(deepcopy(S.Amount),deepcopy(S.Vertices))

"""
    +(S,p)

Return the DeloneCell `S` translated by the Vector `p`
"""
function Base.:+(S::start_subdivision.DeloneCell,p::Vector{Float64})
    T = Set([])
    for s in S.Vertices
        push!(T,s+p)
    end
    return start_subdivision.DeloneCell(S.Amount,T)
end

"""
    -(S,p)

Return the DeloneCell `S` translated by the Vector `-p`
"""
Base.:-(S::start_subdivision.DeloneCell,p::Vector{Float64}) = S + (-p)

"""
    isequal(S,T)

Check if the DeloneCells `S` and `T` are equal.
"""
function Base.:isequal(S::start_subdivision.DeloneCell,T::start_subdivision.DeloneCell)
    if S.Amount != T.Amount
        return false
    end
    if S.Vertices != T.Vertices
        return false
    end
    return true
end

"""
    ==(S,T)

Check if the Simplicies `S` and `T` are equal.
"""
Base.:(==)(S::start_subdivision.DeloneCell,T::start_subdivision.DeloneCell) = isequal(S,T)
"""
    hash(S)

Define the hash of `S` as the hash of the Vertex Set.
"""
Base.hash(S::start_subdivision.DeloneCell) = hash(S.Vertices)

"""
    hash(H)

Overload the hash function for a `HalfSpace`. Needed for the `haskey` function.
"""
Base.hash(H::HalfSpace{Rational{BigInt}, Vector{Rational{BigInt}}}) = hash(H.a)+hash(H.β)

"""
    hash(H)

Overload the hash function for a `HyperPlane`. Needed for the `haskey` function.
"""
Base.hash(Z::HyperPlane{Rational{BigInt}, Vector{Rational{BigInt}}}) = hash(Z.a)+hash(Z.β)




########### Utility Functions ###########

"""
    vec2mat(v)

return the symmetric Matrix `M` which is encoded in the vector `v`.
"""
function vec2mat(v)
    n = Int((sqrt(8*size(v)[1] +1)-1)/2)
    M = zeros(typeof(v[1]),(n,n))
    k = 1
    for i in 1:n
        for j in 1:i-1
            M[i,j] = v[k]
            M[j,i] = M[i,j]
            k = k+1
        end
        M[i,i] = v[k]
        k = k+1
    end
    return M
end

"""
    calcF(D)

Calculates the Matrix F for a given subdivision as in M.Sc Thesis Miebach Chapter 3.6.
"""
function calcF(D::start_subdivision.DeloneSub)
    F = zeros(D.Dimension,D.Dimension)

    d = D.Dimension
    H1 = collect(D.Cells)
    for (i, S1) in enumerate(H1) 
        P1 = collect(S1.Vertices)
        for k = 1:S1.Amount , j=1:k-1
            v = P1[k] - P1[j]
            reshape(v, length(v), 1)
            F += (1/S1.Amount)*(v*v')
        end
    end
    F = F*(1/(D.Number_Cells))
    return F
end

"""
    mat2string(M)

Converts a matrix M to a String s.
"""
function mat2string(M)
    n = size(M)[1]
    s = ""
    s *= string(n)*"x0\n"
    for i in 1:n
        for j in 1:i
            s *= string(Int(M[i,j]))*" "
        end
        s*="\n"
    end
    return s
end




########### Functions ###########

"""
    centralform(L::SecondaryCone)

Compute the T-centralform of the T-secondary cone defined by L and G.
"""
function centralform(L::start_subdivision.SecondaryCone,G::MatrixGroup,d::Int64)
    if length(L.Inequalities) !=0 && length(L.Equalities) !=0
        Q = collect(typeof(first(keys(L.Inequalities))),keys(L.Inequalities))
        H = collect(typeof(first(keys(L.Equalities))),keys(L.Equalities))
        P = Polyhedra.polyhedron(hrep(Q) ∩ hrep(H), CDDLib.Library(:exact))
    elseif length(L.Inequalities) !=0 && length(L.Equalities) ==0
        Q = collect(typeof(first(keys(L.Inequalities))),keys(L.Inequalities))
        P = Polyhedra.polyhedron(hrep(Q), CDDLib.Library(:exact))
    elseif length(L.Inequalities) ==0 && length(L.Equalities) !=0
        H = collect(typeof(first(keys(Delta.Equalities))),keys(L.Equalities))
        P = Polyhedra.polyhedron(hrep(H), CDDLib.Library(:exact))
    end

    Sub = invariant_symmetric_forms(G)
    M = zeros(length(Sub)+1,Int(d*(d+1)/2))
    for i =1:length(Sub)
        M[i,:] = start_subdivision.mat2vec(Matrix(Float64.(Sub[i])))
    end

    P1 = Polyhedra.polyhedron(Polyhedra.vrep(M), CDDLib.Library(:exact))

    h = collect(Polyhedra.hyperplanes(hrep(P1)))
    for i=1:length(h)
        P = Polyhedra.polyhedron(hrep(P) ∩ h[i], CDDLib.Library(:exact))
    end
    
    poly = vrep(P)

    centralform = zeros(Rational{BigInt},(d,d))
    for ray in Polyhedra.allrays(poly)
        R = vec2mat(ray.a)

        R *= lcm([a.den for a in R])
        R /= gcd(R)
        centralform += R
    end

    return centralform
end

"""
    GeneralVoronoi(G,d)

Starts the Enumeration of al T-equivalent t-generic triangulations (see Alg. 2, M.Sc Thesis Miebach)
"""
function GeneralVoronoi(G::MatrixGroup,d::Int64)
    println("Searching for start subdivision...")
    M,D = start_subdivision.TgenericDelone(G,d)
    println(" ")
    println("Found a Subdivision with ",D.Number_Cells," Delone polytopes in the origin star.")
    println("With PQF: ",M)

    #Creating the secondary cone
    L_redundant = start_subdivision.SecondaryCone(D)
    L = start_subdivision.reduce(L_redundant)

    Y = Queue{start_subdivision.SecondaryCone}()
    enqueue!(Y,L)

    R = Queue{start_subdivision.SecondaryCone}()
    cones = Set{start_subdivision.TSecondaryCone}()
    count = 1
    while !isempty(Y)
        println("Inequivalent T-secondary cones found so far: ", count)
        count = count + 1
        L_s = dequeue!(Y)

        M_T = centralform(L_s,G,d)
        M1 = Oscar.matrix(QQ, M_T)
        println(" ")
        println("T-central form: ", M_T)

        enqueue!(R,L_s)

        println(" ")
        #Creating the T-secondary cone with the flips for each facet
        D_ex = extendedSubdivision(M1,d)
        D_ex_int = start_subdivision.toIntegerLattice(D_ex,M1)
        L_ex_redundant = start_subdivision.SecondaryCone(D_ex_int)
        L_ex = start_subdivision.reduce(L_ex_redundant)
        P,P1 = start_subdivision.createTSecondaryCone(L_ex,G,d)
        L_T_ex = start_subdivision.TSecondaryCone(P,L_ex,L_s.DeloneSub)

        push!(cones, L_T_ex)
        (Polyhedra.dim(P) == 1) && break

        Dk = Queue{start_subdivision.SecondaryCone}()
        count1 =1
        for (H,PolytopePairs) in L_T_ex.Inequalities
            println(" ")
            println("Neighbor number: ",count1)
            count1 += 1
            PolytopeChains = RepartitioningPolytopes(PolytopePairs)
            Di = deepcopy(L_s.DeloneSub)
            D1,D2 = flip(L_s.DeloneSub,PolytopeChains,M_T,d)
            D1.Cells == L_s.DeloneSub.Cells ? Di = D2 : (D2.Cells == L_s.DeloneSub.Cells && (Di = D1) )

            L_redundant1 = start_subdivision.SecondaryCone(Di)
            L1 = start_subdivision.reduce(L_redundant1)
            M_new = centralform(L1,G,d)
            println("Found form: ",M_new)
    
            P,P1=start_subdivision.createTSecondaryCone(L1,G,d)

            if start_subdivision.TrigidityIndex(P,P1) == false
                println(" ")
                println("T-Secondary Cone not T-generic")
                continue
            end

            if start_subdivision.ValidFlipCheck(L1,L_s) == false
                println(" ")
                println("Secondary cone dimensions do not coincide") 
                continue
            end

            lower_bound_pre_check = lower_bound(Di,R ∪ Dk,d)
            break_check = false
            for L_check ∈ lower_bound_pre_check
                equivalent(L1,L_check,G,d) && ((break_check = true) && break)
            end
            if break_check == true
                println(" ")
                println("Equivalent")
                continue
            end

            println(" ")
            println("Inequivalent")
            enqueue!(Y,L1)
            enqueue!(Dk,L1)
        end
        empty!(Dk)
    end
    return cones
end

"""
    equivalent(L_D1,L_D2,G,d)

Check for GL_n(Z) equivalence between two secondary coney by their central forms.
"""
function equivalent(L_D1,L_D2,G,d)

    M = centralform(L_D1,G,d)
    N = centralform(L_D2,G,d)

    PATH_ISOM = "ISOM_and_AUTO/"
    ISOM_C = PATH_ISOM * "ISOM.c"
    ISOM_EXE = PATH_ISOM * "ISOM"
    FILE = PATH_ISOM * "tmp"

    # Convert matrices to string
    stringM = mat2string(M)
    stringN = mat2string(N)

    # Prepare input for the C program
    isom_input = stringM * '\n' * stringN

    # Write input to a temporary file
    write(FILE, isom_input)

    # Run the compiled executable and read the output
    out = read(`$ISOM_EXE $FILE`, String)

    # Remove the temporary file
    rm(FILE)
    
    if startswith(out,"The lattices are not isomorphic")
        return false
    else
        return true
    end
end


"""
    RepartitioningPolytopes(H,D)

Creates the equvalence classes described in [1] Section 3.3 by construction all "Polytope Chains" for a facet.

[1] A generalization of Voronoi's reduction theory and its application, Dutour Sikirić, Mathieu and Schürmann, Achill and Vallentin, Frank.
"""
function RepartitioningPolytopes(P::Set{Tuple{start_subdivision.DeloneCell,start_subdivision.DeloneCell}})
    Pairs = copy(P)
    Pairs_iter = collect(Pairs)
    RepartPolytopes = Set{Set{start_subdivision.DeloneCell}}()
    while !isempty(Pairs_iter)
        Pair = pop!(Pairs_iter)
        for (i,Pair1) in enumerate(Pairs_iter)
            if Pair1[1] == reverse(Pair)[1] && Pair1[2] == reverse(Pair)[2]
                popat!(Pairs_iter,i)
            end
        
        end
        A = Pair[1]
        B = Pair[2]

        RepartPolytope = Set{start_subdivision.DeloneCell}()
        push!(RepartPolytope,A)
        push!(RepartPolytope,B)
        while true
            breakcheck = true
            for (i,Pair1) in enumerate(Pairs_iter)
                if Pair1[1] == A
                    push!(RepartPolytope,Pair1[2])
                    A = deepcopy(Pair1[2])
                    popat!(Pairs_iter,i)
                    for (j,Pair_pop) in enumerate(Pairs_iter)
                        if Pair_pop[1] == reverse(Pair1)[1] && Pair_pop[2] == reverse(Pair1)[2]
                            popat!(Pairs_iter,j)
                        end
                    end
                    breakcheck = false
                    break
                elseif Pair1[2] == B
                    push!(RepartPolytope,Pair1[1])
                    B = deepcopy(Pair1[1])
                    popat!(Pairs_iter,i)
                    for (j,Pair_pop) in enumerate(Pairs_iter)
                        if Pair_pop[1] == reverse(Pair1)[1] && Pair_pop[2] == reverse(Pair1)[2]
                            popat!(Pairs_iter,j)
                        end
                    end
                    breakcheck = false
                    break
                elseif Pair1[2] == A
                    push!(RepartPolytope,Pair1[1])
                    A = deepcopy(Pair1[1])
                    popat!(Pairs_iter,i)
                    for (j,Pair_pop) in enumerate(Pairs_iter)
                        if Pair_pop[1] == reverse(Pair1)[1] && Pair_pop[2] == reverse(Pair1)[2]
                            popat!(Pairs_iter,j)
                        end
                    end
                    breakcheck = false
                    break
                elseif Pair1[1] == B
                    push!(RepartPolytope,Pair1[2])
                    B = deepcopy(Pair1[2])
                    popat!(Pairs_iter,i)
                    for (j,Pair_pop) in enumerate(Pairs_iter)
                        if Pair_pop[1] == reverse(Pair1)[1] && Pair_pop[2] == reverse(Pair1)[2]
                            popat!(Pairs_iter,j)
                        end
                    end
                    breakcheck = false
                    break
                end
            end
            if breakcheck == true
                push!(RepartPolytopes,RepartPolytope)
                break
            end
        end
    end
    return RepartPolytopes
end

"""
    extendedSubdivision(Q,d,n) 
      
Calculates the extended subdivision for a given positive definite quadratic form. [1](Main algorithm)

[1] Complexity and algorithms for computing Voronoi cells of lattices, Dutour Sikirić, Mathieu and Schürmann, Achill and Vallentin, Frank
"""
function extendedSubdivision(Q::QQMatrix,d::Int64)
    V = Set([])
    iter = 0
    while true
        iter += 1
        x,V = start_subdivision.startVertex(Q,d)
        length(V) >= d+1 && break
        if iter >= 1000
            return 0
        end
    end
    T = Queue{start_subdivision.DeloneCell}()
    enqueue!(T,start_subdivision.DeloneCell(length(V),V))
    M = Queue{start_subdivision.DeloneCell}()
    facet_size = 0

    for j =1:1000 #Has to be manualy adjusted because it is not clear form the begining how for one has to claculate the extended subdivision e.g 1500 doesnt work
        D = dequeue!(T)
        enqueue!(M,D)
        ℱ = start_subdivision.getFacets(D,d)
        ℱ1 = collect(ℱ.Facets)
        for (i,F) in enumerate(ℱ1)
            facet_size = length(F.Vertices)
            D_bar = start_subdivision.getNeighbor(F,Q,d)
            if !any([D_bar == D1.Vertices for D1 ∈ T ∪ M])
                enqueue!(T,start_subdivision.DeloneCell(length(D_bar),D_bar))
            end
        end
    end
    return start_subdivision.DeloneSub(T ∪ M,d,facet_size)
end

"""
    flip(D,PolytopeChains,M,d)

Calcutales the flip of a PolytopeChain (repartitioning polytope) by calculating the lifting map and changing upper to lower facets.
"""
function flip(D::start_subdivision.DeloneSub,PolytopeChains::Set{Set{start_subdivision.DeloneCell}},M::Matrix,d::Int64)
    D1 = deepcopy(D)
    D2 = deepcopy(D)
    upper_facets = Set{Set{}}()
    lower_facets =Set{Set{}}()
    bool = true
    P_Chains = collect(PolytopeChains)
    for (i,Chain) in enumerate(P_Chains)
        RepartitioningPolytope = Set{}()
        Chain1 = collect(Chain)
        for (j,Cell) in enumerate(Chain1)
            for (k,Vertex) in enumerate(Cell.Vertices)
                push!(RepartitioningPolytope,Vertex)
            end
        end
        Pol = collect(RepartitioningPolytope)
        Pol_Mat = zeros(length(Pol),d)
        for j =1:length(Pol)
            for k=1:d
                Pol_Mat[j,k] = Pol[j][k]
            end
        end
        P = Polyhedra.polyhedron(Polyhedra.vrep(Pol_Mat))
        V_c = LatticePoints(P,Pol,d)
        Lift_Mat = zeros(length(V_c),d+1)
        for j = 1:length(V_c)
            for k =1:d
                Lift_Mat[j,k] = V_c[j][k]
            end
            Lift_Mat[j,d+1] = transpose(V_c[j])*M*V_c[j]        
        end
        P_lift = Polyhedra.vrep(Lift_Mat)
        P_lift_verts = collect(Polyhedra.points(P_lift))
        H_R = collect(Polyhedra.halfspaces(Polyhedra.hrep(Polyhedra.polyhedron(P_lift))))
        
        for (j,H) in enumerate(H_R)
            if H.a[d+1] < -(10^-5)
                lower_cell = Set{}()
                for (k,V) in enumerate(P_lift_verts)
                    if abs.(dot(H.a,V) - H.β) <= 10^-5
                        push!(lower_cell,V[1:d])
                    end
                end
                push!(lower_facets,lower_cell)
            end
            if H.a[d+1] > (10^-5)
                upper_cell = Set{}()
                for (k,V) in enumerate(P_lift_verts)
                    if abs.(dot(H.a,V) - H.β) <= 10^-5
                        push!(upper_cell,V[1:d])
                    end
                end
                push!(upper_facets,upper_cell)
            end
        end
    end
    Upper_facets1 = collect(upper_facets)
    Lower_facets1 = collect(lower_facets)
    for (i,facet) in enumerate(Upper_facets1)
        pop!(D1.Cells,start_subdivision.DeloneCell(length(facet),facet),0)
        pop!(D2.Cells,start_subdivision.DeloneCell(length(facet),facet),0)
    end
    for (i,facet) in enumerate(Lower_facets1)
        pop!(D1.Cells,start_subdivision.DeloneCell(length(facet),facet),0)
        pop!(D2.Cells,start_subdivision.DeloneCell(length(facet),facet),0)
    end
    for (i,facet) in enumerate(Upper_facets1)
        if zeros(d) ∈ facet
            push!(D1.Cells,start_subdivision.DeloneCell(length(facet),facet))
        end
    end
    for (i,facet) in enumerate(Lower_facets1)
        if zeros(d) ∈ facet
            push!(D2.Cells,start_subdivision.DeloneCell(length(facet),facet))
        end
    end
    return D1,D2
end

"""
    LatticePoints(P,vertices,d)

    Finds all lattice points of a polyhedron "P" and its vertices "verticies".
    For the enumeration of the cartesian coordiantes in rang of "sz" see [https://juliapackages.com/p/cartesian] 
    bottom of the page.
"""
function LatticePoints(P::Polyhedra.Polyhedron,vertices::Vector{},d::Int64)
    lb = copy(vertices[1][1])
    ub = copy(vertices[1][1])
    for (i,vertex) in enumerate(vertices)
        for j = 1:length(vertex)
            vertex[j] < lb && (lb = vertex[j])
            vertex[j] > ub && (ub = vertex[j]) 
        end
    end
    sz = (ub-lb+1)*ones(d)
    
    LP = Set{}()
    if !(isempty(sz) || prod(sz) == 0)
        N = length(sz)
        c = ones(Int, N)
        sz1 = sz[1]
        isdone = false
        while !isdone
            c1 = (c.-((1-lb)*ones(d)))
            in(c1,Polyhedra.hrep(P)) && push!(LP,c1)         # This is whatever code we put inside the "loop"
            if (c[1]+=1) > sz1
                idim = 1
                while c[idim] > sz[idim] && idim < N
                    c[idim] = 1
                    idim += 1
                    c[idim] += 1
                end
                isdone = c[end] > sz[end]
            end
        end
    end
    return collect(LP)
end

"""
    lower_bound(D1,Q,d)

    Checks if the lower bound of a triangulation D1 coincides with one that has already been found in Q.
"""
function lower_bound(D1,Q,d)
    R = Queue{start_subdivision.SecondaryCone}()

    max_amount_d1 = 0
    D1_vec = collect(D1.Cells)
    for (i,Cell) in enumerate(D1_vec)
        Cell.Amount > max_amount_d1 && (max_amount_d1 = Cell.Amount) 
    end
    theta1 = sqrt(((d/(max_amount_d1))^d)*det(calcF(D1)))

    for L2 ∈ Q
        max_amount_d2 = 0
        D2_vec = collect(L2.DeloneSub.Cells)
        for (i,Cell) in enumerate(D2_vec)
            Cell.Amount > max_amount_d2 && (max_amount_d2 = Cell.Amount) 
        end

        theta2 = sqrt(((d/(max_amount_d2))^d)*det(calcF(L2.DeloneSub)))
        #println("New form: ",theta1)
       #println("compared to: ",theta2)
        if abs(theta1-theta2) <= 10^-7
            enqueue!(R,L2)
        end
    end
    return R
end

end