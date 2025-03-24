module scc
########### Imports ###########

using Permutations
using Polyhedra
using CDDLib
using LinearAlgebra
using DataStructures
using Oscar
using Base.Threads




########### Structs ###########

"""
    Simplex

Store a Simplex as a set of Vertices which are Vectors of `Rational{BigInt}`.
"""
struct Simplex
    Dimension::Int
    Vertices::Set{Vector{Rational{BigInt}}}
end

"""
    Delone

Store a Delone Triangulation as a set of `Simplex`'s.
"""
struct Delone
    Dimension::Int
    Number_Simplices::Int
    Simplices::Set{Simplex}
end

"""
    SecondaryCone

Struct containing the `Delone` Triangulation and a Dict contaning `HalfSpaces` (i.e. Regulators) and a set of all coressponding Pairs of Simplices.
"""
struct SecondaryCone
    Dimension::BigInt
    Delone::Delone
    Inequalities::Dict{HalfSpace,Set{Tuple{Simplex,Simplex}}}
end




########### Redefined functions ###########

"""
    deepcopy(S::Simplex)

Create a deep copy of `S`.
"""
Base.:deepcopy(S::Simplex) = Simplex(deepcopy(S.Dimension),deepcopy(S.Vertices))

"""
    +(S,p)

Return the Simplex `S` translated by the Vector `p`
"""
function Base.:+(S::Simplex,p::Vector{Rational{BigInt}})
    T = Set([])
    for s in S.Vertices
        push!(T,s+p)
    end
    return Simplex(S.Dimension,T)
end

"""
    -(S,p)

Return the Simplex `S` translated by the Vector `-p`
"""
Base.:-(S::Simplex,p::Vector{Rational{BigInt}}) = S + (-p)

"""
    isequal(S,T)

Check if the Simplices `S` and `T` are equal.
"""
function Base.:isequal(S::Simplex,T::Simplex)
    if S.Dimension != T.Dimension
        return false
    end
    if S.Vertices != T.Vertices
        return false
    end
    return true
end

"""
    ==(S,T)

Check if the Simplices `S` and `T` are equal.
"""
Base.:(==)(S::Simplex,T::Simplex) = isequal(S,T)
"""
    hash(S)

Define the hash of `S` as the hash of the Vertex Set.
"""
Base.hash(S::Simplex) = hash(S.Vertices)




########### Functions ###########

"""
    regulator(S1,S2)

Return the Regulator of the adjacent Simplices `S1` and `S2`.

First a affine relationship `alpha` is computed and then the regulator is computed via Baranovskiĭ formula (cf [1]).

[1] Baranovskiĭ, E. P. , Simplexes of L-subdivisions of Euclidean spaces. Math. Notes 10, 827--834 (1972; Zbl 0234.50008)
"""
function regulator(S1::Simplex,S2::Simplex)
    # Create Arrays for a fixed ordering
    dimension = S1.Dimension
    P1 = collect(S1.Vertices)
    P2 = collect(S2.Vertices)

    # Calculate an affine relationship alpha
    B = Array{Rational{BigInt}}(undef, dimension+1, dimension+1)
    for (i,vector) in enumerate(P1)
        for j in 1:dimension
            B[j,i] = vector[j]
        end
    end
    for i in 1:dimension+1
        B[dimension+1,i] = 1
    end
    v = append!(deepcopy(first(setdiff(P2,P1))),1)
    alpha = B\v

    # Compute regulator
    N = fill(Rational{BigInt}(0), (dimension, dimension))
    for i in 1:dimension+1, j in 1:i-1
        v = P1[i]-P1[j]
        reshape(v, length(v), 1)
        N += alpha[i]*alpha[j]*v*v'
    end
    return N
end

"""
    Delone(n)

Return the Delone Triangulation of Voronoĭ's principal form of the first type in dimension `n`.
"""
function Delone(n::Int)
    E = hcat(I(n),-ones(Rational{BigInt},n))
    D = Delone(n,factorial(n+1),Set{Simplex}())
    for pi in PermGen(n+1)
        L = Simplex(n,Set{Vector{Rational{BigInt}}}())
        for j ∈ 1:n+1
            push!(L.Vertices,sum([E[:,pi(i)] for i ∈ 1:j]))
        end
        push!(D.Simplices,L)
    end
    return D
end

"""
    adjacent(S::Simplex,D::Delone)

Return a Set of all Simplices `T` that are adjacent over a common facet with `S`.
Used for the repartitioning.
"""
function adjacent(S::Simplex, D::Delone)
    Inzident = Set{Simplex}([])
    for vertex in S.Vertices
        for Simplex in D.Simplices
            TranslatedSimplex = Simplex + vertex
            if length(intersect(S.Vertices,TranslatedSimplex.Vertices)) == S.Dimension
                push!(Inzident,deepcopy(TranslatedSimplex))
            end
        end
    end
    return Inzident
end

"""
    adjacent(S::Simplex,D::Delone)

Return a Set of all Simplices `T` in the origin star that are adjacent over a common facet with `S`.
Used for creating the secondary cone. (Faster because one less loop)
"""
function Staradjacent(S::Simplex, D::Delone)
    Inzident = Set{Simplex}([])
    for Simplex in D.Simplices
        if length(intersect(S.Vertices,Simplex.Vertices)) == S.Dimension
            push!(Inzident,deepcopy(Simplex))
        end
    end
    return Inzident
end

"""
    adjacent(D::Delone)

Return a set of tuples `(S,adj)` where `S` is a Simplex in `D` and `adj` contains all Simplices `T` that are adjacent over a common facet with `S`.
"""
function adjacent(D::Delone)
    n = nthreads()  # Get number of threads
    adj_chunks = [Set{Tuple{Simplex, Set{Simplex}}}() for _ in 1:n]  # Each thread gets its own Set

    cole = collect(D.Simplices)  # Convert Set to a Vector

    @threads for i in eachindex(cole)
        S = cole[i]
        Set_adjacent_S = adjacent(S, D)
        push!(adj_chunks[threadid()], (S, Set_adjacent_S))  # No lock needed
    end

    return union(adj_chunks...)  # Merge all results at the end
end

"""
    Staradjacent(D::Delone)

Return a set of tuples `(S,adj)` where `S` is a Simplex in `D` and `adj` contains all Simplices `T` in the origin star that are adjacent over a common facet with `S`.
"""
function Staradjacent(D::Delone)
    n = nthreads()  # Get number of threads
    adj_chunks = [Set{Tuple{Simplex, Set{Simplex}}}() for _ in 1:n]  # Each thread gets its own Set

    cole = collect(D.Simplices)  # Convert Set to a Vector

    @threads for i in eachindex(cole)
        S = cole[i]
        Set_adjacent_S = Staradjacent(S, D)
        push!(adj_chunks[threadid()], (S, Set_adjacent_S))  # No lock needed
    end

    return union(adj_chunks...)  # Merge all results at the end
end

"""
    flip(D,SimplexPairs)

Perform a bistellar flip of the Delone Triangulation `D` with a given list of all pairs of Simplices whose regulator defines a facet of the SecondaryCone.
"""
function flip(D,SimplexPairs)
    dimension = D.Dimension

    FlippablePairs = unique([Set([S,T]) for (S,T) in SimplexPairs])
    D1 = deepcopy(D)
    D2 = deepcopy(D)

    for (S,T) in FlippablePairs

        P = collect(S.Vertices)
        Q = collect(T.Vertices)

        B = Array{Rational{BigInt}}(undef, dimension+1, dimension+1)
        for (i,vector) in enumerate(P)
            for j in 1:dimension
                B[j,i] = vector[j]
            end
        end
        for i in 1:dimension+1
            B[dimension+1,i] = 1
        end
        x = deepcopy(first(setdiff(Q,P)))
        v = append!(x,1)
        alpha = B\v
        pop!(x)
        append!(P,[x])
        append!(alpha,-1)

        # Create Indexsets and corresponding Simplices
        Ipos = [p for (index,p) in enumerate(P) if alpha[index] < 0]
        Ineg = [p for (index,p) in enumerate(P) if alpha[index] > 0]
        Spos = [Simplex(dimension,Set(setdiff(P,[v]))) for v ∈ Ipos]
        Sneg = [Simplex(dimension,Set(setdiff(P,[v]))) for v ∈ Ineg]

        for S in Spos
            pop!(D1.Simplices,S,0)
            push!(D2.Simplices,S)
        end

        for T in Sneg
            push!(D1.Simplices,T)
            pop!(D2.Simplices,T,0)
        end
    end

    # Remove Simplices that are not incident with 0.
    zero = zeros(Rational{BigInt},dimension)
    for S in D1.Simplices
        zero ∉ S.Vertices && pop!(D1.Simplices,S,0)
    end
    for S in D2.Simplices
        zero ∉ S.Vertices && pop!(D2.Simplices,S,0)
    end

    return(D1,D2)
end


"""
    equivalent(D1,D2)

Check if the Delone Triangulation `D1,D2` are GL_n(Z) equivialent.

This is done by checking if there exists a GL_n(Z) matrix between the two triangulations.
"""
#Equivalence check with permutations
#=
function equivalent(D1,D2)
    #Imput: two triangulation which need to be checked for algebraic equality.
    #Reutrn: true if equal to a triangulation in D false if not

    d = D1.Dimension

    if D1.Number_Simplices != D2.Number_Simplices
        return false
    end

    
    pt =  Vector{Matrix{Rational{BigInt}}}(undef,length(D1.Simplices))
    D =  Vector{Matrix{Rational{BigInt}}}(undef,length(D2.Simplices))


    H1 = collect(D1.Simplices)
    for (i, S1) in enumerate(H1)
        P1 = collect(S1.Vertices)
        MatTemp = zeros(D1.Dimension,D1.Dimension+1)
        for (j,V) in enumerate(P1)
            MatTemp[:,j] = V
        end
        pt[i] = MatTemp
    end

    H2 = collect(D2.Simplices)
    for (i, S2) in enumerate(H2)
        P2 = collect(S2.Vertices)
        MatTemp = zeros(D2.Dimension,D2.Dimension+1)
        for (j,V) in enumerate(P2)
            MatTemp[:,j] = V
        end
        D[i] = MatTemp
    end
 
    #Remove the zeros out of the matrices which describe the triangles to get a basis represantation of the potential triangulation
    base_pt = Matrix{Float64}[]
    for i = 1:D1.Number_Simplices
        index = 0
        for j = 1:d+1
            if pt[i][:,j] == zeros(d)
                index = j
                break
            end
        end
        push!(base_pt,pt[i][:,1:end .!=index])
    end

    #Remove the zeros out of the matrices which describe the triangles to get a basis represantation of the current triangulation D[k]
    base_D = Matrix{Float64}[]
    for i = 1:D1.Number_Simplices 
        index = 0
        for j = 1:d+1
            if D[i][:,j] == zeros(d)
                index = j
                break
            end
        end
        push!(base_D,D[i][:,1:end .!=index])
    end


    #Check if there exists a base change matrix B in GL_d(Z) s.t. B * base_triangle = base_pt[i] for all i
    # Where B is calculated by permuting the base triangles and inverting them 
    base_triangle = base_D[1]
    potential_base_change_matrix = Matrix{Float64}
    P = PermGen(d)
    for i = 1:D1.Number_Simplices
        test_pt = copy(base_pt)
        for y in P
            y_1 = two_row(y)
            base_triangle_perm = zeros(d,d)
            for l = 1:d
                base_triangle_perm[:,l] = base_triangle[:,y_1[2,l]]
            end


            potential_base_change_matrix = base_pt[i] * inv(base_triangle_perm) 


            if abs(Int(round(det(potential_base_change_matrix)))) != 1
                continue
            end
            if potential_base_change_matrix.%1 != zeros(d,d)
                continue
            end

            #Check if for a base change matrix the triangulations are equal 
            for q in P
                q_1 = two_row(q)
                for j = 1:D1.Number_Simplices
                    H = potential_base_change_matrix * base_D[j]
                    H_perm = zeros(d,d)
                    for l = 1:d
                        H_perm[:,l] = H[:,q_1[2,l]]
                    end
                    test_pt = filter!(!=(H_perm),test_pt)
                    if length(test_pt) == 0
                        return true
                    end
                end
            end
        end
        #if length(test_pt) == 0
        #    return true
        #end
    end
    return false
end=#

"""
    equivalent(L_D1,L_D2)

Check if the Secondary cones `L_D1,L_D2` are GL_n(Z) equivialent.

This is done by computing the centralform for each and then checking those for equivalence with the ISOM programm. cf. [1] and [2].

[1] M. Dutour Sikirić et al., Acta Crystallogr., Sect. A 72, No. 6, 673--683 (2016; Zbl 1370.82063)
[2] http://www.math.uni-rostock.de/~waldmann/ISOM_and_AUTO.zip
"""
function equivalent(L_D1,L_D2)

    M = centralform(L_D1)
    N = centralform(L_D2)

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

"""
    mat2vec(M)

Return a a vector of size n*(n+1)/2 containing the rows of the lower triangular part of `M`.

See also [`vec2mat`](@ref)
"""
function mat2vec(M)
    n = size(M)[1]
    v = zeros(typeof(M[1,1]),Int((n*(n+1)/2)))
    k = 1
    for i in 1:n
        for j in 1:i-1
            v[k] = M[j,i]+M[i,j]
            k = k+1
        end
        v[k] = M[i,i]
        k = k+1
    end
    return v
end

"""
    SecondaryCone(D::Delone)

Create the `SecondaryCone` for a given `Delone` Triangulation.(Used for repartitioning)
"""
function SecondaryCone(D::Delone)
    Inequalities = Dict{HalfSpace{Rational{BigInt}, Vector{Rational{BigInt}}},Set{Tuple{Simplex,Simplex}}}()
    for (S,listS) in adjacent(D), T in listS
        H = HalfSpace(mat2vec(regulator(S,T)),0)
        # If the Inequality exits append the Tuple `(S,T)` if not create a new Set.
        haskey(Inequalities, H) ? push!(Inequalities[H],(S,T)) : Inequalities[H] = Set([(S,T)])
    end
    return SecondaryCone(D.Dimension,D,Inequalities)
end

"""
    StarSecondaryCone(D::Delone)

Create the `SecondaryCone` for a given `Delone` Triangulation.(Used for checking equivalence)
"""
function StarSecondaryCone(D::Delone)
    Inequalities = Dict{HalfSpace{Rational{BigInt}, Vector{Rational{BigInt}}},Set{Tuple{Simplex,Simplex}}}()
    for (S,listS) in Staradjacent(D), T in listS
        H = HalfSpace(mat2vec(regulator(S,T)),0)
        # If the Inequality exits append the Tuple `(S,T)` if not create a new Set.
        haskey(Inequalities, H) ? push!(Inequalities[H],(S,T)) : Inequalities[H] = Set([(S,T)])
    end
    return SecondaryCone(D.Dimension,D,Inequalities)
end

""" 
    reduce(L::SecondaryCone)

Create a new `SecondaryCone` w/o redundent inequalities.
"""
function reduce(L::SecondaryCone)
    Q = collect(typeof(first(keys(L.Inequalities))),keys(L.Inequalities))
    P = Polyhedra.polyhedron(hrep(Q), CDDLib.Library(:exact))
    removehredundancy!(P)
    
    Facets = Dict{HalfSpace, Set{Tuple{Simplex, Simplex}}}()
    for (Halfspace,SetSimplices) in L.Inequalities
        Halfspace ∈ halfspaces(P) ? Facets[Halfspace] = SetSimplices : continue
    end
    return SecondaryCone(L.Dimension, L.Delone,Facets)
    
end

"""
    Voronoi(n)

    Enumeration of all non-equivalent Delone triangulations in dimension `d`.
"""
function Voronoi(d)
    D = Delone(d)

    Y = Queue{SecondaryCone}()
    enqueue!(Y,reduce(SecondaryCone(D)))

    R = Queue{SecondaryCone}()
    cones = Set{SecondaryCone}()
    count = 1
    while !isempty(Y)
        println("Inequivalent cones found so far: ", count)
        count = count +1
        L = dequeue!(Y)
        L = reduce(SecondaryCone(L.Delone))
        enqueue!(R,L)
        push!(cones, L)
        
        Dk = Queue{SecondaryCone}()
        for (H,SimplexPairs) in L.Inequalities
            Di = deepcopy(L.Delone)
            D1,D2 = flip(L.Delone,SimplexPairs)
            D1.Simplices == L.Delone.Simplices ? Di = D2 : (D2.Simplices == L.Delone.Simplices && (Di = D1) )
            triang_pre_check = lower_bound(Di,R ∪ Dk,d)
            L_Di = reduce(StarSecondaryCone(Di))
            break_check = false
            for L_check ∈ triang_pre_check
                (equivalent(L_Di,L_check) && ((break_check = true) && break))
            end
            (break_check == true) && continue
            enqueue!(Y,L_Di)
            enqueue!(Dk,L_Di)
        end
        empty!(Dk)
    end
    return cones
end

"""
    lower_bound(D1,Q,d)

    Checks if the lower bound of a triangulation D1 coincides with one that has already been found in Q.
"""
function lower_bound(D1,Q,d)
    R = Queue{SecondaryCone}()
    theta1 = sqrt(((d/(d+1))^d)*det(calcF(D1)))
    for D2 ∈ Q
        theta2 = sqrt(((d/(d+1))^d)*det(calcF(D2.Delone)))
        if abs(theta1-theta2) <= 10^-7
            enqueue!(R,D2)
        end
    end
    return R
end

"""
    calcF(D)

Calculates the Matrix F for a given triangulation as in PHD Vallentin 8.3.1
"""
function calcF(D::Delone)
    F = zeros(D.Dimension,D.Dimension)

    d = D.Dimension
    H1 = collect(D.Simplices)
    for (i, S1) in enumerate(H1) 
        P1 = collect(S1.Vertices)
        for k = 1:d+1 , j=1:k-1
            v = P1[k] - P1[j]
            reshape(v, length(v), 1)
            F += v*v'
        end
    end
    F = F*(1/(D.Number_Simplices * (d+1)))
    return F
end

"""
    centralform(L::SecondaryCone)

Compute the centralform of the `SecondaryCone` `L`. cf [1]

[1] M. Dutour Sikirić et al., Acta Crystallogr., Sect. A 72, No. 6, 673--683 (2016; Zbl 1370.82063)
"""
function centralform(L::SecondaryCone)
    n = L.Dimension
    Ine = collect(L.Inequalities)
    poly = vrep(Polyhedra.polyhedron(hrep([H for (H,SP) in Ine])))
    centralform = zeros(Rational{BigInt},(n,n))
    time = @elapsed for ray in Polyhedra.rays(poly)
        R = vec2mat(ray.a)
        R *= lcm([a.den for a in R])
        R /= gcd(R)
        centralform += R
    end
    return centralform
end


"""
    mat2string(M)

Return a string encoding the lower triangular part of M for use with ISOM. cf [1].

[1] http://www.math.uni-rostock.de/~waldmann/ISOM_and_AUTO.zip
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
#End of Module
end