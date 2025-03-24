module start_subdivision


########### Imports ###########

using Oscar
using LinearAlgebra
using Polyhedra
using DataStructures
using CDDLib
using StatsBase
using Base.Threads




########### Structs ###########

"""
    DeloneCell

Store a DeloneCell as a set of Vertices which are Vectors of `Rational{BigInt}`.
"""
struct DeloneCell
    Amount::Int
    Vertices::Set{}
end

"""
    DeloneSub

Store a Delone subdivision as a set of 'DeloneCells'.
"""
struct DeloneSub
    Dimension::Int
    Number_Cells::Int
    Amount_facet_vertices::Int
    Cells::Set{DeloneCell}
end

"""
    Facet

Stores a Facet of a Delone Cell as a Set of vectors and as a Halfspace.
"""
struct Facet
    a::Vector{Rational{BigInt}}
    b::Rational{BigInt}
    Vertices::Set{}
end

"""
    SecondaryCone

Struct containing the `Delone` Subdivision and a Dict contaning `HalfSpaces` (i.e. Regulators) and a set of all coressponding Pairs of Delone Cells.
"""
struct SecondaryCone
    Dimension::BigInt
    DeloneSub::DeloneSub
    Inequalities::Dict{HalfSpace,Set{Tuple{DeloneCell,DeloneCell}}}

    #Its not really a facet but we have the struct already there so we use it anyway
    Equalities::Dict{HyperPlane,Set{DeloneCell}}
end

"""
    TSecondaryCone

Struct containing the `Delone` Subdivision and a Dict contaning `HalfSpaces` (i.e. Regulators) and a set of all coressponding Pairs of Delone Cells.
"""
struct TSecondaryCone
    Dimension::BigInt
    DeloneSub::DeloneSub
    Inequalities::Dict{HalfSpace,Set{Tuple{DeloneCell,DeloneCell}}}
    Equalities::Set{HyperPlane}
end


"""
    Facets

Stores a Facets of a Delone Cell as a Set of Facets.
"""
struct Facets
    Number_Facets::Int
    Facets::Set{Facet}
end




########### Redefined functions ###########

"""
    deepcopy(S::DeloneCell)

Create a deep copy of `S`.
"""
Base.:deepcopy(S::DeloneCell) = DeloneCell(deepcopy(S.Amount),deepcopy(S.Vertices))

"""
    +(S,p)

Return the DeloneCell `S` translated by the Vector `p`
"""
function Base.:+(S::DeloneCell,p::Vector{Float64})
    T = Set([])
    for s in S.Vertices
        push!(T,s+p)
    end
    return DeloneCell(S.Amount,T)
end

"""
    -(S,p)

Return the DeloneCell `S` translated by the Vector `-p`
"""
Base.:-(S::DeloneCell,p::Vector{Float64}) = S + (-p)

"""
    isequal(S,T)

Check if the DeloneCells `S` and `T` are equal.
"""
function Base.:isequal(S::DeloneCell,T::DeloneCell)
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
Base.:(==)(S::DeloneCell,T::DeloneCell) = isequal(S,T)
"""
    hash(S)

Define the hash of `S` as the hash of the Vertex Set.
"""
Base.hash(S::DeloneCell) = hash(S.Vertices)

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
    toQLattice(D,Q)

    Transforms a Delone subdivision form the lattice Z to the lattice Q
"""
function toQLattice(D::DeloneSub,Q::QQMatrix)
    D_Q = DeloneSub(D.Dimension,D.Number_Cells,D.Amount_facet_vertices,Set{DeloneCell}())
    C = cholesky(Matrix(Float64.(Q)))
    D1 = collect(D.Cells)
    for (i,Cell) in enumerate(D1)
        Cell_Q = Set{}()
        C1 = collect(Cell.Vertices)
        for (j,v) in enumerate(C1)
            push!(Cell_Q,Vector(Float64.(C.U*v)))
        end
        push!(D_Q.Cells, DeloneCell(length(C1),Cell_Q))
    end
    return D_Q
end

"""
    toIntegerLattice(D,Q)

    Transforms a Delone subdivision from the lattice Q to the lattice Z
"""
function toIntegerLattice(D::DeloneSub,Q::QQMatrix)
    D_integer = DeloneSub(D.Dimension,D.Number_Cells,D.Amount_facet_vertices,Set{DeloneCell}())
    C = cholesky(Matrix(Float64.(Q)))
    D1 = collect(D.Cells)
    for (i,Cell) in enumerate(D1)
        Cell_integer = Set{}()
        C1 = collect(Cell.Vertices)
        for (j,v) in enumerate(C1)
            a = Vector(round.(Float64.(inv(C.U)*v')))
            for index=1:length(a)
                if a[index] == -0.0
                    a[index] = 0.0
                end
            end
            push!(Cell_integer,a)
        end
        push!(D_integer.Cells, DeloneCell(length(C1),Cell_integer))
    end
    return D_integer
end

"""
    closest_vectors

x::A vector in the lattice defined by Q. Returns the set of closest vectors regarding Q. Yet Oscar.close_vectors requiers the integer lattice.
"""
function closest_vectors(Q::QQMatrix,x::Vector)
    C = cholesky(Matrix(Float64.(Q)))
    L = integer_lattice(gram = Q)
    Cv = []
    temp = inv(C.U)*x
    temp[abs.(temp) .< 10^-5] .= 0

    j=1
    while isempty(Cv)
        Cv = close_vectors(L,Vector(Rational.(temp)),j)
        j=j+1
    end
    #One extra iteration to ensure that all closest points are in Cv (tollerance)
    Cv = close_vectors(L,Vector(Rational.(temp)),j)
    vec = [Cv[i][2] for i=1:length(Cv)]
    min = findmin(vec)[1]
    close=Set{}()
    for i=1:length(Cv)
        if Float64.(abs(vec[i].-min)) <= 10^-5
            push!(close,transpose(C.U*Vector(Float64.(Cv[i][1]))))
        end
    end
    return close
end

"""
    Random_G_psd(G,d)

Generates a random PQF Q invariant under Gf
"""
function Random_G_psd(G::MatrixGroup,d::Int)
    Gen_set = invariant_symmetric_forms(G)
    #println(Gen_set)
    while true
        random = rand(-4:4,length(Gen_set))
    
        M = Oscar.matrix(QQ,Rational.(zeros(d,d)))

        for i=1:length(Gen_set)
            M = M + Gen_set[i].*random[i]
        end
        if minimum(eigvals(Matrix(Float64.(M)))) >0 && det(M) != 0
            return M
        end
    end
end

"""
    TrigidityIndex(D)

compares the T-rigidity Index of a Secondary cone C and dim on a subspace dim(S)
"""
function TrigidityIndex(C::Polyhedra.Polyhedron,S::Polyhedra.Polyhedron)
    println("Dimension closure T-secondary cone: ",Polyhedra.dim(C))
    println("Dimension T: ",Polyhedra.dim(S))
    return Polyhedra.dim(C) == Polyhedra.dim(S)
end

function ValidFlipCheck(L::SecondaryCone,L1::SecondaryCone)
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

    if length(L1.Inequalities) !=0 && length(L1.Equalities) !=0
        Q1 = collect(typeof(first(keys(L1.Inequalities))),keys(L1.Inequalities))
        H1 = collect(typeof(first(keys(L1.Equalities))),keys(L1.Equalities))
        P1 = Polyhedra.polyhedron(hrep(Q1) ∩ hrep(H1), CDDLib.Library(:exact))
    elseif length(L1.Inequalities) !=0 && length(L.Equalities) ==0
        Q1 = collect(typeof(first(keys(L1.Inequalities))),keys(L1.Inequalities))
        P1 = Polyhedra.polyhedron(hrep(Q1), CDDLib.Library(:exact))
    elseif length(L1.Inequalities) ==0 && length(L1.Equalities) !=0
        H1 = collect(typeof(first(keys(Delta.Equalities))),keys(L1.Equalities))
        P1 = Polyhedra.polyhedron(hrep(H1), CDDLib.Library(:exact))
    end

    return Polyhedra.dim(P) == Polyhedra.dim(P1)
end

"""
    affine_independent_set()

    for a given set of points return a affinly independet n-point set containt init
"""
function affine_independent_set(v::Set{},d::Int,n::Int)
    c = collect(v)
    counter = 0
    while true
        s = sample(1:length(c), n, replace = false)
        mat = ones(n,d+1)
        i=1
        for t in s
            mat[i,1:d] = c[t]
            i+=1
        end
        U, S, V = svd(mat)
        tol = 10 * eps(maximum(S))
        r = sum(S .> tol)
        if r == n
            return c[s]
        end
        #The counter is a filter for lower dimensional faces. The enumeration sill findes Polytope Pairs
        #,that share a lower dimensional face. Yet there cant be d affine independent vertices.
        #So just Ignore them. (propably a smater way so that you dont enumerate them in the fist place IDK)
        #if counter == 10^d
        #    return 0
        #end
        if counter == 500
            return 0
        end
        counter +=1
    end
end


"""
    sphere(c,d)

Calculates the d-circumsphere of the the points v using the Cayley-Menger matrix CM.
"""
function sphere(v::Set{}, d::Int64)
    CM1=zeros(d+1,d+1)
    #find d+1 affine independent points in d-space
    p=affine_independent_set(v,d,d+1)
    for (i,p1) in enumerate(p)
        for (j,k1) in enumerate(p)
            CM1[i,j] = norm(p1-k1,2)^2
        end
    end
    up = ones(d+2)
    up[1] = 0
    left = ones(d+1)
    CM = hcat(left,CM1)
    CM = vcat(up',CM)
    inv_CM = inv(CM)
    scale = inv_CM[:,1]
    c = zeros(d)
    for (i,p1) in enumerate(p)
        c = c.+scale[i+1]*p1'
    end
    return c, sqrt(-inv_CM[1,1]/2)
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
            v[k] = M[i,j]
            k = k+1
        end
        v[k] = M[i,i]
        k = k+1
    end
    return v
end

"""
    mat2vec(M)

Return a a vector of size n*(n+1)/2 containing the rows of the lower triangular part of `M`.

See also [`vec2mat`](@ref)
"""
function mat2vecdouble(M)
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




########### Functions ###########
"""
    TSecondaryCone(Delta,G,d)

Converts the T-secondary cone Polyhedron in a TSecondaryCone.
"""
function TSecondaryCone(P::Polyhedra.Polyhedron,L::SecondaryCone, D::DeloneSub)
    
    Facets = Dict{HalfSpace, Set{Tuple{DeloneCell, DeloneCell}}}()
    for (Halfspace,SetSimplicies) in L.Inequalities
        Halfspace ∈ halfspaces(P) ? Facets[Halfspace] = SetSimplicies : continue
    end

    Eq = Set{HyperPlane}()
    Hp = collect(Polyhedra.hyperplanes(P)) 
    for (i,H) in enumerate(Hp)
        push!(Eq,H)
    end
    return TSecondaryCone(Polyhedra.dim(P),D,Facets,Eq)
end


"""
    createTSecondaryCone(Delta,G,d)

Creates the T-secondary cone as a Polyhedron.
"""
function createTSecondaryCone(Delta::SecondaryCone,G::MatrixGroup,d::Int)
    if length(Delta.Inequalities) !=0 && length(Delta.Equalities) !=0
        Q = collect(typeof(first(keys(Delta.Inequalities))),keys(Delta.Inequalities))
        H = collect(typeof(first(keys(Delta.Equalities))),keys(Delta.Equalities))
        P = Polyhedra.polyhedron(hrep(Q) ∩ hrep(H), CDDLib.Library(:exact))
    elseif length(Delta.Inequalities) !=0 && length(Delta.Equalities) ==0
        Q = collect(typeof(first(keys(Delta.Inequalities))),keys(Delta.Inequalities))
        P = Polyhedra.polyhedron(hrep(Q), CDDLib.Library(:exact))
    elseif length(Delta.Inequalities) ==0 && length(Delta.Equalities) !=0
        H = collect(typeof(first(keys(Delta.Equalities))),keys(Delta.Equalities))
        P = Polyhedra.polyhedron(hrep(H), CDDLib.Library(:exact))
    else
        return 0
    end

    Sub = invariant_symmetric_forms(G)
    M = zeros(length(Sub)+1,Int(d*(d+1)/2))
    for i =1:length(Sub)
        M[i,:] = mat2vec(Matrix(Float64.(Sub[i])))
    end

    P1 = Polyhedra.polyhedron(Polyhedra.vrep(M), CDDLib.Library(:exact))

    h = collect(Polyhedra.hyperplanes(hrep(P1)))

    println("Dimension secondary cone: ",Polyhedra.dim(P))

    for i=1:length(h)
        P = Polyhedra.polyhedron(hrep(P) ∩ h[i], CDDLib.Library(:exact))
    end

    removehredundancy!(P)
    return P, P1
end

"""
    SecondaryCone(D::DeloneSub)

Create the `SecondaryCone` for a given `Delone` Subdivision.
"""
function SecondaryCone(D::DeloneSub)
    d = D.Dimension
    Inequalities = Dict{HalfSpace{Rational{BigInt},Vector{Rational{BigInt}}},Set{Tuple{DeloneCell,DeloneCell}}}()
    for (S,listS) in adjacent(D), T in listS
        H = Polyhedra.HalfSpace(-mat2vecdouble(regulator(S,T,D.Dimension)),0)
        # If the Inequality exits append the Tuple `(S,T)` if not create a new Set.
        multiple = checkmultiple(Inequalities, H)
        if multiple == false
            Inequalities[H] = Set([(S,T)])
        else
            push!(Inequalities[Polyhedra.HalfSpace(multiple,0)],(S,T))
        end
    end

    Equalities = Dict{HyperPlane{Rational{BigInt}, Vector{Rational{BigInt}}}, Set{DeloneCell}}()
    C = collect(D.Cells)
    for (S,Cell) in enumerate(C)
        W = collect(Cell.Vertices)
        P1 = affine_independent_set(Cell.Vertices,d,d+1)
        for (i,w) in enumerate(W)
            H = Polyhedra.HyperPlane(mat2vecdouble(eq_regulator(P1,w,D.Dimension)),0)

            haskey(Equalities, H) ? push!(Equalities[H],(Cell)) : Equalities[H] = Set([(Cell)])
        end
    end
    return SecondaryCone(D.Dimension,D,Inequalities,Equalities)
end

""" 
    reduce(L::SecondaryCone)

Create a new `SecondaryCone` w/o redundent inequalities.
"""
function reduce(L::SecondaryCone)
    Q = collect(typeof(first(keys(L.Inequalities))),keys(L.Inequalities))
    H = collect(typeof(first(keys(L.Equalities))),keys(L.Equalities))
    P = Polyhedra.polyhedron(hrep(Q) ∩ hrep(H), CDDLib.Library(:exact))
    removehredundancy!(P)
    #Checks if two halfspaces yield the same inequality when projected on to the hyperplanes
    #and adds their respective respective repartitioning polytopes to the set.
    for (Halfspace,SetSimplicies) in L.Inequalities
        for (Halfspace1,SetSimplicies1) in L.Inequalities
            pole1 = Polyhedra.polyhedron(Halfspace1 ∩ hrep(H) ∩ Halfspace, CDDLib.Library(:exact))
            removehredundancy!(pole1)
            if length(halfspaces(pole1)) == 1
                L.Inequalities[Halfspace]  = L.Inequalities[Halfspace] ∪ SetSimplicies1
            end
        end
    end
    Facets = Dict{HalfSpace, Set{Tuple{DeloneCell, DeloneCell}}}()
    for (Halfspace,SetSimplicies) in L.Inequalities
        Halfspace ∈ halfspaces(P) ? Facets[Halfspace] = SetSimplicies : continue
    end

    solids = Dict{HyperPlane, Set{DeloneCell}}()
    for (Hyperplane,Simplicies) in L.Equalities
        Hyperplane ∈ Polyhedra.hyperplanes(P) ? solids[Hyperplane] = Simplicies : continue
    end
    return SecondaryCone(L.Dimension, L.DeloneSub,Facets,solids)
    
end

"""
        Checkmultiple()

    checks if a multiple of the halfspace H already exists in the Inequality set.
"""
function checkmultiple(Inequalities,H::Polyhedra.HalfSpace)
    if length(Inequalities) == 0
        return false
    end
    Q = collect(keys(Inequalities))
    for (i,Hs) in enumerate(Q)
        a=0
        for j = 1:length(H.a)
            if H.a[j] != 0 && Hs.a[j] != 0
                a = H.a[j] / Hs.a[j]
                if a!=0
                    break
                end
            end
        end
        #if a <= 0
        #    continue
        #end
        if a.* Hs.a == H.a
            return copy(Hs.a)
        end
    end
    return false
end

"""
    startVertex(Q)

Calculates the Voronoi Cell and the start vertex.[DSV, Coplexity and algorithms for computing voronoi cells of lattices](Algorithm 3.3)
"""
function startVertex(Q::QQMatrix,d::Int64)
    C = cholesky(Matrix(Float64.(Q)))
    c = rand!(zeros(d))
    B = Set([])
    b = []
    for i = 1:d
        push!(B,transpose(C.U[:,i]))
        push!(B,-transpose(C.U[:,i]))
    end
    
    while true
        Ɛ = Set([])
        l = length(B)
        A = zeros(l,d)
        b = zeros(l)
        B1 = collect(B)
        for (i,b1) in enumerate(B1)
            b1[abs.(b1) .< 10^-5] .= 0
            A[i,:] = Rational.(b1)
            b[i] = Rational.(0.5*(dot(b1,b1)))
        end
        P = Oscar.polyhedron(A,b)
        LP = linear_program(P,c,k = 0, convention = :max)
        x = optimal_vertex(LP) 
        Ɛ = closest_vectors(Q,Vector(Float64.(x)))
        zeros(d)' ∈ (Ɛ) && (return Vector(Float64.(x)),Ɛ)
        B = B ∪ Ɛ
    end
end

"""
    getFacets(D,d)

Returns all the facets of the Delone cell.
"""
function getFacets(D::DeloneCell,d::Int64)
    D1 = collect(D.Vertices)
    D_Mat = zeros(D.Amount,d)
    for i =1:D.Amount
        D_Mat[i,:] = D1[i]
    end
    P = Oscar.convex_hull(D_Mat)
    A,b = halfspace_matrix_pair(facets(P))
    F = Set{Facet}()
    for i = 1:length(b)
        #Get all the Vertices in F
        V = Set{}()
        for (j,v) in enumerate(D1)
            if abs(dot(Vector(Float64.(A[i,:])),v)-b[i]) <= 10^-5
                push!(V,v)
            end
        end
        #Facet containing a halfspace and the vertices
        push!(F,Facet(A[i,:],b[i],V))
    end
    return Facets(length(b),F)
end

"""
    getNeighbor(F,D)

Calculates the Delone cell D' which shares the facet F with D.[DSV, Coplexity and algorithms for computing voronoi cells of lattices](Algorithm 3.4)
"""
function getNeighbor(F::Facet,Q::QQMatrix,d::Int64)
    L = integer_lattice(gram = Q)
    C = cholesky(Matrix(Float64.(Q)))
    vert = collect(F.Vertices)
    v = vert[1]
    i=1
    #Calculate a vertex on the other side of the hyperplane defined by the facets.
    while dot(F.a,v) - F.b <= 10^-5
        #Note that you have to use the integer lattice in the oscar.close_vectors function for that mulitiply the lattice vector by the inverse of the basis matrix.
        temp = inv(C.U)*vert[1]'
        temp[abs.(temp) .< 10^-5] .= 0
        clo = close_vectors(L,Vector(Rational.(temp)),i)
        for j = 1:length(clo)
            if dot(F.a,C.U*clo[j][1])-F.b > 10^-5
                v = Vector(Float64.(C.U*clo[j][1]))
            end
        end
        i=i+1
    end
    #While loop of the algorithm
    while true
        S = copy(F.Vertices)
        push!(S,transpose(v))
        c,r = sphere(S,d)
        V_bar = closest_vectors(Q,c)
        V_bar1 = collect(V_bar)
        for (i,v_iter) in enumerate(V_bar1)
            if abs(norm(v_iter'-c,2) - r) <= 10^-5
                return V_bar
            end
        end
        #pop!(S,transpose(v))
        v = copy(transpose(V_bar1[1]))
    end
end


"""
    getSubdivision(Q) 
      
Calculates the Triangulation for a given positive definite quadratic form. [DSV, Coplexity and algorithms for computing voronoi cells of lattices](Main algorithm)
"""
function getSubdivision(Q::QQMatrix,d::Int64)
    V = Set([])
    iter = 0
    while true
        iter += 1
        x,V = startVertex(Q,d)
        length(V) >= d+1 && break
        if iter >= 1000
            return 0
        end
    end
    T = Queue{DeloneCell}()
    enqueue!(T,DeloneCell(length(V),V))
    M = Queue{DeloneCell}()
    facet_size = 0

    while !isempty(T)
        D = dequeue!(T)
        enqueue!(M,D)
        ℱ = getFacets(D,d)
        ℱ1 = collect(ℱ.Facets)
        for (i,F) in enumerate(ℱ1)
            facet_size = length(F.Vertices)
            D_bar = getNeighbor(F,Q,d)
            if zeros(d)' ∈ D_bar && !any([D_bar == D1.Vertices for D1 ∈ T ∪ M])
                enqueue!(T,DeloneCell(length(D_bar),D_bar))
            end
        end
    end
    return DeloneSub(T ∪ M,d,facet_size)
end

"""
    DeloneSub(M,d)

    Creates a Delone Subdivision out of a Queue of Delone cells
"""
function DeloneSub(M::Vector{DeloneCell},d::Int64,facet_size::Int64)

    l = length(M)
    D = DeloneSub(d,l,facet_size,Set{DeloneCell}())
    for (i,Cell) in enumerate(M)
        push!(D.Cells, Cell)
    end

    return D
end


"""
    regulator(S1,S2)

Return the Regulator of the adjacent DeloneCells `S1` and `S2`.
"""
function regulator(S1::DeloneCell,S2::DeloneCell,d::Int64)
    # Create Arrays for a fixed ordering
    dimension = d
    P1 = collect(S1.Vertices)
    P2 = collect(S2.Vertices)
    #find d affine indepentent points in a d-1-facet in d-space
    facet_verts = affine_independent_set(S1.Vertices ∩ S2.Vertices,d,d)

    #If there are no affine independent points in the set that means that the two vertices only share
    #a face not a facet
    if facet_verts == 0
        return zeros(dimension,dimension)
    end

    w_1 = rand(setdiff(P1,P2))
    w_2 = rand(setdiff(P2,P1))

    # Calculate an affine relationship alpha
    B = Array{Rational{BigInt}}(undef, dimension+1, dimension+1)
    for i = 1:dimension
        for j = 1:dimension
            B[j,i] = facet_verts[i][j]
        end
    end
    for i = 1:dimension
        B[i,dimension+1] = w_1[i]
    end
    for i in 1:dimension+1
        B[dimension+1,i] = 1
    end
    v = append!(copy(vec(w_2)),1)

    alpha = B\v

    # Compute regulator
    reshape(w_2, length(w_2), 1)
    N = w_2*w_2'
    temp = zeros(dimension, dimension)
    for i in 1:dimension+1
        v =  B[1:dimension,i]
        reshape(v, length(v), 1)
        temp += alpha[i]*v*v'
    end
    Matrix1 = N-temp
    Matrix1[abs.(Matrix1) .< 10^-5] .= 0.0
    return round.(Matrix1)
end


"""
    eq_regulator(Cell,w,d)

Return the Regulator for the equalities.
"""
function eq_regulator(Vert,w,d)
    dimension = d
    P1 = copy(Vert)

    B = Array{Rational{BigInt}}(undef, dimension+1, dimension+1)
    for i = 1:dimension+1
        for j = 1:dimension
            B[j,i] = P1[i][j]
        end
    end
    
    for i in 1:dimension+1
        B[dimension+1,i] = 1
    end
    v = append!(copy(vec(w)),1.0)
    alpha = B\v

    # Compute regulator
    reshape(w, length(w), 1)
    N = w*w'
    temp = zeros(dimension, dimension)
    for i in 1:dimension+1
        v =  B[1:dimension,i]
        reshape(v, length(v), 1)
        temp += alpha[i]*v*v'
    end
    Matrix1 = N-temp
    Matrix1[abs.(Matrix1) .< 10^-5] .= 0.0
    return round.(Matrix1)
end


"""
    adjacent(S::DeloneCell,D::DeloneSub)

Return a Set of all Simplicies `T` that are adjacent over a common facet with `S`.
"""
function adjacent(S::DeloneCell, D::DeloneSub)
    Inzident = Set{DeloneCell}([])
    #The Problem is that is is not sufficent to check the translates on S for beeing neighboring 
    #because the neighbor does not need to be the same 
    for Cell in D.Cells
        A = intersect(S.Vertices,Cell.Vertices)
        if A == S.Vertices && A == Cell.Vertices
            continue
        end
        #If two cells have a d-1 dimensional intersection they have to share a facet no?
        #checks if the polytopes are ajacent. It is not sufficent to  only check for amount of verticies because
        #lower dimensional faces could have the same amount of verticies. thats why we chek for affine independnece
        if length(A) >=D.Dimension
            B = affine_independent_set(A,D.Dimension,D.Dimension)
            if B != 0
                push!(Inzident,deepcopy(Cell))
            end
        end
    end
    return Inzident
end

"""
    adjacent(D::Delone)

Return a set of tuples `(S,adj)` where `S` is a Simplex in `D` and `adj` contains all Simplicies `T` that are adjacent over a common facet with `S`.
"""
function adjacent(D::DeloneSub)
    n = nthreads()  # Get number of threads
    adj_chunks = [Set{Tuple{DeloneCell, Set{DeloneCell}}}() for _ in 1:n]  # Each thread gets its own Set

    cole = collect(D.Cells)  # Convert Set to a Vector

    @threads for i in eachindex(cole)
        S = cole[i]
        Set_adjacent_S = adjacent(S, D)
        push!(adj_chunks[threadid()], (S, Set_adjacent_S))  # No lock needed
    end

    return union(adj_chunks...)  # Merge all results at the end

    return adj
end

"""
    TgenericDelone(T)

Calculates T-generic start triangulation.[DSV, A generalization of Voronois reduction theory and its applications](Algorithm 2)
"""
function TgenericDelone(G::MatrixGroup,d::Int64)
    mat1 = Oscar.matrix(QQ,Rational.(zeros(d,d)))
    D_1 = DeloneSub(0,0,0,Set{DeloneCell}())
    bool = true
    while bool
        mat1 = Random_G_psd(G,d)
        D = getSubdivision(mat1,d)
        if D == 0
            continue
        end
        D_1 = toIntegerLattice(D,mat1)
        cone = SecondaryCone(D_1)

        P,P1=createTSecondaryCone(reduce(cone),G,d)
        TrigidityIndex(P,P1) && (bool = false)
        println(" ")
    end
    return mat1,D_1
end

end #End of Module