########### Imports ###########

using Permutations
using Polyhedra
using CDDLib
using LinearAlgebra
using DataStructures
using Oscar


include("./scc.jl")
include("./Tscc.jl")
include("./Maxdet_scc.jl")
include("./Maxdet_Tscc.jl")

using .scc
using .Tscc
using .Maxdet_scc
using .Maxdet_Tscc




########### Functions ###########

"""
    start(G,d)

Starts the enumeration process
"""
function start(G::MatrixGroup,d::Int64)
    if isfinite(G) == false
        println(" ")
        println("Error: The group is infinite")
    elseif length(G) == 1
        println(" ")
        println("You have entered the trivial group!")
        println("Starting the enumeration of inequivalent secondary cones with T = S^($d) in Dimension $d ...")
        println(" ")
        time = @elapsed (cones = scc.Voronoi(d))
        println(length(cones)," inequivalent secondary cones have been found in, ",time," seconds.")
        println("Starting Maxdet ...")
        Tuple = collect(Maxdet_scc.solve(cones,d))
        println("Optimal solution found!")
        println("Normalized covering density: ≈",Tuple[1])
        println("PQF: ")
        display(Tuple[2])
    else
        println(" ")
        println("You have entered a group of size ",length(G),"!")
        println("Starting the enumeration of inequivalent T-secondary cones in Dimension $d ...")
        println(" ")
        time = @elapsed (cones = Tscc.GeneralVoronoi(G,d))
        println(" ")
        println(length(cones)," inequivalent T-secondary cones have been found in, ",time," seconds.")
        println("Starting Maxdet ...")
        Tuple = collect(Maxdet_Tscc.solve(cones,d))
        println("Optimal solution found!")
        println("Normalized covering density: ≈",Tuple[1])
        println("PQF: ")
        display(Tuple[2])
    end
end

let 
    PATH_ISOM = "ISOM_and_AUTO/"
    ISOM_C = PATH_ISOM * "ISOM.c"
    ISOM_EXE = PATH_ISOM * "ISOM"
    if !isfile(ISOM_EXE)
        run(`gcc -o $ISOM_EXE $ISOM_C ISOM_and_AUTO/auttools.c ISOM_and_AUTO/bachtools.c ISOM_and_AUTO/iotools.c ISOM_and_AUTO/isotools.c ISOM_and_AUTO/lattools.c ISOM_and_AUTO/mattools.c ISOM_and_AUTO/orbtools.c ISOM_and_AUTO/preproc.c ISOM_and_AUTO/sorttools.c`)
    end

    #Set the dimension
    d = 4

    #Scc input

    #mat = Oscar.matrix(ZZ, I(d))






    #Gscc input

    #Eisenstein group for d divisable by 2

    #=Id = Diagonal(ones(Int(d/2)))
    eisenstein_mat = [0 -1; 1 -1]   
    mat = Oscar.matrix(QQ, Rational.(kron(Id,eisenstein_mat)))=#



    #Gaussian group for d divisable by 2

    Id = Diagonal(ones(Int(d/2)))
    gaussian_mat = [0 -1; 1 0]
    mat = Oscar.matrix(QQ, Rational.(kron(Id,gaussian_mat)))
    


    #Hurwitz quaternionic group 2*A4 only for d divisable by 4
    
    #This gives a nicer lattice basis

    #=Id = Diagonal(ones(Int(d/4)))
    twoA4 = Set{QQMatrix}()
    push!(twoA4, Oscar.matrix(QQ, Rational.(kron(Id,[0 0 1 0; 1 1 1 -2; -1 0 0 0; 0 1 1 -1]))))
    push!(twoA4, Oscar.matrix(QQ, Rational.(kron(Id,[0 0 0 1; 0 1 1 -1; -1 -1 0 1; -1 0 0 1]))))
    mat = collect(twoA4)=#

    #Described in the paper

    #=Id = Diagonal(ones(Int(d/4)))
    twoA4 = Set{QQMatrix}()
    push!(twoA4, Oscar.matrix(QQ, Rational.(kron(Id,[-1 -1 -1 -1; 1 1 0 1; 1 1 1 0; 1 0 1 1]))))
    push!(twoA4, Oscar.matrix(QQ, Rational.(kron(Id,[-1 -2 0 0; 1 1 0 0; 0 1 0 -1; 1 1 1 0]))))
    mat = collect(twoA4)=#

    G = matrix_group(mat)
    start(G,d)

end
