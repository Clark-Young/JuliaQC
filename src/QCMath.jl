# Description: This Julia code is a submodule of JuliaQC, here we write some algorithms for quantum chemistry calculation.
# Author: Chunyu Yang
# Email: chunyu.yang@duke.edu, chyyangustc@outlook.com (permanent)
# Organization: Department of Chemistry, Duke University
# Date: 2024-12-15
# Reference: Szabo and Ostlund's Modern Quantum Chemistry; Helgaker's Molecular Electronic Structure Theory.
#= Here we mainly use
- Davidson Diagonalization and Jocobi-Davidson Diagonalization
- Romberg Integration
=#


using LinearAlgebra

M = [1 2 3 4 ; 4 3 2 1 ; 1 2 3 4 ; 4 3 2 1]

function DavidsonDiag(A::Matrix, m::Int, epsilon::Float64=1e-12, maxIter::Int=1000)
    n = size(A)[1]
    if m > n
        println("m should be less than n")
        return
    end
    
    
end

function JacobiDavidsonDiag(A::Matrix, m::Int, epsilon::Float64=1e-12, maxIter::Int=1000)
    
end

function JacobiDiag(A::Matrix, epsilon::Float64=1e-12, maxIter::Int=1000)
    
end


