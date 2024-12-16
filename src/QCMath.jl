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


