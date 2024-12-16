#!/usr/bin/env julia

using LinearAlgebra

struct atomList
    atom::String
    x::Float64
    y::Float64
    z::Float64
    index::Int
end

function symmetricOrthogonalization(S::Matrix)
    eigval, eigvec = eigen(S)
    X = eigvec * diagm(0 => 1 ./ sqrt.(eigval)) * eigvec'
    return X
end

function canonicalOrthogonalization(S::Matrix)
    eigval, eigvec = eigen(S)
    X = eigvec * diagm(0 => 1 ./ sqrt.(eigval))
    return X
    
end

function totalEnergy(P::Matrix, h_core::Matrix, Fock::Matrix)
    return sum(P .* (h_core + Fock))/2
    
end

function SCF(h_core::Matrix, S::Matrix, integral2e::Array{Float64, 4}, N_e::Int, max_iter::Int, tol::Float64)
    n = size(h_core, 1)
    Fock = h_core
    P = ones(n, n)
    for i in 1:max_iter
        # println(P)
        G = zeros(n, n)
        for i in 1:n
            for j in 1:n
                for k in 1:n
                    for l in 1:n
                        G[i, j] += P[k, l] * (integral2e[i, j, k, l] - 0.5 * integral2e[i, l, k, j])
                    end
                end
            end
        end

        Fock = h_core + G
        # println(Fock)
        X = symmetricOrthogonalization(S)
        # X = canonicalOrthogonalization(S)

        FockPrime = X' * Fock * X
        eigval, eigvec = eigen(FockPrime)
        println("Iteration: ", i, " Energy: ", totalEnergy(P, h_core, Fock))
        # println(eigvec)
        C = X * eigvec
        # println(C)
        P_new = 2 * C[:, 1:N_e÷2] * C[:, 1:N_e÷2]'
        # println(P_new)
        # println(P)
        if norm(P_new - P) < tol
            return totalEnergy(P_new, h_core, Fock), P_new, eigval
        end
        P = copy(P_new)
        # println("Iteration: ", i, " Energy: ", eigval[1])
    end
end

function MullikenPopulationAnalysis(atomlist::Vector,atom2basis::Vector,element2Z::Dict, P::Matrix,S::Matrix)
    n = size(P, 1)
    N = size(atomlist, 1)
    MPA = zeros(N)
    PS= P*S
    for i in atomlist
        index=atom2basis[i.index]
        # println(index)
        for j in index
            MPA[i.index] += PS[j,j]
        end
        MPA[i.index] = element2Z[i.atom] - MPA[i.index]
    end
    return MPA
end

element2Z = Dict("H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10)
atom2basis = [[1],[2]]
atomlist = [atomList("H", 0.0, 0.0, 0.0, 1), atomList("H", 0.0, 0.0, 0.74, 2)]

h_core = [-1.1204 -0.9584; -0.9584 -1.1204]
S = [1.0 0.6593; 0.6593 1.0]

integral2e=zeros(Float64, 2, 2, 2, 2)
integral2e[1, 1, 1, 1] = 0.7746
integral2e[2, 2, 2, 2] = 0.7746
integral2e[1, 1, 2, 2] = 0.5697
integral2e[2, 2, 1, 1] = 0.5697
integral2e[1, 2, 1, 2] = 0.2970
integral2e[2, 1, 2, 1] = 0.2970
integral2e[1, 2, 2, 1] = 0.2970
integral2e[2, 1, 1, 2] = 0.2970
integral2e[1, 2, 2, 2] = 0.4441
integral2e[2, 1, 2, 2] = 0.4441
integral2e[2, 2, 1, 2] = 0.4441
integral2e[2, 2, 2, 1] = 0.4441
integral2e[1, 1, 1, 2] = 0.4441
integral2e[1, 1, 2, 1] = 0.4441
integral2e[1, 2, 1, 1] = 0.4441
integral2e[2, 1, 1, 1] = 0.4441
E0,P,ε=SCF(h_core, S, integral2e, 2, 100, 1e-8)
println("The total energy of H2 is: ", E0, " Hartree.")
MPA=MullikenPopulationAnalysis(atomlist,atom2basis,element2Z, P, S)
println("The atomic charge of on H1 is: ", MPA[1], "; and the atomic charge of on H2 is: ", MPA[2], ".")

h_core_HeH = [-2.6527 -1.3472;  -1.3472 -1.7318]
S_HeH = [1.0 0.4508; 0.4508 1.0]
integral2e_HeH=zeros(Float64, 2, 2, 2, 2)
integral2e_HeH[1, 1, 1, 1] = 1.3072
integral2e_HeH[2, 2, 2, 2] = 0.7746
integral2e_HeH[1, 1, 2, 2] = 0.6057
integral2e_HeH[2, 2, 1, 1] = 0.6057
integral2e_HeH[1, 2, 1, 2] = 0.1773
integral2e_HeH[2, 1, 2, 1] = 0.1773
integral2e_HeH[1, 2, 2, 1] = 0.1773
integral2e_HeH[2, 1, 1, 2] = 0.1773
integral2e_HeH[1, 2, 2, 2] = 0.3118
integral2e_HeH[2, 1, 2, 2] = 0.3118
integral2e_HeH[2, 2, 1, 2] = 0.3118
integral2e_HeH[2, 2, 2, 1] = 0.3118
integral2e_HeH[1, 1, 1, 2] = 0.4373
integral2e_HeH[1, 1, 2, 1] = 0.4373
integral2e_HeH[1, 2, 1, 1] = 0.4373
integral2e_HeH[2, 1, 1, 1] = 0.4373

atomlist_HeH = [atomList("He", 0.0, 0.0, 0.0, 1), atomList("H", 0.0, 0.0, 1.4632, 2)]
# SCF(h_core_HeH, S_HeH, integral2e_HeH, 2, 100, 1e-8)
E0,P_HeH,ε=SCF(h_core_HeH, S_HeH, integral2e_HeH, 2, 100, 1e-8)
MPA_HeH=MullikenPopulationAnalysis(atomlist_HeH,atom2basis,element2Z, P_HeH, S_HeH)

println("The total energy of HeH is: ", E0, " Hartree.")
println("The atomic charge of on He is: ", MPA_HeH[1], "; and the atomic charge of on H is: ", MPA_HeH[2], ".")



