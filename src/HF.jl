# Description: This Julia Code is used to perform Hartree-Fock (HF) calculation.
# Author: Chunyu Yang
# Email: chunyu.yang@duke.edu, chyyangustc@outlook.com (permanent)
# Organization: Department of Chemistry, Duke University
# Date: 2024-12-15
# Reference: Szabo and Ostlund's Modern Quantum Chemistry;
#= Here we mainly use
- Roothaan-Hall equation
- Pople-Nesbet equation
- Mulliken Population Analysis
=#

module HF

export atomList, SCF, MullikenPopulationAnalysis, element2Z

include("Integral.jl")
using .Integral

using LinearAlgebra

element2Z = Dict("H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10, "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20, "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, "Mn" => 25, "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30, "Ga" => 31, "Ge" => 32, "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36, "Rb" => 37, "Sr" => 38, "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42, "Tc" => 43, "Ru" => 44, "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48, "In" => 49, "Sn" => 50, "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54, "Cs" => 55, "Ba" => 56, "La" => 57, "Ce" => 58, "Pr" => 59, "Nd" => 60, "Pm" => 61, "Sm" => 62, "Eu" => 63, "Gd" => 64, "Tb" => 65, "Dy" => 66, "Ho" => 67, "Er" => 68, "Tm" => 69, "Yb" => 70, "Lu" => 71, "Hf" => 72, "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77, "Pt" => 78, "Au" => 79, "Hg" => 80, "Tl" => 81, "Pb" => 82)

struct atom
    atomName::String
    x::Float64
    y::Float64
    z::Float64
    index::Int
end

struct molecule
    atomlist::Vector{atom}
    charge::Int
    multiplicity::Int
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

#= The procedure of Roothaan-Hall SCF calculation:
1. Specify a molecule (atom, geometry, charge, multiplicity, and basis set)
2. Calculate all required molecular integrals (overlap, kinetic, nuclear attraction, electron repulsion)
3. Diaognalize the overlap matrix S to get X (symmetric orthogonalization or canonical orthogonalization)
4. Set the initial guess of the density matrix P
5. Calculate the G matrix and Fock matrix (Fock = h_core + G)
6. Transform the Fock matrix to the orthogonal basis FockPrime = X' * Fock * X
7. Diagonalize the Fock matrix to get the orbital energies (ε) and coefficients (C')
8. Calculate C = X * C', and then the new density matrix P_new
9. If the difference between P_new and P is less than a certain threshold, return the total energy and the density matrix
10. Otherwise, update P and repeat the above steps
=#

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

function MullikenPopulationAnalysis(ml::molecule,atom2basis::Vector,element2Z::Dict, P::Matrix,S::Matrix)
    n = size(P, 1)
    N = size(atomlist, 1)
    MPA = zeros(N)
    PS= P*S
    for i in ml.atomlist
        index=atom2basis[i.index]
        # println(index)
        for j in index
            MPA[i.index] += PS[j,j]
        end
        MPA[i.index] = element2Z[i.atom] - MPA[i.index]
    end
    return MPA
end


end
