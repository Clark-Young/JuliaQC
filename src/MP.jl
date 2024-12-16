# Description: This Julia Code is used to calculate the second order Moller-Plesset perturbation theory (MP2) energy.
# Author: Chunyu Yang
# Email: chunyu.yang@duke.edu, chyyangustc@outlook.com (permanent)
# Organization: Department of Chemistry, Duke University
# Date: 2024-12-13
# Reference: Szabo and Ostlund's Modern Quantum Chemistry; 

module MP
include("HF.jl")
include("Integral.jl")
using .HF
using .Integral


function MP2(atomlist::atomList, basis::String, N_e::Int, max_iter::Int, tol::Float64)
    n = size(atomlist, 1)
    S, h_core, integral2e = Integral.integrals(atomlist, basis)
    X = HF.symmetricOrthogonalization(S)
    G = zeros(n, n, n, n)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                for l in 1:n
                    G[i, j, k, l] = sum(integral2e[i, j, k, l] .* X)
                end
            end
        end
    end
    E = HF.SCF(h_core, S, integral2e, N_e, max_iter, tol)
    E_MP2 = 0
    for i in 1:n
        for j in 1:n
            for a in 1:n
                for b in 1:n
                    E_MP2 += G[i, j, a, b] * (2 * G[i, j, a, b] - G[i, b, a, j]) / (E[i] + E[j] - E[a] - E[b])
                end
            end
        end
    end
    return E_MP2 + E
    
end



end