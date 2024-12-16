# Description: This Julia Code is used to calculate the integrals of the Gaussian Type Orbitals
# Author: Chunyu Yang
# Email: chunyu.yang@duke.edu, chyyangustc@outlook.com (permanent)
# Organization: Department of Chemistry, Duke University
# Date: 2024-12-13
# Reference: Wim Klopper's Lecture Notes on esqc.org; Szabo and Ostlund's Modern Quantum Chemistry; Helgaker's Molecular Electronic Structure Theory.
#= Here we mainly use
- GTO product theorem;
- Obara–Saika (OS) recurrence relations (for overlap and kinetic integral);
- $1/r_{c} = ∫_{-∞}^{+∞} \exp (- (r-r_{c})^2 t^2 \mathrm{d} t$;
- After recurrance relation, the four center integrals can become two center integrals.
=#


using LinearAlgebra
using SpecialFunctions
include("HF.jl")
using .HF

struct gto # Define the struct of Gaussian Type Orbital
    center::Tuple{Float64,Float64,Float64} # The center of the GTO
    exponent::Float64 # The exponent of the GTO
    power::Tuple{Int,Int,Int} # The power of the GTO
    contractedCoeff::Float64 # The contracted coefficient of the GTO
    # Primary constructor
    function gto(center::Tuple{Float64,Float64,Float64}, 
        exponent::Float64,
        power::Tuple{Int,Int,Int},
        contractedCoeff::Float64=1.0)
        new(center, exponent, power, contractedCoeff)
    end

    # Moving constructor
    function gto(gto1::gto, center::Tuple{Float64,Float64,Float64})
        new(center, gto1.exponent, gto1.power, gto1.contractedCoeff)
    end
end
# The gto is set to be immutable, to avoid the modification of the struct.


struct cgto
    center::Tuple{Float64,Float64,Float64}
    gtos::Array{gto,1}
    function cgto(center::Tuple{Float64,Float64,Float64}, gtos::Array{gto,1})
        new(center, gtos)
    end
    function cgto(cgto1::cgto, center::Tuple{Float64,Float64,Float64})
        gtos = Array{gto,1}(undef, length(cgto1.gtos))
        for i in 1:length(cgto1.gtos)
            gtos[i] = gto(cgto1.gtos[i], center)
        end
        return new(center, gtos)
    end
    
end


function Boys(m::Int64,T::Float64)::Float64 # Boys function, used in the calculation of the electron repulsion integral
    return T^(-(m+0.5)) * (gamma(m+0.5) - gamma(m+0.5, T))/2  
end



# Helper function to sum over an interval with a given step size, used in the Romberg integration
function funcSum(func::Function, start::Float64, step::Float64, stop::Float64)::Float64
    s = 0.0
    x = start
    while x <= stop
        s += func(x)
        x += step
    end
    return s
end


# Romberg integration, the function is used to calculate Boys function, when the distance is small, we need to use numerical integration instead of Gamma function.
function romberg(a::Float64, b::Float64, M::Int, func::Function, varepsilon::Float64)::Float64
    h = b - a
    r1 = Vector{Float64}(undef, M)
    r2 = similar(r1)

    # Initial estimate using trapezoidal rule
    r1[1] = (func(a) + func(b)) * h / 2

    for k in 1:M-1
        h_k = h / 2^k
        # Estimate the sum of the function evaluated at new points
        sum_a = sum(func, a + h_k/2:h_k:b-h_k)
        
        # Calculate the next level of Romberg integration
        r2[1] = 0.5 * (r1[1] + sum_a * h_k)

        for j in 1:k
            r2[j+1] = r2[j] + (r2[j] - r1[j]) / (4^j - 1)
        end

        if abs(r2[k+1] - r1[k]) < 0.001 * varepsilon
            return r2[k+1]
        end
        
        # Copy r2 to r1 for the next iteration
        copyto!(r1, r2)
    end

    return r2[end]
end



function distance(atom1::Tuple{Float64,Float64,Float64}, atom2::Tuple{Float64,Float64,Float64})
    return sqrt((atom1[1] - atom2[1])^2 + (atom1[2] - atom2[2])^2 + (atom1[3] - atom2[3])^2)
end

function gto_m1(gto1::gto, i::Int) # Create a new GTO with the power of i-th direction minus 1
    if i == 1
        return gto(gto1.center, gto1.exponent, (gto1.power[1] - 1, gto1.power[2], gto1.power[3]))
    elseif i == 2
        return gto(gto1.center, gto1.exponent, (gto1.power[1], gto1.power[2] - 1, gto1.power[3]))
    elseif i == 3
        return gto(gto1.center, gto1.exponent, (gto1.power[1], gto1.power[2], gto1.power[3] - 1))
    else
        println("Invalid i")
    end
end

function gto_p1(gto1::gto, i::Int) # Create a new GTO with the power of i-th direction plus 1
    if i == 1
        return gto(gto1.center, gto1.exponent, (gto1.power[1] + 1, gto1.power[2], gto1.power[3]))
    elseif i == 2
        return gto(gto1.center, gto1.exponent, (gto1.power[1], gto1.power[2] + 1, gto1.power[3]))
    elseif i == 3
        return gto(gto1.center, gto1.exponent, (gto1.power[1], gto1.power[2], gto1.power[3] + 1))
    else
        println("Invalid i")
    end
end

# gto_m1(gto1,4)

function cal_p(gto1::gto, gto2::gto) # Calculate the new center according to GTO product rule, see the usage of the p in the overlap integral,
    newcenter = zeros(3)
    for i in 1:3
        newcenter[i] = (gto1.center[i] * gto1.exponent + gto2.center[i] * gto2.exponent) / (gto1.exponent + gto2.exponent)
    end
    return Tuple(newcenter)

end

function normalizeConstant(gto1::gto) # Calculate the normalization constant of the GTO
    return (2 * gto1.exponent / π)^0.75 * sqrt((8 * gto1.exponent)^sum(gto1.power) * factorial(gto1.power[1]) * factorial(gto1.power[2]) * factorial(gto1.power[3]) / factorial(2 * gto1.power[1]) / factorial(2 * gto1.power[2]) / factorial(2 * gto1.power[3]))
end



gto1 = gto((2.0, 3.0, 1.0), 2.4, (0, 0, 2))
gto2 = gto((2.0, 3.0, 2.0), 1.8, (0, 2, 0))

function overlapIntegral(gto1::gto, gto2::gto)::Float64 # Calculate the overlap integral of two GTOs
    ζ = gto1.exponent + gto2.exponent
    ξ = gto1.exponent * gto2.exponent / ζ
    p = cal_p(gto1, gto2)
    integral = 0.0
    if sum(gto1.power) == 0 && sum(gto2.power) == 0
        integral = (π / ζ)^1.5 * exp(-ξ * distance(gto1.center, gto2.center)^2)
        return integral
    else
        for i in 1:3
            if gto1.power[i] != 0
                # println(i)

                integral += (p[i] - gto1.center[i]) * overlapIntegral(gto_m1(gto1, i), gto2)
                # println((p[i] - gto1.center[i]))
                if gto1.power[i] > 1
                    # println("power")
                    integral += 0.5 * (gto1.power[i] - 1) / ζ * overlapIntegral(gto_m1(gto_m1(gto1, i), i), gto2)
                    # println(overlapIntegral(gto_m1(gto_m1(gto1, i),i), gto2))
                end
                if gto2.power[i] != 0
                    integral += 0.5 * gto2.power[i] / ζ * overlapIntegral(gto_m1(gto1, i), gto_m1(gto2, i))
                end
            end
            if sum(gto1.power) == 0
                if gto2.power[i] != 0
                    integral += (p[i] - gto2.center[i]) * overlapIntegral(gto1, gto_m1(gto2, i))
                    if gto2.power[i] > 1
                        integral += 0.5 * (gto2.power[i] - 1) / ζ * overlapIntegral(gto1, gto_m1(gto_m1(gto2, i), i))
                    end
                end
            end
        end
        return integral
    end
end

function overlapIntegral(cgto1::cgto, cgto2::cgto)::Float64 # Calculate the overlap integral of two Contracted GTOs
    integral = 0.0
    for i in 1:length(cgto1.gtos)
        for j in 1:length(cgto2.gtos)
            integral += cgto1.gtos[i].contractedCoeff * cgto2.gtos[j].contractedCoeff * overlapIntegral(cgto1.gtos[i], cgto2.gtos[j])
        end
    end
    return integral

end

function overlapIntegral(basis::Array{cgto,1})::Array{Float64,2} # Calculate the overlap matrix of the basis set
    n = length(basis)
    S = zeros(n, n)
    for i in 1:n
        for j in 1:n
            S[i, j] = overlapIntegral(basis[i], basis[j])
        end
    end
    return S
    
end

function kineticIntegral(gto1::gto, gto2::gto) # Calculate the kinetic integral of two GTOs
    integral = 0.0
    for i in 1:3
        integral -= 2 * gto2.exponent * (2 * gto2.power[i] + 1) * overlapIntegral(gto1, gto2)
        integral += 4 * gto2.exponent * gto2.exponent * overlapIntegral(gto1, gto_p1(gto_p1(gto2, i), i))
        if gto2.power[i] > 1
            integral += gto2.power[i] * (gto2.power[i] - 1) * overlapIntegral(gto1, gto_p1(gto_m1(gto2, i), i))
        end
    end
    return integral
end

function kineticIntegral(cgto1::cgto, cgto2::cgto) # Calculate the kinetic integral of two Contracted GTOs
    integral = 0.0
    for i in 1:length(cgto1.gtos)
        for j in 1:length(cgto2.gtos)
            integral += cgto1.gtos[i].contractedCoeff * cgto2.gtos[j].contractedCoeff * kineticIntegral(cgto1.gtos[i], cgto2.gtos[j])
        end
    end
    return integral
end


function kineticIntegral(basis::Array{cgto,1})::Array{Float64,2} # Calculate the kinetic matrix of the basis set
    n = length(basis)
    T = zeros(n, n)
    for i in 1:n
        for j in 1:n
            T[i, j] = kineticIntegral(basis[i], basis[j])
        end
    end
    return T
    
end



function nuclearAttractionIntegral(gto1::gto, gto2::gto, center::Tuple{Float64,Float64,Float64},m::Int8)::Float64
    α = gto1.exponent 
    β = gto2.exponent
    ζ = gto1.exponent + gto2.exponent
    ξ = gto1.exponent * gto2.exponent / ζ
    p = cal_p(gto1, gto2)
    if sum(gto1.power) == 0 && sum(gto2.power) == 0
        T = ζ * distance(p, center)^2
        return Boys(m, T) * 2 * π / ζ
    else
        integral = 0.0
        for i in 1:3
            if gto1.power[i] != 0
                integral += β/ζ* (gto2.center[i] - gto1.center[i]) * nuclearAttractionIntegral(gto_m1(gto1, i), gto2, center, m)
                integral += (center[i] - gto1.center[i] - β/ζ* (gto2.center[i] - gto1.center[i])) * nuclearAttractionIntegral(gto_m1(gto1, i), gto2, center, m+1)
                if gto1.power[i] > 1
                    integral += 0.5 * gto1.power[i] / ζ * (nuclearAttractionIntegral(gto_m1(gto_m1(gto1, i), i), gto2, center, m) - nuclearAttractionIntegral(gto_m1(gto_m1(gto1, i), i), gto2, center, m+1))
                end
                if gto2.power[i] != 0
                    integral += 0.5 * gto2.power[i] / ζ * (nuclearAttractionIntegral(gto_m1(gto1, i), gto_m1(gto2, i), center, m) - nuclearAttractionIntegral(gto_m1(gto1, i), gto_m1(gto2, i), center, m+1))
                end
            else 
                integral += α/ζ* (gto1.center[i] - gto2.center[i]) * nuclearAttractionIntegral(gto1, gto_m1(gto2, i), center, m)
                integral += (center[i] - gto2.center[i] - α/ζ* (gto1.center[i] - gto2.center[i])) * nuclearAttractionIntegral(gto1, gto_m1(gto2, i), center, m+1)
                if gto2.power[i] > 1
                    integral += 0.5 * gto2.power[i] / ζ * (nuclearAttractionIntegral(gto1, gto_m1(gto_m1(gto2, i), i), center, m) - nuclearAttractionIntegral(gto1, gto_m1(gto_m1(gto2, i), i), center, m+1))
                end
            end
        end
    end

end


function nuclearAttractionIntegral(cgto1::cgto, cgto2::cgto, center::Tuple{Float64,Float64,Float64}, charge::Float64)::Float64
    integral = 0.0
    for i in 1:length(cgto1.gtos)
        for j in 1:length(cgto2.gtos)
            integral += cgto1.gtos[i].contractedCoeff * cgto2.gtos[j].contractedCoeff * nuclearAttractionIntegral(cgto1.gtos[i], cgto2.gtos[j], center, 0)
        end
    end
    return integral*charge
end

function nuclearAttractionIntegral(basis::Array{cgto,1}, )
    
end


function coulombIntegral(gto1::gto, gto2::gto, gto3::gto, gto4::gto, τ::Float64)::Float64


end

function coulombIntegral(cgto1::cgto, cgto2::cgto, cgto3::cgto, cgto4::cgto, τ::Float64)::Float64
    integral = 0.0
    for i in 1:length(cgto1.gtos)
        for j in 1:length(cgto2.gtos)
            for k in 1:length(cgto3.gtos)
                for l in 1:length(cgto4.gtos)
                    integral += cgto1.gtos[i].contractedCoeff * cgto2.gtos[j].contractedCoeff * cgto3.gtos[k].contractedCoeff * cgto4.gtos[l].contractedCoeff * coulombIntegral(cgto1.gtos[i], cgto2.gtos[j], cgto3.gtos[k], cgto4.gtos[l], τ)
                end
            end
        end
    end
    return integral
    
end

function exchangeIntegral(gto1::gto, gto2::gto, gto3::gto, gto4::gto, τ::Float64)::Float64

end

function exchangeIntegral(cgto1::cgto, cgto2::cgto, cgto3::cgto, cgto4::cgto, τ::Float64)::Float64
    integral = 0.0
    for i in 1:length(cgto1.gtos)
        for j in 1:length(cgto2.gtos)
            for k in 1:length(cgto3.gtos)
                for l in 1:length(cgto4.gtos)
                    integral += cgto1.gtos[i].contractedCoeff * cgto2.gtos[j].contractedCoeff * cgto3.gtos[k].contractedCoeff * cgto4.gtos[l].contractedCoeff * exchangeIntegral(cgto1.gtos[i], cgto2.gtos[j], cgto3.gtos[k], cgto4.gtos[l], τ)
                end
            end
        end
    end
    return integral
    
end







overlapIntegral(gto1, gto2) * normalizeConstant(gto1) * normalizeConstant(gto2)
