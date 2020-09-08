# 2020 Ágoston Sipos

module IpatchApprox

using LinearAlgebra

import ..ImplicitGeometry

function Ipatch(constraints :: Array{<:Tuple{<:ImplicitGeometry.ImplicitCurve, <:ImplicitGeometry.ImplicitCurve}, 1},
                weights = ones(length(constraints)+1), exponent = 2) :: ImplicitGeometry.PolynomialCurve
    toPolynomial(c :: ImplicitGeometry.ImplicitCurve) = typeof(c) == ImplicitGeometry.PolynomialCurve ? c : ImplicitGeometry.PolynomialCurve(c)
    newconstraints = map(c->(toPolynomial(c[1]),toPolynomial(c[2])), constraints)
    n = length(constraints)
    boundaries = Set([b for (_,b) in newconstraints])
    patch = sum((i)->weights[i]*newconstraints[i][1] * prod(setdiff(boundaries, [newconstraints[i][2]]))^exponent, 1:n) + weights[end]*prod(boundaries)^exponent
end


rotate90(p :: Vector) = [-p[2], p[1]]

function approxIpatch(P1, N1, P2, N2, points)
    primaries = [ImplicitGeometry.Line(P1, N1), ImplicitGeometry.Line(P2, N2)]
    boundaries = [ImplicitGeometry.lineThrough(P1, P1+N1), ImplicitGeometry.lineThrough(P2, P2+N2)]
    constraints = collect(zip(primaries, boundaries))

    n = 2
    exponent = 2

    F = (p) -> prod(b->b(p), boundaries)^exponent / sum(b->b(p)^2, boundaries)
    M(i) = i>n ? F : (p) -> constraints[i][1](p) * prod(b->b(p), setdiff(boundaries, [constraints[i][2]]))^exponent / sum(b->b(p)^2, boundaries)
    weightsToIPatch(w :: Vector{<:Number}) = IPatch2D.Ipatch(constraints, w)

    if isempty(points)
        w = repeat([1], n)
        return Ipatch(constraints, w)
        #elseif length(points) == 1
        #return ipatch(constraints, points[1], exponent, repeat([1], n), false)
    else
        w = Array{Float64,1}(undef, n+1)
        A = Array{Float64,2}(undef, n+1, n+1)
        b = Array{Float64,1}(undef, n+1)
        for k=1:n+1
            for j=1:n+1
                A[k,j] = sum((p)->M(j)(p)*M(k)(p), points)
            end
            b[k] = sum((p)->F(p)*M(k)(p), points)
        end
        vals, vecs = LinearAlgebra.eigen(A)
        #println(vecs)
        min = argmin(vals)
        w = vecs[:,1]
        pvalues = [1, 1]
        try
            singularPoint = ImplicitGeometry.intersection(constraints[1][2], constraints[2][2])
            pvalues = [constraints[1][1](singularPoint), constraints[2][1](singularPoint)]
        catch e # if boundaries are parallel, we check in ideal point
            singularPoint = [rotate90(constraints[1][2].n)..., 0]
            pvalues = [ImplicitGeometry.evalHomogenous(ImplicitGeometry.PolynomialCurve(constraints[1][1]), singularPoint), ImplicitGeometry.evalHomogenous(ImplicitGeometry.PolynomialCurve(constraints[2][1]), singularPoint)]
        finally
            #println(pvalues)
            pvalues2 = [constraints[1][1]((P1+P2)/2), constraints[2][1]((P1+P2)/2)]
            if prod(w[1:2] .* pvalues) < 0# || (w[1]*pvalues2[1] > 0 && w[2]*pvalues2[2] > 0 && w[3] > 0) || (w[1]*pvalues2[1] < 0 && w[2]*pvalues2[2] < 0 && w[3] < 0)
                w = vecs[:,2]
                #println(vecs[:,1])
                #println(vecs[:,2])
                while prod(w[1:2] .* pvalues) < 0# || (w[1]*pvalues2[1] > 0 && w[2]*pvalues2[2] > 0 && w[3] > 0) || (w[1]*pvalues2[1] < 0 && w[2]*pvalues2[2] < 0 && w[3] < 0)
                    x = 2*rand()-1
                    y = 2*rand()-1
                    w = (x*vecs[:,1] + y*vecs[:,2])
                    w /= norm(w)
                end
            end
        end
        #println(w)# .* vcat(pvalues, 1))
        #w = A\b
        #println(vecs[:,1])
        #println(vecs[:,2])
        #println(cross(vecs[:,1], vecs[:,2]))
        #println(w)
        #println()
    end


    Ipatch([(ImplicitGeometry.Line(P1, N1), ImplicitGeometry.lineThrough(P1, P1+N1)),
            (ImplicitGeometry.Line(P2, N2), ImplicitGeometry.lineThrough(P2, P2+N2))],
           w)
end

function coeffsToVector(c :: ImplicitGeometry.PolynomialCurve)
    result = []
    degree = 4
    for i in 0:degree, j in 0:degree-i
        push!(result, c.coeffs[i+1,j+1])
    end
    result
end

end
