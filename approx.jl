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

function originalBoundings(P1, P2, cell)
    Q1 = Q2 = []
    minval1 = minval2 = norm(cell[2]-cell[1])
    for i = 1:2
        if abs(P1[1] - cell[i][1]) < minval1
            Q1 = cell[i]
            minval1 = P1[1] - cell[i][1]
        end
        if abs(P1[2] - cell[i][2]) < minval1
            Q1 = cell[i]
            minval1 = P1[2] - cell[i][2]
        end
        if abs(P2[1] - cell[i][1]) < minval2
            Q2 = cell[i]
            minval2 = P2[1] - cell[i][1]
        end
        if abs(P2[2] - cell[i][2]) < minval2
            Q2 = cell[i]
            minval2 = P2[2] - cell[i][2]
        end
    end
    boundings = [ImplicitGeometry.lineThrough(P1, Q1), ImplicitGeometry.lineThrough(P2, Q2)]
    if boundings[1](P2) < 0
        boundings[1] = -boundings[1]
    end
    if boundings[2](P1) < 0
        boundings[2] = -boundings[2]
    end
    boundings
end

function approxIpatch(P1, N1, P2, N2, points, cell, modifiedBoundings = false)
    primaries = [ImplicitGeometry.Line(P1, N1), ImplicitGeometry.Line(P2, N2)]
    origBoundings = originalBoundings(P1, P2, cell)
    if modifiedBoundings
        boundings = [ImplicitGeometry.lineThrough(P1, P1+N1), ImplicitGeometry.lineThrough(P2, P2+N2)]
    else
        boundings = origBoundings
    end
    constraints = collect(zip(primaries, boundings))

    #println(boundings)

    n = 2
    exponent = 2

    insideCell = (p) -> (cell[1] <= p <= cell[2])

    F = (p) -> prod(b->b(p), boundings)^exponent / sum(b->b(p)^2, boundings)
    M(i) = i>n ? F : (p) -> constraints[i][1](p) * prod(b->b(p), setdiff(boundings, [constraints[i][2]]))^exponent / sum(b->b(p)^2, boundings)
    weightsToIPatch(w :: Vector{<:Number}) = IPatch2D.Ipatch(constraints, w)

    if isempty(points)
        ip = ImplicitGeometry.intersection(primaries[1], primaries[2])
        w1 = w2 = 1
        if origBoundings[1](ip) > 0 && origBoundings[2](ip) > 0
            #println("3szog")
            P3 = ip
            a = norm(P2 - P3)
            b = norm(P1 - P3)
            c = norm(P2 - P1)
            #w1 = b
            #w2 = a
			p = (P1 * a + P2 * b + P3 * c) / (a + b + c)
            #println(p)
            eps = (cell[2][1] - cell[1][1]) / 30
            if norm(p-P2) < eps || norm(p-P1) < eps
                p = (P1+P2)/2
            end
        else
            #println("4szog")
            P3 = ImplicitGeometry.intersection(primaries[1], boundings[2])
            P4 = ImplicitGeometry.intersection(boundings[1], primaries[2])
            p = ImplicitGeometry.intersection(ImplicitGeometry.lineThrough(P1,P2), ImplicitGeometry.lineThrough(P3,P4))
        end
		w1 = 1 / M(1)(p)
		w2 = 1 / M(2)(p)
        w0 = -(w1 * M(1)(p) + w2 * M(2)(p)) / F(p)
        w = [w1,w2,w0]
        try
            singularPoint = ImplicitGeometry.intersection(constraints[1][2], constraints[2][2])
            pvalues = [constraints[1][1](singularPoint), constraints[2][1](singularPoint)]
        catch e # if boundaries are parallel, we check in ideal point
            singularPoint = [rotate90(constraints[1][2].n)..., 0]
            pvalues = [ImplicitGeometry.evalHomogenous(ImplicitGeometry.PolynomialCurve(constraints[1][1]), singularPoint), ImplicitGeometry.evalHomogenous(ImplicitGeometry.PolynomialCurve(constraints[2][1]), singularPoint)]
        finally
            if prod(w[1:2] .* pvalues) < 0
                w0 = -(w1*M(1)(p) - w2*M(2)(p)) / F(p)
                w = [w1,-w2,w0]
            end
        end
        #println(w)
        return Ipatch(constraints, w)
    elseif length(points) == 1
        p = points[1]
        w0 = -(M(1)(p) + M(2)(p)) / F(p)
        w = [1,1,w0]
        try
            singularPoint = ImplicitGeometry.intersection(constraints[1][2], constraints[2][2])
            pvalues = [constraints[1][1](singularPoint), constraints[2][1](singularPoint)]
        catch e # if boundaries are parallel, we check in ideal point
            singularPoint = [rotate90(constraints[1][2].n)..., 0]
            pvalues = [ImplicitGeometry.evalHomogenous(ImplicitGeometry.PolynomialCurve(constraints[1][1]), singularPoint), ImplicitGeometry.evalHomogenous(ImplicitGeometry.PolynomialCurve(constraints[2][1]), singularPoint)]
        finally
            if prod(w[1:2] .* pvalues) < 0
                w0 = -(M(1)(p) - M(2)(p)) / F(p)
                w = [1,-1,w0]
            end
        end
        return Ipatch(constraints, w)
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
