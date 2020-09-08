"""
Copyright 2019 Ágoston Sipos

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

module ImplicitGeometry

using LinearAlgebra, PaddedViews, Polynomials
import Base.+, Base.-, Base.*, Base.^, LinearAlgebra.adjoint
export ImplicitCurve, Line, PolynomialCurve, adjoint

"""
Abstract type for storing an implicit plane curve

Can have multiple representations but any of them should be convertible
to `PolynomialCurve` (see below example with `Line`)
"""
abstract type ImplicitCurve <: Function end

"""
    Line(p, n)

Is an implicit representation for the line
interpolating `p` and having normal vector `n`.
"""
struct Line <: ImplicitCurve
    p :: Vector{Float64}
    n :: Vector{Float64}
end

"""
Evaluating the line in a 2D point
"""
function (l :: Line)(q :: Vector) :: Float64
    return dot(q-l.p, l.n)
end

"""
Flipping a line
"""
(-)(l :: Line) = Line(l.p, -l.n)

"""
Line interpolating two given points
"""
lineThrough(p1 :: Vector, p2 :: Vector) = Line(p1, normalize([p1[2]-p2[2], p2[1]-p1[1]]))

"""
Point where two lines intersect
"""
function intersection(l1 :: Line, l2 :: Line) :: Vector
    A = [l1.n'; l2.n']
    b = [dot(l1.n, l1.p); dot(l2.n, l2.p)]
    return A\b
end

"""
    PolynomialCurve

An object representing a two-variable polynomial, i.e. an algebraic curve in plane
`coeffs` is an (N+1)×(N+1) matrix, where N is the degree. The bottom-right triangle is always zero.
"""
struct PolynomialCurve <: ImplicitCurve
    coeffs :: Array{Float64, 2}
end

"""
Evaluating the polynomial in a point
"""
function (c :: PolynomialCurve)(q :: Vector)
    N = size(c.coeffs, 1)-1
    return sum(c.coeffs .* hcat(repeat([q[1].^(0:N)],N+1)...) .* vcat(repeat([q[2].^(0:N)]',N+1)...))
end

"""
Evaluating the polynomial in a point in homogenous coordinates (R^3)
"""
function evalHomogenous(c :: PolynomialCurve, q :: Vector)
    N = size(c.coeffs, 1)-1
    homComp = ones(N+1,N+1)
    for i=1:N
        for j=1:i
            homComp[i,j] = q[3]^(N+2-i-j)
        end
    end
    return sum(c.coeffs .* hcat(repeat([q[1].^(0:N)],N+1)...) .* vcat(repeat([q[2].^(0:N)]',N+1)...) .* homComp)
end

"""
Do-nothing constructor, so that a constructor recall does not cause an error
"""
PolynomialCurve(c :: PolynomialCurve) = PolynomialCurve(c.coeffs)

"""
Converting a line to polynomial form
"""
function PolynomialCurve(l :: Line)
    coeffs = [-dot(l.n, l.p) l.n[2]; l.n[1] 0]
    PolynomialCurve(coeffs)
end

"""
Adding two polynomials together (degree will be max of degrees)
"""
function (+)(c1 :: PolynomialCurve, c2 :: PolynomialCurve) :: PolynomialCurve
    N = max(size(c1.coeffs, 1), size(c2.coeffs, 1))
    PolynomialCurve(PaddedView(0, c1.coeffs, (N,N)) .+ PaddedView(0, c2.coeffs, (N,N)))
end

"""
Multiplying two polynomials (degree will be sum of degrees)
"""
function (*)(c1 :: ImplicitCurve, c2 :: ImplicitCurve) :: PolynomialCurve
    c1 = typeof(c1) == PolynomialCurve ? c1 : PolynomialCurve(c1)
    c2 = typeof(c2) == PolynomialCurve ? c2 : PolynomialCurve(c2)
    N = size(c1.coeffs, 1) - 1
    M = size(c2.coeffs, 1) - 1
    coeffs = zeros(N+M+1, N+M+1)
    for i = 1:N+M+1
        for j = 1:N+M+2-i
            for k = max(1,i-M):min(i,N+1)
                for l = max(1,j-M):min(j,N+1)
                    coeffs[i,j] += c1.coeffs[k,l] * c2.coeffs[i-k+1,j-l+1]
                end
            end
        end
    end
    PolynomialCurve(coeffs)
end

"""
Multiplying a polynomial with a number
"""
*(a :: Number, c :: ImplicitCurve) = PolynomialCurve(a .* PolynomialCurve(c).coeffs)
*(c :: ImplicitCurve, a :: Number) = PolynomialCurve(a .* PolynomialCurve(c).coeffs)

"""
Raising polynomial to an integer power
"""
^(c :: ImplicitCurve, n :: Integer) = prod(repeat([PolynomialCurve(c)],n))

"""
This is the gradient operator for implicit curves, use it as `c'`
"""
function adjoint(c :: ImplicitCurve) :: Function
    c = PolynomialCurve(c)
    N = size(c.coeffs, 1)
    if N <= 1
        return p->zeros(2)
    end
    grady = PolynomialCurve(c.coeffs[1:N-1, 2:N])
    gradx = PolynomialCurve(c.coeffs[2:N, 1:N-1])
    for i=1:N-1
        grady.coeffs[:,i] *= i
        gradx.coeffs[i,:] *= i
    end
    return p->[gradx(p), grady(p)]
end

"""
This is the hessian operator for implicit curves
"""
function hessian(c :: ImplicitCurve) :: Function
    c = PolynomialCurve(c)
    N = size(c.coeffs, 1)
    if N <= 2
        return p->zeros(2,2)
    end
    grady = PolynomialCurve(c.coeffs[1:N-1, 2:N])
    gradx = PolynomialCurve(c.coeffs[2:N, 1:N-1])
    for i=1:N-1
        grady.coeffs[:,i] *= i
        gradx.coeffs[i,:] *= i
    end
    return p->[(gradx'(p))'; (grady'(p))']
end

function laplace(c :: ImplicitCurve) :: PolynomialCurve
    c = PolynomialCurve(c)
    N = size(c.coeffs, 1)
    if N <= 2
        return PolynomialCurve([0])
    end
    grady = PolynomialCurve(c.coeffs[1:N-1, 2:N])
    gradx = PolynomialCurve(c.coeffs[2:N, 1:N-1])
    for i=1:N-1
        grady.coeffs[:,i] *= i
        gradx.coeffs[i,:] *= i
    end
    N -= 1
    gradyy = PolynomialCurve(grady.coeffs[1:N-1, 2:N])
    gradxx = PolynomialCurve(gradx.coeffs[2:N, 1:N-1])
    for i=1:N-1
        gradyy.coeffs[:,i] *= i
        gradxx.coeffs[i,:] *= i
    end
    return gradxx + gradyy
end

function normal(c :: ImplicitCurve) :: Function
    return p -> normalize(c'(p))
end

function curvature(c :: ImplicitCurve) :: Function
    function (p)
        G = c'(p)
        H = hessian(c)(p)
        return (G[2]^2 * H[1,1] - 2 * G[1] * G[2] * H[1,2] + G[1]^2 * H[2,2]) / (G[1]^2 + G[2]^2)^(3/2)
    end
end

function integralxy(c :: ImplicitCurve) :: PolynomialCurve
    c = PolynomialCurve(c)
    N = size(c.coeffs, 1)
    coeffs = zeros(N+2, N+2)
    coeffs[2:N+1, 2:N+1] = c.coeffs
    for i=2:N+1
        coeffs[:,i] /= (i-1)
        coeffs[i,:] /= (i-1)
    end
    return PolynomialCurve(coeffs)
end

"""
Point where axis-oriented line and conic intersect
"""
function intersection(l :: Line, c :: ImplicitCurve) :: Vector
    N = size(c.coeffs, 1)
    isinside = (x)->0<=x<=1
    rank = (x) -> abs(x-.5)
    realroots(p) = convert(Array{Float64}, filter((z)->imag(z)==0, roots(p)))
    bestrealroot(p) = isinside(sort(realroots(p), by=rank)[1]) ? sort(realroots(p), by=rank)[1] : (abs(polyval(p,0)) < abs(polyval(p,1)) ? 0.0 : 1.0)
    if l.n == [1,0] || l.n == [-1,0] # x = c type line
        x = l.p[1]
        coeffs = zeros(N)
        for i=1:N
            coeffs[i] = sum(j->c.coeffs[j,i]*x^(j-1), 1:N)
        end
        y = bestrealroot(Poly(coeffs))
        return [x,y]
    elseif l.n == [0,1] || l.n == [0,-1] # y = c type line
        y = l.p[2]
        coeffs = zeros(N)
        for i=1:N
            coeffs[i] = sum(j->c.coeffs[i,j]*y^(j-1), 1:N)
        end
        x = bestrealroot(Poly(coeffs))
        return [x,y]
    else
        throw("Non-axis-oriented line and conic intersection not implemented.")
    end
end

end
