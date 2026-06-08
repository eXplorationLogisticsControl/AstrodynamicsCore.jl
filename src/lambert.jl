"""Lambert's problem algorithm

Adopted from Izzo's algorithm in pykep
"""

abstract type AbstractLambertOut end

"""
    LambertResults

Single-branch Lambert solution returned when `m == 0`.

Fields: `μ`, `tof`, `r1`, `r2`, `v1`, `v2`, `exitflag` (`1` on success).
"""
struct LambertResults <: AbstractLambertOut
    μ::Float64
    tof::Real
    r1::Vector
    r2::Vector
    v1::Vector
    v2::Vector
    exitflag::Int
end

"""
    LambertMultiResults

Multi-branch Lambert solution returned when `m > 0`.

Contains up to `2 * Nmax + 1` velocity pairs (`v1`, `v2` as vectors of vectors),
together with the solved `x` variable, revolution count `revs`, and branch label
`branch` (`:zero`, `:left`, or `:right`) per solution.
"""
struct LambertMultiResults <: AbstractLambertOut
    μ::Float64
    tof::Real
    r1::Vector
    r2::Vector
    v1::Vector{Vector{Float64}}
    v2::Vector{Vector{Float64}}
    x::Vector{Float64}
    revs::Vector{Int}
    branch::Vector{Symbol}
    exitflag::Int
end

function _show_lambert_single(io::IO, r1, v1, r2, v2)
    @printf(io, "       r1 : %1.4e %1.4e %1.4e\n", r1[1], r1[2], r1[3])
    @printf(io, "       v1 : %1.4e %1.4e %1.4e\n", v1[1], v1[2], v1[3])
    println(io, "   Arrival : ")
    @printf(io, "       r2 : %1.4e %1.4e %1.4e\n", r2[1], r2[2], r2[3])
    @printf(io, "       v2 : %1.4e %1.4e %1.4e\n", v2[1], v2[2], v2[3])
end

function Base.show(io::IO, out::LambertResults)
    println(io, "Lambert problem solution structure")
    @printf(io, "   Exitflag  : %d\n", out.exitflag)
    @printf(io, "   TOF       : %1.4e\n", out.tof)
    println(io, "   Departure : ")
    _show_lambert_single(io, out.r1, out.v1, out.r2, out.v2)
end

function Base.show(io::IO, out::LambertMultiResults)
    println(io, "Lambert multi-revolution solution structure")
    @printf(io, "   Exitflag    : %d\n", out.exitflag)
    @printf(io, "   TOF         : %1.4e\n", out.tof)
    @printf(io, "   Solutions   : %d\n", length(out.v1))
    @printf(io, "   Revolutions : %d\n", out.revs[1])
    for k in eachindex(out.v1)
        @printf(io, "   [%d] rev = %d, branch = %s, x = %1.6e\n", k, out.revs[k], out.branch[k], out.x[k])
        @printf(io, "       v1 : %1.4e %1.4e %1.4e\n", out.v1[k][1], out.v1[k][2], out.v1[k][3])
        @printf(io, "       v2 : %1.4e %1.4e %1.4e\n", out.v2[k][1], out.v2[k][2], out.v2[k][3])
    end
end


function _hypergeometric_f(z::Float64, tol::Float64)
    Sj = 1.0
    Cj = 1.0
    err = 1.0
    j = 0
    while err > tol
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1)
        Sj1 = Sj + Cj1
        err = abs(Cj1)
        Sj = Sj1
        Cj = Cj1
        j += 1
    end
    return Sj
end

function _x2tof2(lambda::Float64, x::Float64, N::Int)
    a = 1.0 / (1.0 - x * x)
    if a > 0
        alfa = 2.0 * acos(x)
        beta = 2.0 * asin(sqrt(lambda * lambda / a))
        if lambda < 0.0
            beta = -beta
        end
        return (a * sqrt(a) * ((alfa - sin(alfa)) - (beta - sin(beta)) + 2.0 * π * N)) / 2.0
    else
        alfa = 2.0 * acosh(x)
        beta = 2.0 * asinh(sqrt(-lambda * lambda / a))
        if lambda < 0.0
            beta = -beta
        end
        return -a * sqrt(-a) * ((beta - sinh(beta)) - (alfa - sinh(alfa))) / 2.0
    end
end

function _x2tof(lambda::Float64, x::Float64, N::Int)
    battin = 0.01
    lagrange = 0.2
    dist = abs(x - 1)
    if dist < lagrange && dist > battin
        return _x2tof2(lambda, x, N)
    end
    K = lambda * lambda
    E = x * x - 1.0
    ρ = abs(E)
    z = sqrt(1 + K * E)
    if dist < battin
        η = z - lambda * x
        S1 = 0.5 * (1.0 - lambda - x * η)
        Q = _hypergeometric_f(S1, 1e-11)
        Q = 4.0 / 3.0 * Q
        return (η^3 * Q + 4.0 * lambda * η) / 2.0 + N * π / ρ^1.5
    else
        y = sqrt(ρ)
        g = x * z - lambda * E
        if E < 0
            l = acos(g)
            d = N * π + l
        else
            f = y * (z - lambda * x)
            d = log(f + g)
        end
        return (x - lambda * z - d / y) / E
    end
end

function _dtdx(lambda::Float64, x::Float64, T::Float64)
    l2 = lambda * lambda
    l3 = l2 * lambda
    umx2 = 1.0 - x * x
    y = sqrt(1.0 - l2 * umx2)
    y2 = y * y
    y3 = y2 * y
    DT = 1.0 / umx2 * (3.0 * T * x - 2.0 + 2.0 * l3 * x / y)
    DDT = 1.0 / umx2 * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - l2) * l3 / y3)
    DDDT = 1.0 / umx2 * (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - l2) * l2 * l3 * x / y3 / y2)
    return DT, DDT, DDDT
end

function _householder(
    lambda::Float64,
    T::Float64,
    x0::Float64,
    N::Int,
    eps::Float64,
    iter_max::Int,
)
    it = 0
    err = 1.0
    x = x0
    while err > eps && it < iter_max
        tof = _x2tof(lambda, x, N)
        DT, DDT, DDDT = _dtdx(lambda, x, tof)
        delta = tof - T
        DT2 = DT * DT
        xnew = x - delta * (DT2 - delta * DDT / 2.0) / (DT * (DT2 - delta * DDT) + DDDT * delta^2 / 6.0)
        err = abs(x - xnew)
        x = xnew
        it += 1
    end
    return x, it, err
end

function _lambert_geometry(r0::Vector, r1::Vector, cw::Bool)
    c = norm(r1 - r0)
    Rs = norm(r0)
    Rf = norm(r1)
    s = (c + Rs + Rf) / 2.0

    irs = r0 / Rs
    irf = r1 / Rf
    ih = cross(irs, irf)
    ih_norm = norm(ih)
    if ih_norm == 0.0 || ih[3] == 0.0
        throw(ArgumentError(
            "lambert: angular momentum has no z component; cannot define transfer direction automatically",
        ))
    end
    ih = ih / ih_norm

    lambda2 = 1.0 - c / s
    lambda = sqrt(lambda2)

    its = cross(ih, irs)
    itf = cross(ih, irf)
    its = its / norm(its)
    itf = itf / norm(itf)

    if ih[3] < 0.0
        lambda = -lambda
        its = -its
        itf = -itf
    end
    if cw
        lambda = -lambda
        its = -its
        itf = -itf
    end

    return (;
        c,
        s,
        lambda,
        lambda2,
        Rs,
        Rf,
        irs,
        irf,
        its,
        itf,
    )
end

function _lambert_nmax(lambda::Float64, lambda2::Float64, T::Float64, multi_revs::Int)
    Nmax = min(multi_revs, Int(floor(T / π)))
    T00 = acos(lambda) + lambda * sqrt(1.0 - lambda2)
    T0 = T00 + Nmax * π
    if Nmax > 0 && T < T0
        T_min = T0
        x_old = 0.0
        for it in 0:12
            DT, DDT, DDDT = _dtdx(lambda, x_old, T_min)
            x_new = if DT != 0.0
                x_old - DT * DDT / (DDT * DDT - DT * DDDT / 2.0)
            else
                x_old
            end
            err = abs(x_old - x_new)
            if err < 1e-13
                break
            end
            T_min = _x2tof(lambda, x_new, Nmax)
            x_old = x_new
        end
        if T_min > T
            Nmax -= 1
        end
    end
    return min(multi_revs, Nmax)
end

function _lambert_initial_x0(lambda::Float64, lambda2::Float64, T::Float64)
    lambda3 = lambda * lambda2
    T00 = acos(lambda) + lambda * sqrt(1.0 - lambda2)
    T1 = 2.0 / 3.0 * (1.0 - lambda3)
    if T >= T00
        return -(T - T00) / (T - T00 + 4)
    elseif T <= T1
        return T1 * (T1 - T) / (2.0 / 5.0 * (1 - lambda2 * lambda3) * T) + 1
    else
        return (T / T00)^(0.69314718055994529 / log(T1 / T00)) - 1.0
    end
end

function _lambert_velocities_from_x(
    x::Float64,
    geom,
    μ::Float64,
)
    gamma = sqrt(μ * geom.s / 2.0)
    ρ = (geom.Rs - geom.Rf) / geom.c
    σ = sqrt(1 - ρ^2)
    y = sqrt(1.0 - geom.lambda2 + geom.lambda2 * x^2)
    vrs = gamma * ((geom.lambda * y - x) - ρ * (geom.lambda * y + x)) / geom.Rs
    vrf = -gamma * ((geom.lambda * y - x) + ρ * (geom.lambda * y + x)) / geom.Rf
    vt = gamma * σ * (y + geom.lambda * x)
    vts = vt / geom.Rs
    vtf = vt / geom.Rf
    v0 = vrs * geom.irs + vts * geom.its
    v1 = vrf * geom.irf + vtf * geom.itf
    return v0, v1
end

function _solve_lambert_branches(
    geom,
    μ::Float64,
    tof::Real,
    multi_revs::Int;
    maxiter::Int = 15,
    tol::Float64 = 1e-8,
)
    T = sqrt(2.0 * μ / geom.s^3) * tof
    Nmax = _lambert_nmax(geom.lambda, geom.lambda2, T, multi_revs)

    nsol = Nmax * 2 + 1
    xs = Vector{Float64}(undef, nsol)
    iters = Vector{Int}(undef, nsol)
    errs = Vector{Float64}(undef, nsol)
    revs = Vector{Int}(undef, nsol)
    branches = Vector{Symbol}(undef, nsol)

    x0_init = _lambert_initial_x0(geom.lambda, geom.lambda2, T)
    xs[1], iters[1], errs[1] = _householder(geom.lambda, T, x0_init, 0, tol, maxiter)
    revs[1] = 0
    branches[1] = :zero

    for i in 1:Nmax
        tmp = ((i * π + π) / (8.0 * T))^(2.0 / 3.0)
        xs[2i] = (tmp - 1) / (tmp + 1)
        xs[2i], iters[2i], errs[2i] = _householder(geom.lambda, T, xs[2i], i, tol, maxiter)
        revs[2i] = i
        branches[2i] = :left

        tmp = ((8.0 * T) / (i * π))^(2.0 / 3.0)
        xs[2i + 1] = (tmp - 1) / (tmp + 1)
        xs[2i + 1], iters[2i + 1], errs[2i + 1] =
            _householder(geom.lambda, T, xs[2i + 1], i, tol, maxiter)
        revs[2i + 1] = i
        branches[2i + 1] = :right
    end

    v1_list = Vector{Vector{Float64}}(undef, nsol)
    v2_list = Vector{Vector{Float64}}(undef, nsol)
    for k in 1:nsol
        eps_k = k == 1 ? tol : tol
        if errs[k] <= eps_k
            v1_list[k], v2_list[k] = _lambert_velocities_from_x(xs[k], geom, μ)
        else
            v1_list[k] = fill(NaN, 3)
            v2_list[k] = fill(NaN, 3)
        end
    end

    exitflag = errs[1] <= tol ? 1 : 0
    return Nmax, xs, v1_list, v2_list, revs, branches, exitflag
end


"""
    lambert(
        r1_vec::Vector,
        r2_vec::Vector,
        tof::Real,
        m::Int,
        μ::Float64,
        cw::Bool = false;
        tol::Float64 = 1.e-12,
        maxiter::Int = 20,
    )

Fast Lambert algorithm via Dario Izzo's method (kep3 port).

# Args
- `r1_vec::Vector`: initial position vector
- `r2_vec::Vector`: final position vector
- `tof::Real`: time of flight
- `m::Int`: maximum number of full revolutions (kep3 `multi_revs`)
- `μ::Float64`: gravitational parameter
- `cw::Bool`: retrograde (clockwise) transfer when `true`
- `tol::Float64`: unused legacy keyword (kep3 tolerances are fixed internally)
- `maxiter::Int`: maximum Householder iterations per branch

# Returns
- `m == 0`: `LambertResults` with a single `(v1, v2)` pair
- `m > 0`: `LambertMultiResults` with up to `2*Nmax + 1` solution branches
  (0-rev plus left/right pairs for each revolution `1..Nmax`)
"""
function lambert(
    r1_vec::Vector,
    r2_vec::Vector,
    tof::Real,
    m::Int,
    μ::Float64,
    cw::Bool = false;
    tol::Float64 = 1.e-12,
    maxiter::Int = 20,
)
    if tof <= 0
        throw(ArgumentError("lambert: time of flight must be positive"))
    end
    if μ <= 0
        throw(ArgumentError("lambert: gravitational parameter must be positive"))
    end
    if m < 0
        throw(ArgumentError("lambert: maximum revolutions m must be non-negative"))
    end

    geom = _lambert_geometry(r1_vec, r2_vec, cw)
    Nmax, xs, v1_list, v2_list, revs, branches, exitflag =_solve_lambert_branches(geom, μ, tof, m; maxiter=maxiter, tol=tol)

    if m == 0
        if exitflag == 0
            return LambertResults(
                μ,
                tof,
                copy(r1_vec),
                copy(r2_vec),
                fill(NaN, 3),
                fill(NaN, 3),
                exitflag,
            )
        end
        return LambertResults(μ, tof, copy(r1_vec), copy(r2_vec), v1_list[1], v2_list[1], exitflag)
    end

    return LambertMultiResults(
        μ,
        tof,
        copy(r1_vec),
        copy(r2_vec),
        v1_list,
        v2_list,
        xs,
        revs,
        branches,
        exitflag,
    )
end


"""
    lambert_jac(
        r1_vec::Vector,
        r2_vec::Vector,
        tof::Float64,
        m::Int,
        μ::Float64,
        cw::Bool = false;
        tol::Float64 = 1.e-12,
        maxiter::Int = 20,
    )

Solve Lambert problem and get its sensitivities.
Fast lambert algorithm via Dario Izzo's method.
Sensitivities from Arora et al, 2015. JGCD.
*Partial derivatives of the solution to the lambert boundary value problem*.
Adopted from Python implementation by @kdricemt

Jacobian support is limited to the single 0-revolution solution (`m == 0`).

# Returns
- `LambertResults`: Lambert solution
- `Matrix{Float64}`: 6×2 matrix of dv/dt
- `Matrix{Float64}`: 6×6 matrix of dv/dr
"""
function lambert_jac(
    r1_vec::Vector,
    r2_vec::Vector,
    tof::Float64,
    m::Int,
    μ::Float64,
    cw::Bool = false;
    tol::Float64 = 1.e-12,
    maxiter::Int = 20,
)
    m == 0 || throw(ArgumentError("lambert_jac: only m == 0 is supported"))

    res = lambert(r1_vec, r2_vec, tof, m, μ, cw; tol=tol, maxiter=maxiter)
    if res.exitflag == 0
        return res, NaN, NaN
    end

    v1_vec = res.v1
    v2_vec = res.v2

    r1n = norm(r1_vec)
    r2n = norm(r2_vec)
    v1n = norm(v1_vec)
    v2n = norm(v2_vec)
    k = [0, 0, 1]
    smu = sqrt(μ)

    cos_dnu = dot(r1_vec, r2_vec) / r1n / r2n
    sin_dnu = sign(dot(cross(r1_vec, r2_vec), k)) * sqrt(1 - cos_dnu^2)

    energy = v1n^2 / 2 - μ / r1n
    a = -μ / 2 / energy
    h = cross(r1_vec, v1_vec)
    e = sqrt(1 - norm(h)^2 / μ / a)
    p = a * (1 - e^2)

    f = 1 - r2n / p * (1 - cos_dnu)
    fdot = sqrt(μ / p) * ((1 - cos_dnu) / sin_dnu) * (1 / p * (1 - cos_dnu) - 1 / r1n - 1 / r2n)
    g = r1n * r2n * sin_dnu / sqrt(μ * p)
    gdot = 1 - r1n / p * (1 - cos_dnu)

    alpha = 1 / a
    sig1 = dot(r1_vec, v1_vec) / smu
    sig2 = dot(r2_vec, v2_vec) / smu
    chi = alpha * smu * tof + sig2 - sig1

    U1 = -r2n * r1n * fdot / smu
    U2 = r1n * (1 - f)
    U3 = smu * (tof - g)
    U4 = U1 * U3 - 1 / 2 * (U2^2 - alpha * U3^2)
    U5 = (chi^3 / 6 - U3) / alpha
    Cbar = 1 / smu * (3 * U5 - chi * U4 - smu * tof * U2)

    dv = v2_vec - v1_vec
    dr = r2_vec - r1_vec

    A = r2n / μ * dv * transpose(dv) +
        1 / r1n^3 * (r1n * (1 - f) * r2_vec * transpose(r1_vec) + Cbar * v2_vec * transpose(r1_vec)) +
        f * I(3)
    B = r1n / μ * (1 - f) * (dr * transpose(v1_vec) - dv * transpose(r1_vec)) +
        Cbar / μ * v2_vec * transpose(v1_vec) + g * I(3)
    C = -1 / r1n^2 * dv * transpose(r1_vec) - 1 / r2n^2 * r2_vec * transpose(dv) +
        fdot * (I(3) - 1 / r2n^2 * r2_vec * transpose(r2_vec) +
                1 / (μ * r2n) * (r2_vec * transpose(v2_vec) - v2_vec * transpose(r2_vec)) * r2_vec * transpose(dv)) -
        μ * Cbar / r2n^3 / r1n^3 * r2_vec * transpose(r1_vec)
    D = r1n / μ * dv * transpose(dv) +
        1 / r2n^3 * (r1n * (1 - f) * r2_vec * transpose(r1_vec) - Cbar * r2_vec * transpose(v1_vec)) +
        gdot * I(3)

    r1dot = v1_vec
    r2dot = v2_vec
    v1dot = -μ * r1_vec / r1n^3
    v2dot = -μ * r2_vec / r2n^3

    BinvA = inv(B) * A
    CDBinvA = C - D * BinvA
    DBinv = transpose((transpose(B) \ transpose(D)))

    dv1_dt1 = (v1dot + BinvA * r1dot)
    dv1_dt2 = -dv1_dt1
    dv2_dt1 = -CDBinvA * r1dot
    dv2_dt2 = -dv2_dt1

    dv_dt = [dv1_dt1 dv1_dt2; dv2_dt1 dv2_dt2]

    dv1_dr1 = -transpose(BinvA)
    dv1_dr2 = -transpose(CDBinvA)
    dv2_dr1 = CDBinvA
    dv2_dr2 = DBinv

    dv_dr = [dv1_dr1 dv1_dr2; dv2_dr1 dv2_dr2]

    return res, dv_dt, dv_dr
end
