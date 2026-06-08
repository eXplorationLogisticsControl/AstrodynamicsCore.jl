"""Lambert's problem algorithm"""
abstract type AbstractLambertOut end

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
Overload method for showing OptimalControlSCPProblem
"""
function Base.show(io::IO, LambertOut::AbstractLambertOut)
    println("Lambert problem solution structure")
    @printf("   Exitflag  : %d\n", LambertOut.exitflag)
    @printf("   TOF       : %1.4e\n", LambertOut.tof)
    println("   Departure : ")
    @printf("       r1 : %1.4e %1.4e %1.4e\n", LambertOut.r1[1], LambertOut.r1[2], LambertOut.r1[3])
    @printf("       v1 : %1.4e %1.4e %1.4e\n", LambertOut.v1[1], LambertOut.v1[2], LambertOut.v1[3])
    println("   Arrival : ")
    @printf("       r2 : %1.4e %1.4e %1.4e\n", LambertOut.r2[1], LambertOut.r2[2], LambertOut.r2[3])
    @printf("       v2 : %1.4e %1.4e %1.4e\n", LambertOut.v2[1], LambertOut.v2[2], LambertOut.v2[3])
end


"""
    lambert(
        r1_vec::Vector,
        r2_vec::Vector,
        tof::Real,
        m::Int,
        μ::Float64,
        cw::Bool = false,
        tol::Float64 = 1.e-12,
        maxiter::Int = 20,
    )

Fast lambert algorithm via Dario Izzo's method

# Args
- `r1_vec::Vector`: initial position vector
- `r2_vec::Vector`: final position vector
- `tof::Real`: time of flight
- `m::Int` number of revolutions
- `μ::Float64`: gravitational parameter
- `cw::Bool`: use clockwise-motion, default is `false`
- `tol::Float64`: tolerance on Newton-Raphson
- `maxiter::Int`: max iteration on Newton-Raphson

# Returns
	`LambertResults`: object containing output of lambert problem
"""
function lambert(
    r1_vec::Vector,
    r2_vec::Vector,
    tof::Real,
    m::Int,
    μ::Float64,
    cw::Bool = false,
    tol::Float64 = 1.e-12,
    maxiter::Int = 20,
)

    # initialize
    exitflag = 0

    # re-scale w.r.t. initial position vector
    r1 = norm(r1_vec)
    vstar = sqrt(μ / r1)
    tstar = r1 / vstar

    # re-scale
    r1_vec = r1_vec / r1
    r2_vec = r2_vec / r1
    tof = tof / tstar

    # geometry parameters
    mr2_vec = norm(r2_vec)
    r1crossr2 = cross(r1_vec, r2_vec)
    dθ = acos_safe(dot(r1_vec, r2_vec) / mr2_vec)
    mcr = norm(r1crossr2)
    nrmunit = r1crossr2 / mcr 

    # --------------------------------------------------------- #
    # decide whether to use the left or right branch (for multi-revolution
    # problems), and the long or short way
    leftbranch = sign(m)
    longway = sign(tof)
    m = abs(m)
    tof = abs(tof)

    # check for direction
    if r1crossr2[3] >= 0
        cw_is_short = true
    else
        cw_is_short = false
    end

    if (cw_is_short == false) && (dθ < π)
        dθ = 2π - dθ
        longway = -1.0
    elseif (cw_is_short == true) && (dθ > π)
        dθ = 2π - dθ
        longway = -1.0
    end
    #println("Corrected dθ: $dθ")

    # chord
    c = sqrt(1 + mr2_vec^2 - 2 * mr2_vec * cos(dθ))

    # semi-perimeter
    s = (1 + mr2_vec + c) / 2   # non-dimensional semi-perimeter
    
    # min.energy semi-major axis
    a_min = s / 2

    # Λ parameter (c.f. BATTIN)
    Λ = sqrt(mr2_vec) * cos(dθ / 2) / s               

    # --------------------------------------------------------- #
    # initialization for Newton-Raphson
    logt = log(tof)

    # single revolution (1 solution)
    if (m == 0)

        # initial values
        inn1 = -0.5233      # first initial guess
        inn2 = +0.3233      # second initial guess
        x1 = log(1 + inn1)# transformed first initial guess
        x2 = log(1 + inn2)# transformed first second guess

    # multiple revolutions (0, 1 or 2 solutions)
    # the returned soltuion depends on the sign of [m]
    else
        # select initial values
        if (leftbranch < 0)
            inn1 = -0.5234 # first initial guess, left branch
            inn2 = -0.2234 # second initial guess, left branch
        else
            inn1 = +0.7234 # first initial guess, right branch
            inn2 = +0.5234 # second initial guess, right branch
        end
        x1 = tan(inn1 * π / 2) # transformed first initial guess
        x2 = tan(inn2 * π / 2) # transformed first second guess
    end

    # since (inn1, inn2) < 0, initial estimate is always ellipse
    xx = [inn1, inn2]
    aa = [a_min / (1 - xx[1]^2), a_min / (1 - xx[2]^2)]
    bβ = [
        longway * 2 * asin(sqrt(((s - c) / 2) / aa[1])),
        longway * 2 * asin(sqrt(((s - c) / 2) / aa[2])),
    ]
    aα = 2 * acos_safe(minimum(xx))
    #println("aα: $aα")
    #println("bβ: $bβ")

    # evaluate the time of flight via Lagrange expression
    y12 = [
        aa[1] * sqrt(aa[1]) * ((aα - sin(aα)) - (bβ[1] - sin(bβ[1])) + 2π * m),
        aa[2] * sqrt(aa[2]) * ((aα - sin(aα)) - (bβ[2] - sin(bβ[2])) + 2π * m),
    ]
    # aa.*sqrt_vector(aa) .*((aα - sin(aα)) .- (bβ-sin_vec(bβ)) .+ 2π*m)
    #println("y12: $y12")

    # initial estimates for y
    if m == 0
        y1 = log(y12[1]) - logt
        y2 = log(y12[2]) - logt
    else
        y1 = y12[1] - tof
        y2 = y12[2] - tof
    end

    # initialize error
    err = Inf
    # storage for x
    x = 0.0

    # --------------------------------------------------------- #
    # Newton-Raphson iterations
    for iter in 1:maxiter
        # update xnew
        xnew = (x1 * y2 - y1 * x2) / (y2 - y1)
        # copy-pasted code (for performance)
        if m == 0
            x = exp(xnew) - 1
        else
            x = atan(xnew) * 2 / π
        end
        a = a_min / (1 - x^2)

        if (x < 1) # ellipse
            β = longway * 2 * asin(sqrt((s - c) / 2 / a))
            α = 2 * acos_safe(x)
        else # hyperbola
            α = 2 * acosh(x)
            β = longway * 2 * asinh(sqrt((s - c) / (-2 * a)))
        end

        # evaluate tof via Lagrange expression
        if (a > 0)
            tof_iter = a * sqrt(a) * ((α - sin(α)) - (β - sin(β)) + 2π * m)
        else
            tof_iter = -a * sqrt(-a) * ((sinh(α) - α) - (sinh(β) - β))
        end

        # update y
        if m == 0
            ynew = log(tof_iter) - logt
        else
            ynew = tof_iter - tof
        end

        # update last two iterations (preventing bouncing)
        x1 = x2
        x2 = xnew
        y1 = y2
        y2 = ynew

        # update error
        err = abs(x1 - xnew)

        # break if error is within tolerance
        if err <= tol
            exitflag = 1
            break
        end
    end

    # if failed, return here
    if exitflag == 0
        return LambertResults(
            μ,
            tof * tstar,
            r1_vec * r1,
            r2_vec * r1,
            [NaN, NaN, NaN],
            [NaN, NaN, NaN],
            exitflag,
        )
    end

    # solution for semi-major axis
    a = a_min / (1 - x^2)

    # calculate ψ
    if (x < 1) # ellipse
        β = longway * 2asin(sqrt((s - c) / 2 / a))
        # ensure acos is well-defined
        α = 2 * acos_safe(x)
        ψ = (α - β) / 2
        η2 = 2a * sin(ψ)^2 / s
        η = sqrt(η2)
    else       # hyperbola
        β = longway * 2asinh(sqrt((c - s) / 2 / a))
        α = 2acosh(x)
        ψ = (α - β) / 2
        η2 = -2a * sinh(ψ)^2 / s
        η = sqrt(η2)
    end

    # normalized normal vector
    ih = longway * nrmunit

    # compute unit-vectors in r1, r2 directions
    r1n = r1_vec / norm(r1_vec)
    r2n = r2_vec / mr2_vec

    # cross-products of ih with unit r1, r2
    ihcrossr1 = cross(ih, r1_vec)
    ihcrossr2 = cross(ih, r2n)

    # radial and tangential components of v1
    vr1 = 1 / η / sqrt(a_min) * (2 * Λ * a_min - Λ - x * η)
    vt1 = sqrt(mr2_vec / a_min / η2 * sin(dθ / 2)^2)

    # radial and tangential components of v2
    vt2 = vt1 / mr2_vec
    vr2 = (vt1 - vt2) / tan(dθ / 2) - vr1

    # obtain velocity vectors
    v1 = (vr1 * r1n + vt1 * ihcrossr1) * vstar
    v2 = (vr2 * r2n + vt2 * ihcrossr2) * vstar

    # construct output
    return LambertResults(μ, tof * tstar, r1_vec * r1, r2_vec * r1, v1, v2, exitflag)
end


"""
    lambert_jac(
        r1_vec::Vector,
        r2_vec::Vector,
        tof::Float64,
        m::Int,
        μ::Float64,
        cw::Bool = false,
        tol::Float64 = 1.e-12,
        maxiter::Int = 20,
    )

Solve Lambert problem and get its sensitivities. 
Fast lambert algorithm via Dario Izzo's method.
Sensitivities from Arora et al, 2015. JGCD. 
*Partial derivatives of the solution to the lambert boundary value problem*.
Adopted from Python implementation by @kdricemt

# Args
	- `r1_vec::Vector`: initial position vector
	- `r2_vec::Vector`: final position vector
	- `tof::Float64`: time of flight
	- `m::Int` number of revolutions
	- `μ::Float64`: gravitational parameter
	- `cw::Bool`: use clockwise-motion, default is `false`
	- `tol::Float64`: tolerance on Newton-Raphson
	- `maxiter::Int`: max iteration on Newton-Raphson

# Returns
	`FastLambertOut`: object containing output of lambert problem
	`Matrix{Float64}: 6x2 matrix of dv/dt
	`Matrix{Float64}: 6x6 matrix of dv/dr
"""
function lambert_jac(
	r1_vec::Vector,
    r2_vec::Vector,
    tof::Float64,
    m::Int,
    μ::Float64,
    cw::Bool = false,
    tol::Float64 = 1.e-12,
    maxiter::Int = 20,
)
	# solve Lambert problem to get velocity vectors
	res = lambert(
	    r1_vec,
	    r2_vec,
	    tof,
	    m,
	    μ,
	    cw,
	    tol,
	    maxiter,
	)
	if res.exitflag == 0
		return res, NaN, NaN
	end

	v1_vec = res.v1
	v2_vec = res.v2

	# define common variables
    r1n = norm(r1_vec)
    r2n = norm(r2_vec)
    v1n = norm(v1_vec)
    v2n = norm(v2_vec)
    k = [0, 0, 1]
    smu = sqrt(μ)

    # true anamaly difference
    cos_dnu = dot(r1_vec, r2_vec)/r1n/r2n
    sin_dnu = sign(dot(cross(r1_vec, r2_vec), k)) * sqrt(1 - cos_dnu^2)

    energy = v1n^2/2 - μ/r1n
    a = -μ/2/energy
    h = cross(r1_vec, v1_vec)
    e = sqrt(1 - norm(h)^2/μ/a)
    p = a * (1 - e^2)

    # "ミッション解析と軌道力学の基礎" p99
    f = 1 - r2n/p * (1 - cos_dnu)
    fdot = sqrt(μ/p) * ((1 - cos_dnu)/sin_dnu) * (1/p * (1 - cos_dnu) - 1/r1n - 1/r2n)
    g = r1n * r2n * sin_dnu/ sqrt(μ * p)
    gdot = 1 - r1n/p * (1 - cos_dnu)

    # kepler STM
    alpha = 1/a
    sig1 = dot(r1_vec, v1_vec) / smu
    sig2 = dot(r2_vec, v2_vec) / smu
    chi = alpha * smu * tof + sig2 - sig1

    # U expressions (eq. 38 - 41 from Arora et al)
    U1 = - r2n * r1n * fdot / smu
    U2 = r1n * (1 - f)
    U3 = smu * (tof - g)
    U4 = U1*U3 - 1/2 * (U2^2 - alpha*U3^2)
    U5 = (chi^3/6 - U3)/alpha
    Cbar = 1/smu * (3 * U5 - chi * U4 - smu * tof * U2)

    # r1 = r1.reshape(3,1)
    # r2 = r2.reshape(3,1)
    # v1 = v1.reshape(3,1)
    # v2 = v2.reshape(3,1)
    dv = v2_vec - v1_vec   # FIXME - check if this is ok for Julia
    dr = r2_vec - r1_vec

    # matrices (eq. 46 - 49 from Arora et al)
    A = r2n/μ * dv * transpose(dv) + 1/r1n^3 * (r1n * (1 - f) * r2_vec*transpose(r1_vec) + Cbar*v2_vec*transpose(r1_vec)) + f*I(3)
    B = r1n/μ * (1- f) * (dr*transpose(v1_vec) - dv*transpose(r1_vec)) + Cbar/μ * v2_vec*transpose(v1_vec) + g*I(3)
    C = -1/r1n^2 * dv*transpose(r1_vec) - 1/r2n^2 * r2_vec*transpose(dv) + fdot*(I(3) - 1/r2n^2 * r2_vec*transpose(r2_vec) + 1/(μ * r2n) * (r2_vec*transpose(v2_vec) - v2_vec*transpose(r2_vec))*r2_vec*transpose(dv)) - μ * Cbar/r2n^3/r1n^3 * r2_vec*transpose(r1_vec)
    D = r1n/μ * dv*transpose(dv) + 1/r2n^3 * (r1n * (1-f) * r2_vec*transpose(r1_vec) - Cbar * r2_vec*transpose(v1_vec)) + gdot * I(3)
    # println("A: $A")
    # println("B: $B")
    # println("C: $C")
    # println("D: $D")

    r1dot = v1_vec
    r2dot = v2_vec
    v1dot = -μ * r1_vec/r1n^3
    v2dot = -μ * r2_vec/r2n^3

    BinvA = inv(B)*A  # FIXME - syntax for Julia?
    CDBinvA = C - D*BinvA
    DBinv = transpose( (transpose(B) \ transpose(D)) )

    # println("BinvA: $BinvA")
    # println("CDBinvA: $CDBinvA")
    # println("DBinv: $DBinv")

    # time

    # from paper -> does not mach result
    # dv1_dt1 = - (v1dot + BinvA @ r1dot)  # 3 x 1
    # dv1_dt2 = CDBinvA @ r2dot
    # dv2_dt1 = CDBinvA @ r1dot
    # dv2_dt2 = v2dot - DBinv @ r2dot

    dv1_dt1 = (v1dot + BinvA*r1dot)  # 3 x 1
    dv1_dt2 = -dv1_dt1
    dv2_dt1 = - CDBinvA*r1dot
    dv2_dt2 = -dv2_dt1

    # construct matrix for time-sensitivities
    dv_dt = [dv1_dt1 dv1_dt2; dv2_dt1 dv2_dt2]  # 6 x 2

    # position and velocity
    dv1_dr1 = - transpose(BinvA)   # 3 x 3
    dv1_dr2 = - transpose(CDBinvA)
    dv2_dr1 = CDBinvA
    dv2_dr2 = DBinv

    # construct matrix for r-sensitivities
    dv_dr = [dv1_dr1 dv1_dr2; dv2_dr1 dv2_dr2]  # 6 x 6

    return res, dv_dt, dv_dr
end