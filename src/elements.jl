"""Elements functions"""


"""
Convert 3-element vector from perifocal to geocentric frame
"""
function perifocal2geocentric(vec_pf::Vector, ω::Real, i::Real, Ω::Real)
    # rotate by ω
    v1 = _rotmat_ax3(-ω) * vec_pf
    # rotate by inclination
    v2 = _rotmat_ax1(-i) * v1
    # rotate by Ω
    v_gec = _rotmat_ax3(-Ω) * v2
    return v_gec
end


function kep2rv(kep::Array{<:Real,1}, μ::Real)
	# unpack
	a, e, i, Ω, ω, θ = kep
    # angular momentum
    if e <= 1.0
        h = sqrt(a * μ * (1 - e^2))
    else
        h = sqrt(abs(a) * μ * (e^2 - 1))   # hyperbolic case
    end
    # perifocal vector
    x = (h^2 / μ) * (1 / (1 + e * cos(θ))) * cos(θ)
    y = (h^2 / μ) * (1 / (1 + e * cos(θ))) * sin(θ)
    vx = (μ / h) * (-sin(θ))
    vy = (μ / h) * (e + cos(θ))
    rpf = [x, y, 0.0]
    vpf = [vx, vy, 0.0]

    r_rot3 = perifocal2geocentric(rpf, ω, i, Ω)
    v_rot3 = perifocal2geocentric(vpf, ω, i, Ω)

    # save inertial state
    state_inr = vcat(r_rot3, v_rot3)[:]
    return state_inr
end


function rv2kep(rv::Array{<:Real,1}, μ::Real)
    rvec = rv[1:3]
    vvec = rv[4:6]
    hvec = cross(rvec, vvec)
    r, v, h = norm(rvec), norm(vvec), norm(hvec)
    evec = (1/μ) * cross(vvec, hvec) - rvec / r
    e = norm(evec)
    i = acos(hvec[3] / h)
    if 0 < i < π
        ndir = cross([0,0,1], hvec/h)
    else
        ndir = [1,0,0]
    end
    Ω = atan(ndir[2], ndir[1])
    a = (r * μ) / (2*μ - r*v^2)
    θ = atan(h*dot(rvec, vvec)/r, h^2/r - μ)
    px = dot(rvec, ndir)
    py = dot(rvec, cross(hvec, ndir))/h
    ω = atan(py, px) - θ
    return [a,e,i,Ω,ω,θ]
end


function kep2mee(kep::Array{<:Real,1})
    # unpack
    a, e, i, Ω, ω, θ = kep
    # compute MEEs
    p = a * (1 - e^2)
    f = e*cos(Ω+ω)
    g = e*sin(Ω+ω)
    h = tan(i/2)*cos(Ω)
    k = tan(i/2)*sin(Ω)
    l = Ω + ω + θ
    return [p,f,g,h,k,l]
end


function mee2kep(mee::Array{<:Real,1})
    # unpack
    p, f, g, h, k, l = mee
    # compute Keplerian elements
    e = sqrt(f^2 + g^2)
    a = p / (1 - e^2)
    i = atan(2*sqrt(h^2+k^2), 1-h^2-k^2)
    Ω = atan(k,h)
    ω = atan(g*h-f*k, f*h+g*k)
    θ = l - Ω - ω
    return [a,e,i,Ω,ω,θ]
end


function mee2rv(mee::Array{<:Real,1}, μ::Real)
    p, f, g, h, k, l = mee
    sinL, cosL = sin(l), cos(l)
    α2 = h^2 - k^2
    s2 = 1 + h^2 + k^2
    w = 1 + f * cosL + g * sinL
    r = p / w
    sqrt_mup = sqrt(μ / p)

    r = [r/s2 * (cosL + α2 * cosL + 2 * h * k * sinL);
         r/s2 * (sinL - α2 * sinL + 2 * h * k * cosL);
         2r/s2 * (h * sinL - k * cosL)]
    v = [-1/s2 * sqrt_mup * ( sinL + α2 * sinL - 2 * h * k * cosL + g - 2*f*h*k + α2 * g);
         -1/s2 * sqrt_mup * (-cosL + α2 * cosL + 2 * h * k * sinL - f + 2*g*h*k + α2 * f);
         2/s2 * sqrt_mup * (h * cosL + k * sinL + f * h + g * k)]
    return [r; v]
end


function rv2mee(rv::Array{<:Real,1}, μ::Real, retrograde::Bool = false)
    rvec = rv[1:3]
    vvec = rv[4:6]
    hvec = cross(rvec, vvec)
    r, v, h = norm(rvec), norm(vvec), norm(hvec)
    w = hvec / h
    evec = (1/μ) * cross(vvec, hvec) - rvec / r

    # compute h, k
    h,k = hvec[1:2]

    # compute equinoctial frame
    is_retro = retrograde ? -1 : 1
    k =  w[1] / (1 + is_retro * w[3]);
    h = -w[2] / (1 + is_retro * w[3]);
    den = k * k + h * h + 1;
    fv = [(1. - k * k + h * h) / den;
          (2. * k * h) / den;
          (-2. * is_retro * k) / den]

    gv = [(2. * is_retro * k * h) / den;
          (1. + k * k - h * h) * is_retro / den;
          (2. * h) / den]

    # compute f, g
    g = dot(evec, gv)
    f = dot(evec, fv)

    # compute L
    det1 = (gv[2]*fv[1]-fv[2]*gv[1])    # xy
    det2 = (gv[3]*fv[1]-fv[3]*gv[1])    # xz
    det3 = (gv[3]*fv[2]-fv[3]*gv[2])    # yz
    max_det = max(abs(det1), abs(det2), abs(det3));

    if abs(det1) == max_det
        X = ( gv[2]*rvec[1] - gv[1] * rvec[2]) / det1
        Y = (-fv[2]*rvec[1] + fv[1] * rvec[2]) / det1
    elseif abs(det2) == max_det
        X = ( gv[3]*rvec[1] - gv[1] * rvec[3]) / det2
        Y = (-fv[3]*rvec[1] + fv[1] * rvec[3]) / det2
    else
        X = ( gv[3]*rvec[2] - gv[2] * rvec[3]) / det3
        Y = (-fv[3]*rvec[2] + fv[2] * rvec[3]) / det3
    end
    L = atan(Y/r, X/r)

    # compute p
    p = r * (1 + f * cos(L) + g * sin(L))
    return [p,f,g,h,k,L]
end


"""Convert mean anomaly to eccentric anomaly"""
function ma2ea(M::Real, e::Real; maxiter::Int=20, tol::Float64=1e-14)
    # Initialize starting values
    M = mod(M, 2.0*pi)
    if e < 0.8
        E = M
    else
        E = pi
    end

    # Initialize working variable
    f = E - e*sin(E) - M
    i = 0

    # Iterate until convergence
    for it in 1:maxiter
        if abs(f) < tol
            break
        end
        f = E - e*sin(E) - M
        E = E - f / (1.0 - e*cos(E))
    end
    return E
end


"""Convert mean anomaly to true anomaly"""
function ma2ta(M::Real, e::Real)
	EA = ma2ea(M, e)  # in radians
	ta = 2*atan(sqrt((1+e)/(1-e))*tan(EA/2))
    return ta
end


"""Convert true anomaly to eccentric anomaly"""
function ta2ea(ta::Real, e::Real)
	E = 2*atan(tan(ta/2) / sqrt((1+e)/(1-e)))
    return E
end


"""Convert true anomaly to mean anomaly"""
function ta2ma(ta::Real, e::Real)
	E = ta2ea(ta, e)
    return E - e*sin(E)
end