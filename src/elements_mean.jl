"""Keplerian mean / osculating element conversions under J2 (Schaub Appendix G)."""


"""
    _osc_mean_first_order_map_j2(direction, kep; Re, J2)

First-order map between mean and osculating Keplerian elements under J2.

Use: 
- `direction = -1`: convert osculating to mean
- `direction = +1`: convert mean to osculating.

`kep = [a,e,i,Ω,ω,M]` where `M` is the mean anomaly
"""
function _osc_mean_first_order_map_j2(
    direction::Real,
    kep::AbstractVector{<:Real},
    Re::Real,
    J2::Real,
)
    @assert length(kep) == 6 "Keplerian elements must be a length-6 vector"
    a, ecc, inc, Om, w, M = kep
    inc = mod(inc, π)
    Om = _wrap_to_pi(Om)
    w = _wrap_to_pi(w)
    M = _wrap_to_pi(M)

    nu = AstrodynamicsCore.ma2ta(M, ecc)

    g2 = sign(direction) * (J2 / 2) * (Re / a)^2
    η = sqrt(1 - ecc^2)
    g2_p = g2 / η^4
    a_r = (1 + ecc * cos(nu)) / η^2

    de1 = (g2_p / 8) * ecc * η^2 * (
        1 - 11 * cos(inc)^2 - 40 * (cos(inc)^4 / (1 - 5 * cos(inc)^2))
    ) * cos(2w)
    de2 = g2 * ((3 * cos(inc)^2 - 1) / η^6) * (
        ecc * η + ecc / (1 + η) + 3 * cos(nu) + 3 * ecc * cos(nu)^2 + ecc^2 * cos(nu)^3
    )
    de3 = g2 * (3 * (1 - cos(inc)^2) / η^6) * (
        ecc + 3 * cos(nu) + 3 * ecc * cos(nu)^2 + ecc^2 * cos(nu)^3
    ) * cos(2w + 2nu)
    de4 = -g2_p * (1 - cos(inc)^2) * (3 * cos(2w + nu) + cos(2w + 3nu))
    de = de1 + η^2 / 2 * (de2 + de3 + de4)

    di1 = -(ecc * de1) / (η^2 * tan(inc))
    di2 = (g2_p / 2) * cos(inc) * sqrt(1 - cos(inc)^2) * (
        3 * cos(2w + 2nu) + 3 * ecc * cos(2w + nu) + ecc * cos(2w + 3nu)
    )
    di = di1 + di2

    dL0 = M + w + Om
    dL1 = (g2_p / 8) * η^3 * (
        1 - 11 * cos(inc)^2 - 40 * (cos(inc)^4 / (1 - 5 * cos(inc)^2))
    )
    dL2 = -(g2_p / 16) * (
        2 + ecc^2 - 11 * (2 + 3 * ecc^2) * cos(inc)^2 -
        40 * (2 + 5 * ecc^2) * (cos(inc)^4 / (1 - 5 * cos(inc)^2)) -
        400 * ecc^2 * (cos(inc)^6 / (1 - 5 * cos(inc)^2)^2)
    )
    dL3 = (g2_p / 4) * (
        -6 * (1 - 5 * cos(inc)^2) * (nu - M + ecc * sin(nu)) +
        (3 - 5 * cos(inc)^2) * (
            3 * sin(2w + 2nu) + 3 * ecc * sin(2w + nu) + ecc * sin(2w + 3nu)
        )
    )
    dL4 = -(g2_p / 8) * ecc^2 * cos(inc) * (
        11 + 80 * (cos(inc)^2 / (1 - 5 * cos(inc)^2)) +
        200 * (cos(inc)^4 / (1 - 5 * cos(inc)^2)^2)
    )
    dL5 = -(g2_p / 2) * cos(inc) * (
        6 * (nu - M + ecc * sin(nu)) -
        3 * sin(2w + 2nu) - 3 * ecc * sin(2w + nu) - ecc * sin(2w + 3nu)
    )
    dL = mod(dL0 + dL1 + dL2 + dL3 + dL4 + dL5, 2π)

    edM1 = (g2_p / 8) * ecc * η^3 * (
        1 - 11 * cos(inc)^2 - 40 * (cos(inc)^4 / (1 - 5 * cos(inc)^2))
    )
    edM2 = 2 * (3 * cos(inc)^2 - 1) * ((a_r * η)^2 + a_r + 1) * sin(nu)
    edM3 = 3 * (1 - cos(inc)^2) * (
        (-(a_r * η)^2 - a_r + 1) * sin(2w + nu) +
        ((a_r * η)^2 + a_r + 1 / 3) * sin(2w + 3nu)
    )
    edM = edM1 - (g2_p / 4) * η^3 * (edM2 + edM3)

    dOm1 = -(g2_p / 8) * ecc^2 * cos(inc) * (
        11 + 80 * (cos(inc)^2 / (1 - 5 * cos(inc)^2)) +
        200 * (cos(inc)^4 / (1 - 5 * cos(inc)^2)^2)
    )
    dOm2 = -(g2_p / 2) * cos(inc) * (
        6 * (nu - M + ecc * sin(nu)) -
        3 * sin(2w + 2nu) - 3 * ecc * sin(2w + nu) - ecc * sin(2w + 3nu)
    )
    dOm = dOm1 + dOm2

    d1 = (ecc + de) * sin(M) + edM * cos(M)
    d2 = (ecc + de) * cos(M) - edM * sin(M)
    d3 = (sin(inc / 2) + cos(inc / 2) * (di / 2)) * sin(Om) + sin(inc / 2) * dOm * cos(Om)
    d4 = (sin(inc / 2) + cos(inc / 2) * (di / 2)) * cos(Om) - sin(inc / 2) * dOm * sin(Om)

    aO = a + a * g2 * (
        (3 * cos(inc)^2 - 1) * (a_r^3 - 1 / η^3) +
        3 * (1 - cos(inc)^2) * a_r^3 * cos(2w + 2nu)
    )
    eO = sqrt(d1^2 + d2^2)
    iO = 2 * asin(sqrt(d3^2 + d4^2))
    OmO = mod(atan(d3, d4), 2π)
    MO = mod(atan(d1, d2), 2π)
    wO = mod(dL - MO - OmO, 2π)

    OmO = OmO > π ? OmO - 2π : OmO
    wO = wO > π ? wO - 2π : wO
    MO = MO > π ? MO - 2π : MO

    return [aO, eO, iO, OmO, wO, MO]
end


"""
    kep_osc2mean(osc_kep, Re, J2)

Convert osculating Keplerian elements to mean elements under J2 using the
first-order map from Schaub Appendix G.
"""
function kep_osc2mean(osc_kep::AbstractVector{<:Real}, Re::Real, J2::Real)
    return _osc_mean_first_order_map_j2(-1, osc_kep, Re, J2)
end


"""
    kep_mean2osc(mean_kep, Re, J2)

Convert mean Keplerian elements to osculating elements under J2 using the
first-order map from Schaub Appendix G.
"""
function kep_mean2osc(mean_kep::AbstractVector{<:Real}, Re::Real, J2::Real)
    return _osc_mean_first_order_map_j2(1, mean_kep, Re, J2)
end
