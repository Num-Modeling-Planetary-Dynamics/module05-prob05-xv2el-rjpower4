# ========================================================================================
# File: keplerian_elements.jl
# Brief: Structure and related methods for keplerian element sets
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

using LinearAlgebra

"""
    KeplerianElements

Container for standard keplerian elements.
"""
struct KeplerianElements
    sma::Float64
    ecc::Float64
    inc::Float64
    aop::Float64
    raan::Float64
    ta::Float64
end

function KeplerianElements(sma; ecc=0.0, inc=0.0, aop=0.0, raan=0.0, ta=0.0)
    if sma == 0.0
        throw(ArgumentError("zero semi-major axis"))
    elseif ecc < 0.0
        throw(ArgumentError("negative eccentricity"))
    elseif (ecc == 1.0 && sma != inf) || (sign(sma * (1 - ecc)) < 0)
        throw(ArgumentError("incompatible eccentricity and semi-major axis"))
    end
    return KeplerianElements(
        sma,
        ecc,
        inc,
        aop,
        raan,
        ta
    )
end

# ----------------------------------------------------------------------------------------
# Basic Getters
# ----------------------------------------------------------------------------------------
semi_major_axis(ke::KeplerianElements) = ke.sma
eccentricity(ke::KeplerianElements) = ke.ecc
inclination(ke::KeplerianElements) = ke.inc
argument_of_periapsis(ke::KeplerianElements) = ke.aop
right_ascension(ke::KeplerianElements) = ke.raan
true_anomaly(ke::KeplerianElements) = ke.ta

# ----------------------------------------------------------------------------------------
# Standard Two-Body Geometric Relations
# ----------------------------------------------------------------------------------------
semi_minor_axis(ke::KeplerianElements) = semi_major_axis(ke) * sqrt(1 - eccentricity(ke)^2)
semi_latus_rectum(ke::KeplerianElements) = semi_major_axis(ke) * (1 - eccentricity(ke)^2)

periapsis_radius(ke::KeplerianElements) = semi_major_axis(ke) * (1 - eccentricity(ke))
apiapsis_radius(ke::KeplerianElements) = semi_major_axis(ke) * (1 + eccentricity(ke))

function apsis_radii(ke::KeplerianElements)
    a = semi_major_axis(ke)
    ae = a * eccentricity(ke)
    return (a - ae, a + ae)
end

# ----------------------------------------------------------------------------------------
# Standard Two-Body Dynamical Relations
# ----------------------------------------------------------------------------------------
mean_motion(gm, ke::KeplerianElements) = sqrt(gm / semi_major_axis(ke)^3)
mean_motion(gm, a::AbstractFloat) = mean_motion(gm, KeplerianElements(a))

orbital_period(gm, ke::KeplerianElements) = 2pi / mean_motion(gm, ke)
orbital_period(gm, a::AbstractFloat) = orbital_period(gm, KeplerianElements(a))

specific_energy(gm, ke::KeplerianElements) = -gm / (2 * semi_major_axis(ke))
specific_energy(gm, a::AbstractFloat) = specific_energy(gm, KeplerianElements(a))

function specific_angular_momentum(gm, ke::KeplerianElements) 
    return sqrt(gm * semi_major_axis(ke) * (1 - eccentricity(ke)^2))
end

function specific_angular_momentum(gm, sma, ecc) 
    return specific_angular_momentum(gm, KeplerianElements(sma, ecc=ecc))
end

function orbital_radius(gm, ke::KeplerianElements; ta=true_anomaly(ke)) 
    return parameter(ke) / (1 + eccentricity(ke) * cos(ta))
end

function orbital_velocity(gm, ke::KeplerianElements; ta=true_anomaly(ke))
    energy = specific_energy(gm, ke)
    radius = orbital_radius(gm, ke; ta=ta)
    return sqrt(2 * (energy + gm / r))
end

# ----------------------------------------------------------------------------------------
# Other
# ----------------------------------------------------------------------------------------
function inclination(angular_momentum::AbstractVector) 
    return acos(angular_momentum[3] / norm(angular_momentum))
end

specific_energy(gm, r::T, v::T) where {T <: AbstractFloat} = (1//2) * v^2 - gm / r

# ----------------------------------------------------------------------------------------
# Conversion from cartestian vector
# ----------------------------------------------------------------------------------------
function KeplerianElements(gm, state::AbstractVector)
    # Helper functions
    wrap_if(val, condition) = condition ? 2pi - val : val
    angle_between(v1, v2) = acos(dot(v1, v2))
    normalize(v) = v ./ norm(v)

    r = state[1:3]
    v = state[4:6]

    r_mag = norm(r)
    r_hat = r ./ r_mag
    v_mag = norm(v)

    # Semi-Major Axis
    spec_energy = specific_energy(gm , r_mag, v_mag)
    sma = -gm / (2 * spec_energy)

    # Inclination
    h = cross(r, v)
    inc = inclination(h)

    # Eccentricity
    e_vector = (1/gm) * (v_mag^2 * r - dot(r, v) * v) - r_hat
    ecc = norm(e_vector)
    e_hat = e_vector ./ ecc

    
    nodal_axis = (-h[2], h[1], 0.0)
    nodal_axis_hat = h[2] == 0 ? (1.0, 0.0, 0.0) : normalize(nodal_axis) 

    # Right Ascension, Argument of Periapsis, and True anomaly
    raan = wrap_if(acos(nodal_axis_hat[1]), nodal_axis_hat[2] < 0.0)
    aop = wrap_if(angle_between(nodal_axis_hat, e_hat), e_vector[3] < 0.0)
    ta = clamp(dot(e_hat, r_hat), -1.0, 1.0) |> acos

    return KeplerianElements(
        sma, 
        ecc=ecc, 
        inc=inc, 
        aop=aop, 
        raan=raan, 
        ta=ta
    )
end

# ----------------------------------------------------------------------------------------
# Conversion to Cartestian Vector
# ----------------------------------------------------------------------------------------
function perifocal_state(gm, ke::KeplerianElements)
    n = mean_motion(gm, ke)
    
    sma = semi_major_axis(ke)
    ecc = eccentricity(ke)
    ta = true_anomaly(ke)

    r = orbital_radius(gm, ke)

    alpha = (sma * n) / sqrt(1 - ecc^2)

    v_x = - alpha * sin(ta)
    v_y = alpha * (ecc + cos(ta))
    v_z = zero(v_x)

    x = r * cos(ta)
    y = r * sin(ta)
    z = zero(x)

    return [x, y, z, v_x, v_y, vz]
end


"""
    perifocal_to_inertial(::KeplerianElements)

Return rotation matrix converting states in perifocal (e, p, h) to inertial (X, Y, Z) frame.
"""
function perifocal_to_inertial(ke::KeplerianElements)
    s_aop, cos_aop = sincos(right_ascension(ke))
    s_inc, cos_inc = sincos(right_ascension(ke))
    s_raan, cos_raan = sincos(right_ascension(ke))

    return [
        (cos_raan * cos_aop - sin_raan * sin_aop * cos_inc) (-cos_raan * sin_aop - sin_raan * cos_aop * cos_inc) (sin_raan * sin_inc);
        (sin_raan * cos_aop + cos_raan * sin_aop * cos_inc) (-sin_raan * sin_aop + cos_raan * cos_aop * cos_inc) (-cos_raan * sin_inc);
        (sin_aop * sin_inc) (cos_aop * sin_inc) cos_inc 
    ]
end