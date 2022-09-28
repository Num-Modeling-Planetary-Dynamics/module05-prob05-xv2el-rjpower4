# ----------------------------------------------------------------------------------------
# Main Type Definition
# ----------------------------------------------------------------------------------------
"""
    KeplerianElements
"""
struct KeplerianElements <: FieldVector{6, Float64}
    sma::Float64
    ecc::Float64
    inc::Float64
    aop::Float64
    raan::Float64
    ta::Float64

    # Inforce structure invariants on valid Keplerian element values
    function KeplerianElements(sma, ecc, inc, aop, raan, ta)
        if sma == 0
            throw(ArgumentError("zero semi-major axis"))
        end

        if ecc < 0.0
            throw(ArgumentError("negatie eccentricity"))
        end

        if ecc == 1 && sma != Inf
            throw(ArgumentError("non-infinite semi-major axis with unit eccentricity"))
        end

        if ecc > 1 && sma > 0.0
            throw(ArgumentError("non-negative semi-major axis with eccentricity > 1"))
        end

        new(sma, ecc, inc, aop, raan, ta)
    end
end

# ----------------------------------------------------------------------------------------
# Extra External Constructors
# ----------------------------------------------------------------------------------------
function KeplerianElements(sma; ecc=0.0, inc=0.0, aop=0.0, raan=0.0, ta=0.0)
    return KeplerianElements(sma, ecc, inc, aop, raan, ta)
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
semi_minor_axis(k::KeplerianElements) = semi_minor_axis(semi_major_axis(k), eccentricity(k))
semi_latus_rectum(k::KeplerianElements) = semi_latus_rectum(semi_major_axis(k), eccentricity(k))
periapsis_radius(k::KeplerianElements) = periapsis_radius(semi_major_axis(k), eccentricity(k))
apoapsis_radius(k::KeplerianElements) = apoapsis_radius(semi_major_axis(k), eccentricity(k))
apsis_radii(k::KeplerianElements) = apsis_radii(semi_major_axis(k), eccentricity(k))

# ----------------------------------------------------------------------------------------
# Standard Two-Body Dynamical Relations
# ----------------------------------------------------------------------------------------
mean_motion(gm, k::KeplerianElements) = mean_motion(gm, semi_major_axis(k))
period(gm, k::KeplerianElements) = mean_motion(gm, semi_major_axis(k))
specific_energy(gm, k::KeplerianElements) = specific_energy(gm, semi_major_axis(k))

function specific_angular_momentum(gm, k::KeplerianElements) 
    return specific_angular_momentum(gm, semi_major_axis(k), eccentricity(k))
end

function radius(gm, k::KeplerianElements; ta=true_anomaly(k))
    return radius(gm, semi_major_axis(k), eccentricity(k), ta)
end

function velocity(gm, k::KeplerianElements; ta=true_anomaly(k))
    return velocity(gm, semi_major_axis(k), eccentricity(k), ta)
end

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

    r = radius(gm, ke)

    alpha = (sma * n) / sqrt(1 - ecc^2)

    v_x = - alpha * sin(ta)
    v_y = alpha * (ecc + cos(ta))
    v_z = zero(v_x)

    x = r * cos(ta)
    y = r * sin(ta)
    z = zero(x)

    return [x, y, z, v_x, v_y, v_z]
end