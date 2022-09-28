"""
Basic two-body functionality.
"""

"""
    semi_minor_axis(sma, ecc)
"""
semi_minor_axis(sma, ecc) = sma * sqrt(1 - ecc^2)

"""
    mean_motion(gm, sma)
"""
mean_motion(gm, sma) = sqrt(gm / sma^3)

"""
    specific_energy(gm, sma)
"""
specific_energy(gm, sma) = -gm / (2sma)

"""
    specific_energy(gm, r, v)
"""
specific_energy(gm, r, v) = (1//2) * v^2 - gm / r

"""
    specific_angular_momentum(gm, sma, ecc)
"""
specific_angular_momentum(gm, sma, ecc) = sqrt(gm * sma * (1 - ecc^2))

"""
    semi_latus_rectum(sma, ecc)
"""
semi_latus_rectum(sma, ecc) = sma * (1 - ecc^2)

"""
    periapsis_radius(sma, ecc)
"""
periapsis_radius(sma, ecc) = sma * (1 - ecc)

"""
    apoapsis_radius(sma, ecc)
"""
apoapsis_radius(sma, ecc) = sma * (1 + ecc)

"""
    apsis_radii(sma, ecc)
"""
apsis_radii(sma, ecc) = (sma - sma * ecc, sma + sma * ecc)

"""
    period(gm, sma)
"""
period(gm, sma) = 2pi / mean_motion(gm, sma)

"""
    radius(gm, sma, ecc, ta)
"""
radius(gm, sma, ecc, ta) = semi_latus_rectum(sma, ecc) / (1 + ecc * cos(ta))

"""
    velocity(gm, sma, ecc, ta)

Compute velocity magnitude on specified orbit.
"""
function velocity(gm, sma, ecc, ta)
    e = specific_energy(gm, sma)
    r = radius(gm, sma, ecc, ta)
    return sqrt(2 * (e + gm / r))
end

"""
    perifocal_to_inertial_rotation_matrix(aop, inc, raan)
"""
function perifocal_to_inertial_rotation_matrix(aop, inc, raan)
    sa, ca = sincos(aop)
    si, ci = sincos(inc)
    sr, cr = sincos(raan)

    return [
        (cr*ca-sr*sa*ci) (-cr*sa-sr*ca*ci) (sr*si)
        (sr*ca+cr*sa*ci) (-sr*sa+cr*ca*ci) (-cr*si)
        (sa*si) (ca*si) ci
    ]
end

"""
    inertial_to_perifocal_rotation_matrix(aop, inc, raan)
"""
function inertial_to_perifocal_rotation_matrix(aop, inc, raan)
    return transpose(perifocal_to_inertial_rotation_matrix(aop, inc, raan))
end

