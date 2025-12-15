
"""
Rotational matrix about axis-1
"""
function _rotmat_ax1(phi::Real)
    return [1.0 0.0 0.0; 0.0 cos(phi) sin(phi); 0.0 -sin(phi) cos(phi)]
end


"""
Rotational matrix about axis-2
"""
function _rotmat_ax2(phi::Real)
    return [cos(phi) 0.0 -sin(phi); 0.0 1.0 0.0; sin(phi) 0.0 cos(phi)]
end


"""
Rotational matrix about axis-3
"""
function _rotmat_ax3(phi::Real)
    return [cos(phi) sin(phi) 0.0; -sin(phi) cos(phi) 0.0; 0.0 0.0 1.0]
end
