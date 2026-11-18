function C = eci2lvlh(r_bar,v_bar)
% ECI -> radial/transverse/normal axes used in Oct. 13 relative-motion notes.
% x: radially outward; y: transverse in direction of motion; z: r cross v.
% For an ellipse, y is not exactly parallel to the full velocity vector.
r_bar = r_bar(:); v_bar = v_bar(:);
ex = r_bar/norm(r_bar);
h = cross(r_bar,v_bar);
ez = h/norm(h);
ey = cross(ez,ex);
C = [ex.';ey.';ez.'];
end

