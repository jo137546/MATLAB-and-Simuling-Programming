function C = CoordinateTransform(axis, angle)
% CoordinateTransform
% Returns a principal coordinate-frame rotation matrix.
%
% axis  = 'x', 'y', or 'z'
% angle = rotation angle in degrees

c = cosd(angle);
s = sind(angle);

switch lower(axis)

    case 'x'
        C = [1  0  0;
             0  c  s;
             0 -s  c];

    case 'y'
        C = [ c  0 -s;
              0  1  0;
              s  0  c];

    case 'z'
        C = [ c  s  0;
             -s  c  0;
              0  0  1];

    otherwise
        error('Axis must be x, y, or z.');

end
end