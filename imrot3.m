 function image = imrot3(im, ang, inttype)
 %   image = imrot2(im, angle, inttype);
 %     im = input image (must be square)
 %     angle = rotation angle in radians
 %     imttype = 'bilinear', etc.
 % This function is similar to the built-in function imrotate
 % but the point of rotation is different here: (N/2+1, N/2+1)
 % and angles are in radians

L = size(im,1);
v = [-L/2:L/2-1];
[xx,yy] = meshgrid(v,v);
xr = xx.*cos(ang) - yy.*sin(ang);
yr = xx.*sin(ang) + yy.*cos(ang);
image = interp2(v,v,im,xr,yr,inttype);
