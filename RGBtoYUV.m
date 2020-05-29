function YUV = RGBtoYUV(RGB)
%
% YUV and RGB are Nx3 uint8 matrices
%

% First convert RGB to range 0.0 to 1.0, and add a column of ones.
RGB1 = [double(RGB)/255, ones(size(RGB,1),1)];

% Define color transform.
Q = [0.29899999,    -0.1687,         0.5,
     0.587,         -0.3313,        -0.4187,
     0.114,          0.5,           -0.0813,
     0,              0.50196078,	 0.50196078];

% Do the transform.
YUV = RGB1 * Q;

% Clip to range 0.0 to 1.0, and convert to uint8.
%YUV = uint8(255*max(min(YUV,ones(size(YUV))),zeros(size(YUV))));
%eduardo, keep it double
YUV = (255*max(min(YUV,ones(size(YUV))),zeros(size(YUV))));