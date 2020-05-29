function RGB = YUVtoRGB(YUV)
%
% RGB and YUV are Nx3 uint8 matrices
%

% First convert YUV to range 0.0 to 1.0, and add a column of ones.
YUV1 = [double(YUV)/255, ones(size(YUV,1),1)];

% Define color transform.
M = [1,             1,          1,
     0,            -0.34414,    1.772,
     1.402,        -0.71414,    0,
    -0.703749019,  0.53121505, -0.88947451];

% Do the transform.
RGB = YUV1 * M;

% Clip to range 0.0 to 1.0, and convert to uint8.
RGB = uint8(255*max(min(RGB,ones(size(RGB))),zeros(size(RGB))));