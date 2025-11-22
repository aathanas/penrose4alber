function [h,x_internal] = AdaptiveSpectralHilbertTransformDxP(P,X, t, tol)
arguments
    P function_handle
    X double
    t double
    tol (1,1) double {mustBePositive} = 1e-4 % Optional, default 1e-6
end

% Compute the continuous Hilbert transform on a specified mesh
% f1 is the function handle for the function of which we want the 
% Hilbert transform
% x is the grid of points on which the Hilbert transform is required
% tol (optional) is the pointwise error tolerance to aim for
% x_internal (optional) retuns the sampling points used to compute
% the Hilbert transforms





N=500; % starting number of sampling points
L=10; % scaling parameter controlling how far out the sampling points reach

% [hold,x_internal] = hilb1(f1, N, L, x);
[hold,x_internal] = hilbDxP(P,X, N, L, t);

err=1;
% error initialization

while err>tol


    L=L*sqrt(2); N=N*2;
    % sampling points are extended and refined at the same time

    [h,x_internal] = hilbDxP(P,X, N, L, t);
    % new guess for the requested Hilbert transform values

    err = max(abs(h-hold));
    % difference  of old values and new values is taken as a proxy for
    % error

    hold=h;

end



end