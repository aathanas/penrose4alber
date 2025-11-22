function [h,x_internal] = AdaptiveSpectralHilbertTransform(f1, x, tol)
arguments
    f1 function_handle
    x double
    tol (1,1) double {mustBePositive} = 1e-3 % Optional, default 1e-6
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

[hold,x_internal] = hilb1(f1, N, L, x);
% initial guess for the requested Hilbert transform values

err=1;
% error initialization

while err>tol 

    if N > 64000
        disp(['Hilbert transform error is ' num2str(err)])
        err = 0;
    else
        L=L*sqrt(2); N=N*2;
        % sampling points are extended and refined at the same time

        [h,x_internal] = hilb1(f1, N, L, x);
        % new guess for the requested Hilbert transform values

        err = max(abs(h-hold));
        % difference  of old values and new values is taken as a proxy for
        % error

        hold=h;
    end

end



end