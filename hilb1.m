function [h,x] = hilb1(F, N, L, y)

%   The function h = hilb1(F, N, L, y) computes the Hilbert transform
%   of a function F(x) defined on the real line, at specified
%   values of y (y could be a scalar, vector, or matrix.)

%   The method is based on an expansion in rational eigenfunctions
%   of the Hilbert transform, described in the reference below.
%   The total number of sampling points in the expansion is 2N, and L is
%   a scalar parameter controlling how far out the sampling points reach.

%   Reference:  J.A.C. Weideman, "Computing the Hilbert Transform on
%               the Real Line", Math. Comp. Vol. 64, pp. 745--762 (1995)
 
%  File modified from https://appliedmaths.sun.ac.za/~weideman/research/hilbert.html
%  with permission from the author


n  = [-N:N-1]';                        %  Compute the collocation points ...
x  = L*tan(pi*(n+1/2)/(2*N));
FF = F(x);                          %  ... and sample the function.



a  = fft(fftshift(FF.*(L-i*x)));       %  These three lines compute the
a  = exp(-i*n*pi/(2*N)).*fftshift(a);  %  expansion coefficients.
a  = flipud(i*(sign(n+1/2).*a))/(2*N);

z  = (L+i*y)./(L-i*y);                 %  The evaluation of the transform
h  = -polyval(a,z)./(z.^N.*(L-i*y));    %  reduces to polynomial evaluation
                                       %  in the variable z.















