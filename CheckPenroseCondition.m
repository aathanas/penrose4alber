function  [Xstar,fig1,S]=CheckPenroseCondition(P,target,plot_flag,X,Cutoff)
arguments
    P function_handle = @(x) 4*exp(-24*x.^2) % function handle for the background spectrum
    target  double = 4*pi    % the target point
    plot_flag logical = 1   % 1 for a figure, 0 for no figure
    X double = linspace(4e-3,1.5,160)   % initial vector of wavenumbers to check for stability
    Cutoff double = 9 % cutoff of the effective support of the spectrum
end

% Performs a check of the Penrose stability condition for the Alber
% equation, for the spectrum P by checking if H[D_X P](z) = 4*pi *p/q
% has solutions. See Athanassoulis et al., "Strong solutions for the
% Alber equation and stability of unidirectional wave spectra" (2020)
% Kinetic and Related Models, 13(4), pp.703-737.
%
% Run with no input arguments for a basic demo.
%
% OUTPUT ARGUMENTS:
% Xstar: approximate largest unstable wavenumber. Returns empty array if
% the spectrum is stable.
% fig1: handle of figure (if produced)
% S: points of the signal transform curve

Cutoff =Cutoff + max(abs(X));


if plot_flag
    fig1=figure;
    plot(target,0,'bs','MarkerSize',9,'MarkerFaceColor','b','DisplayName','target point')
    grid on
    hold on
    set(gca,'FontSize',19)
end



tt=linspace(-Cutoff,Cutoff,10000); % points used in each curve S(tt) = H[D_X P(tt)] - i D_X P(tt). 
                                   % Add more points if the curves don't look well resolved. 
                                   % Use fewer points if the curves are well resolved bu the overall computation too slow.

unstables=[]; % vector of unstable wavenumbers, if any

for ii=1:length(X)
    DxP = @(t) ( P(t+X(ii)/2) - P(t-X(ii)/2) ) / X(ii);


    h = AdaptiveSpectralHilbertTransform(DxP,tt);

    if sum(isnan(h))
        error('Hilbert transform yields NaN')
    end
    



    S= h - 1i * DxP(tt);
    S=[S S(1)];

    if inpolygon(target,0,real(S),imag(S))
        unstables=[unstables X(ii)];
        if plot_flag
            plot(real(S),imag(S),'r-','LineWidth',0.5,'DisplayName','curve for unstable wavenumber')

        end
    else
        if plot_flag
            plot(real(S),imag(S),'k-','LineWidth',0.5,'DisplayName','curve for unstable wavenumber')
        end
    end


end


Xstar = max(unstables);


if ~isempty(Xstar) & plot_flag
    xlabel('Curves for unstable wavenumbers are in red.')
end


end