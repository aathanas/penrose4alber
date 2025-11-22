clc; close all; clear all

%% Demo stability analysis for a single Gaussian

sigma = 0.36;
C=1;
P =@(k)  C^2 * exp(-pi *k.^2 / sigma^2) / sigma;


X=linspace(4e-3,1.5,160);
Cutoff = 3;


thisX= min(X);
DxP = @(t) ( P(t+min(thisX)/2) - P(t-min(thisX)/2) ) / min(thisX);

kk=linspace(-Cutoff,Cutoff,999);


fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

subplot(3,1,1)
plot(kk,P(kk))
hold on
grid on
xlabel('k')
ylabel('P(k)')
title('Gaussian spectrum, P')
set(gca,'Fontsize',19)



subplot(3,1,2)
semilogy(kk,P(kk))
hold on
grid on
xlabel('k')
ylabel('P(k)')
title('log plot of P')
set(gca,'Fontsize',19)




subplot(3,1,3)

thisX= min(X);
DxP = @(t) ( P(t+min(thisX)/2) - P(t-min(thisX)/2) ) / min(thisX);
plot(kk,DxP(kk),'DisplayName',['D_XP for X = ' num2str(thisX)])

hold on

thisX= mean(X);
DxP = @(t) ( P(t+min(thisX)/2) - P(t-min(thisX)/2) ) / min(thisX);
plot(kk,DxP(kk),'DisplayName',['D_XP for X = ' num2str(thisX)])

thisX= max(X);
DxP = @(t) ( P(t+min(thisX)/2) - P(t-min(thisX)/2) ) / min(thisX);
plot(kk,DxP(kk),'DisplayName',['D_XP for X = ' num2str(thisX)])

grid on
xlabel('k')
ylabel('D_XP(k)')
title('Divided difference with increment X, D_XP')

legend

set(gca,'Fontsize',19)

drawnow

[Xtmp,f1,S]=CheckPenroseCondition(P,4*pi,1,X,Cutoff);

if ~isempty(Xtmp)
    disp(['For this Gaussian spectrum, the bandwidth of Penrose-unstable wavenumbers is ' num2str(2*Xtmp) '.'])
else
    disp('This spectrum is Penrose-stable.')
end



%% Stability analysis for Gaussians with progressively larger intensity


Cvec = linspace(0.9,2.2,81);


tic
for ii=1:length(Cvec)

    C=Cvec(ii);
    P =@(k)  C^2 * exp(-pi *k.^2 / sigma^2) / sigma;


    Xtmp=CheckPenroseCondition(P,4*pi,0,X,Cutoff);

    if isempty(Xtmp)
        Xtmp=0;
    end
    Xstar(ii) = Xtmp;


end
toc

f2 = figure;
plot(Cvec, 2*Xstar,'k+-','MarkerSize',8,'LineWidth',1.5)
ylabel('B')
hold on
xlabel('C')
set(gca,'FontSize',19)


% post-processing, quadratic relation of unstable bandwidth and intensity

p = polyfit(2*Xstar, Cvec, 2); % Returns [a, b, c] for Y = a*X^2 + b*X + c




Y_fit = polyval(p, 2*Xstar);

plot(Y_fit, 2*Xstar, 'r-', 'LineWidth', 2)

grid on

legend('Computed bandwidth of unstable wavenumbers, B','Quadratic fit C \approx a_2B^2 + a_1B + a_0','interpreter','tex','Location','southeast')

fprintf('Quadratic coefficients: a = %.4f, b = %.4f, c = %.4f\n', p(1), p(2), p(3));

saveas(f2,'gaussian_UnstBand','png')