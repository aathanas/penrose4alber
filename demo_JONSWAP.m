clc; close all; clear all

%% Stability analysis for a single JONSWAP spectrum

P = ConstructJONSWAP(3e-2,5);

DxP = @(t) ( P(t+0.01/2) - P(t-0.01/2) ) / 0.01;

Cutoff = 40;

kk=linspace(-Cutoff,Cutoff,999);

X = linspace(4e-3,1.5,70);


fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

subplot(3,1,1)
plot(kk,P(kk))
hold on
grid on
xlabel('k')
ylabel('P(k)')
title('JONSWAP spectrum, P')
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




[Xtmp,f1,S]=CheckPenroseCondition(P,1/(4*pi),1,X,Cutoff );




%% Stability analysis for different values of alphs, gamma



alpha = linspace(9e-3,8e-2,30);
gamma = linspace(1,5,30);



tic
for ii=1:length(alpha)
    for jj=1:length(gamma)

        P = ConstructJONSWAP(alpha(ii),gamma(jj));

        Xtmp=CheckPenroseCondition(P,1/(4*pi),0,X,Cutoff );

        if isempty(Xtmp)
            Xtmp=0;
        end
        Xstar(ii,jj) = Xtmp;

    end
end
toc

f1 = figure;
pcolor(gamma,alpha,Xstar')
shading interp
colorbar
title('Unstable bandwidth')
xlabel('\gamma')
ylabel('\alpha')
set(gca,'Fontsize',19)

saveas(f1,'JSPpcolorUB','png')





f2 = figure;
contour(gamma,alpha,Xstar')
shading interp
colorbar
title('Unstable bandwidth')
xlabel('\gamma')
ylabel('\alpha')
set(gca,'Fontsize',19)

saveas(f2,'JSPcontourUB','png')
