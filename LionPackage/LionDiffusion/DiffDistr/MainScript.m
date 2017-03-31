clear all
close all

display('Loading set file ... select the Data folder')
opt=VB3_getOptions('General_set');
%Folder with data
display(sprintf('%s%s','JobID: ',opt.jobID))
load(sprintf('%s%s',opt.dataDirectory,'Results/D.mat'))

Nbins=400;

% ----
lb=0;
ub=2.5;
x=linspace(0,ub,Nbins);
n=4; 
A=0.5;
% ---

D1_constant=0.01; %immobile D
D2_constant=0.2; %mobile D

% -- Control PDFs (one component, two component 50/50

%Single component
pdfone=@(x,D) ((n/D)^n*x.^(n-1).*exp(-n*x/D))/(factorial(n-1));

%Double component (50/50)
pdftwo=@(x,D1,D2) (A*(n/D1)^n*x.^(n-1).*exp(-n*x/D1))/(factorial(n-1))+(1-A)*((n/D2)^n*x.^(n-1).*exp(-n*x/D2))/(factorial(n-1));

% -- Important PDFs

% variables : D1, D2, B (fraction)
pdftwofrac=@(x,D1,D2,B) (B*(n/D1)^n*x.^(n-1).*exp(-n*x/D1))/(factorial(n-1))+(1-B)*((n/D2)^n*x.^(n-1).*exp(-n*x/D2))/(factorial(n-1));

%variables : B (immobile fraction)
pdffrac=@(x,B) (B*(n/D1_constant)^n*x.^(n-1).*exp(-n*x/D1_constant))/(factorial(n-1))+(1-B)*((n/D2_constant)^n*x.^(n-1).*exp(-n*x/D2_constant))/(factorial(n-1));


% Control PDF MLE
[phattwo,pci1]=mle(D,'pdf',pdftwo,'start',[0.01,0.2]);
[phatone,pci2]=mle(D,'pdf',pdfone,'start',0.01);

% Important PDF MLE
% Starting values can induce error!
[phattwofrac,pci3]=mle(D,'pdf',pdftwofrac,'start',[0.01,0.2,0.01]); %[D1,D2,B]
[phatfrac,pci4]=mle(D,'pdf',pdffrac,'start',0.01); % B


Art=pdftwo(x,phattwo(1),phattwo(2));
Brt=pdfone(x,phatone(1));
Crt=pdftwofrac(x,phattwofrac(1),phattwofrac(2),phattwofrac(3));

Drt=pdfone(x,phattwofrac(1));
Ert=pdfone(x,phattwofrac(2));

Frt=pdffrac(x,phatfrac(1));

FracLow=phattwofrac(3);
FracHigh=1-FracLow;

%% Goodness of fit
bins=linspace(lb, ub, Nbins);
Dc=histc(D,bins);


%% Plotting options for available models

hold on
fig1=figure(1);
h=bar(bins,Dc*100/(sum(Dc)));
h.EdgeColor='k';
h.FaceColor=[.5 .5 .5];
set(gca,'fontsize',18)
axis([-0.01 1.5 0 7])
xlabel('Diffusion Coefficient (um^{2}/s)')
ylabel('Probability (-)')


% Single D Component model
% plot(x,Brt/sum(Brt)*100,'b','LineWidth',2);
% t=annotation('textbox',[0.45 0.2 0.2 0.07],...
%     'String',strcat({'D_{single}= '},sprintf('%.3f',phatone(1))),...
%     'FontName','Arial');

% Double Component model, with fraction immobile A=50%
% plot(x,Art/sum(Art)*100,'r','LineWidth',2);
% t=annotation('textbox',[0.45 0.2 0.2 0.07],...
%     'String',strcat({'D_{fast} = '},sprintf('%.3f',max(phattwo(1),phattwo(2)))),...
%     'FontName','Arial');
% s=annotation('textbox',[0.2 0.8 0.2 0.07],...
%     'String',strcat({'D_{slow} = '},sprintf('%.3f',min(phattwo(1),phattwo(2)))),...
%     'FontName','Arial');
%     
% SlowP=annotation('textbox',[0.17 0.6 0.1 0.1],...
%     'String',strcat(sprintf('%2.0f',A*100),'%'),...
%     'FontName','Arial',...
%     'EdgeColor','w','FontSize',25,'Color','r','FontWeight','bold');
% FastP=annotation('textbox',[0.3 0.2 0.1 0.1],...
%     'String',strcat(sprintf('%2.0f',(1-A)*100),'%'),...
%     'FontName','Arial',...
%      'EdgeColor','w','FontSize',25,'Color','b','FontWeight','bold');

% Double component D1,D2,B variable. SEPARATED

plot(x,phattwofrac(3)*Drt/sum(Drt)*100,'r','LineWidth',2);
plot(x,(1-phattwofrac(3))*Ert/sum(Ert)*100,'b','LineWidth',2);
t=annotation('textbox',[0.45 0.2 0.2 0.07],...
    'String',strcat({'D_{fast} = '},sprintf('%.3f',max(phattwofrac(1),phattwofrac(2)))),...
    'FontName','Arial');
s=annotation('textbox',[0.2 0.8 0.2 0.07],...
    'String',strcat({'D_{slow} = '},sprintf('%.3f',min(phattwofrac(1),phattwofrac(2)))),...
    'FontName','Arial');
    
SlowP=annotation('textbox',[0.17 0.6 0.1 0.1],...
    'String',strcat(sprintf('%2.0f',phattwofrac(3)*100),'%'),...
    'FontName','Arial',...
    'EdgeColor','w','FontSize',25,'Color','r','FontWeight','bold');
FastP=annotation('textbox',[0.3 0.2 0.1 0.1],...
    'String',strcat(sprintf('%2.0f',(1-phattwofrac(3))*100),'%'),...
    'FontName','Arial',...
    'EdgeColor','w','FontSize',25,'Color','b','FontWeight','bold');

% Double component D1,D2,B variable. COMBINED
%  plot(x,Crt/sum(Crt)*100,'r','LineWidth',2);
% t=annotation('textbox',[0.45 0.2 0.2 0.07],...
%     'String',strcat({'D_{fast} = '},sprintf('%.3f',max(phattwofrac(1),phattwofrac(2)))),...
%     'FontName','Arial');
% s=annotation('textbox',[0.2 0.8 0.2 0.07],...
%     'String',strcat({'D_{slow} = '},sprintf('%.3f',min(phattwofrac(1),phattwofrac(2)))),...
%     'FontName','Arial');
%     
% SlowP=annotation('textbox',[0.17 0.6 0.1 0.1],...
%     'String',strcat(sprintf('%2.0f',phattwofrac(3)*100),'%'),...
%     'FontName','Arial',...
%     'EdgeColor','w','FontSize',25,'Color','r','FontWeight','bold');
% FastP=annotation('textbox',[0.3 0.2 0.1 0.1],...
%     'String',strcat(sprintf('%2.0f',(1-phattwofrac(3))*100),'%'),...
%     'FontName','Arial',...
%     'EdgeColor','w','FontSize',25,'Color','b','FontWeight','bold');

% Fraction B is variable, D1 and D2 set.
% plot(x,Frt/(sum(Frt))*100,'m','LineWidth',2);
% t=annotation('textbox',[0.45 0.2 0.2 0.07],...
%     'String',strcat({'D_{fast} = '},sprintf('%.3f',D2_constant)),...
%     'FontName','Arial');
% s=annotation('textbox',[0.2 0.8 0.2 0.07],...
%     'String',strcat({'D_{slow} = '},sprintf('%.3f',D1_constant)),...
%     'FontName','Arial');
%     
% SlowP=annotation('textbox',[0.17 0.6 0.1 0.1],...
%     'String',strcat(sprintf('%2.0f',phatfrac(1)*100),'%'),...
%     'FontName','Arial',...
%     'EdgeColor','w','FontSize',25,'Color','r','FontWeight','bold');
% FastP=annotation('textbox',[0.3 0.2 0.1 0.1],...
%     'String',strcat(sprintf('%2.0f',(1-phatfrac(1))*100),'%'),...
%     'FontName','Arial',...
%     'EdgeColor','w','FontSize',25,'Color','b','FontWeight','bold');

%%
title(strcat({opt.jobID },{' n_{tracks}: '},num2str(size(D,1))));


t.FontSize = 16;
s.FontSize = 16;
t.FontWeight='bold';
s.FontWeight='bold';


