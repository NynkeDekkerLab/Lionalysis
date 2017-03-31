function [f2,DwellTime] = PlotDT(Traj)
%PLOTDT Dwell time analysis
%   Detailed explanation goes here

OSslash='/';

display('Select the Data folder')
dataDirectory = uigetdir(pwd,'Select Data Directory');
dataDirectory = sprintf('%s%s', dataDirectory,OSslash);

if nargin<1,
    Traj=load(sprintf('%s%s%s%s',dataDirectory,'Results',OSslash,'Traj.mat'));
    Traj=Traj.Traj;
end

% Exposure time and interval between frames (in s)
ExposureTime=1;

FrameTime=6;

SizeT=size(Traj,2);

for i=1:SizeT
    if size(Traj{i})>1
    DwellTime(i)=ExposureTime+(size(Traj{i},1)-1)*FrameTime; % (in s)
    else 
    DwellTime(i)=0;
    end
end

%Removal of spots that are always present
DwellTime(DwellTime==max(DwellTime))=0;
DwellTime=nonzeros(DwellTime);

Nbins=15;

[histdata,XDT]=hist(DwellTime,Nbins+1);
x=linspace(min(XDT),max(XDT),Nbins+1)';
X=linspace(0,max(XDT),1000)';

F=@(t,a,Tm) a*exp(Tm*t);

f2=fit(x,histdata','exp1');

%Correction for fluorophore bleaching
Tbl=(FrameTime/ExposureTime)*7.1;

C=coeffvalues(f2);

TotalDecayTime=1/(abs(C(2)))

UnloadingTime=1/(1/TotalDecayTime-1/Tbl)

figure(1)
hist(DwellTime,Nbins)
hold on
k=plot(X,F(X,C(1),-1/UnloadingTime),'--g','LineWidth',3);
m=plot(X,F(X,C(1),C(2)),'-r','LineWidth',3);
pp=plot(X,F(X,C(1),-1/Tbl),'--b','LineWidth',3);
% set(l(1),'LineWidth',3);
axis([5 100 0 700]);
hold off
legend('data', 'unloading curve', 'total decay', 'bleaching');
% H = fitdist(DwellTime','exponential');
% histfit(DwellTime,Nbins,'exponential')
% set(gca,'FontSize',18)
title(strcat({'Tus Dwell Time OriZ strain - t_{unload} = '},num2str(UnloadingTime),{' N_{tracks} : '},num2str(size(Traj,2))));
xlabel('Dwell Time (s)')
ylabel('Frequency (-)')
% axis([0 30 0 ])
set(gca,'FontSize',16)




end

