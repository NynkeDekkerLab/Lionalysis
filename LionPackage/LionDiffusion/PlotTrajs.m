function [] = PlotTrajs()
%PLOTTRAJS Summary of this function goes here
%   Detailed explanation goes here
x = input('Chopped Trajectories? (1=yes, 0=no)');
Dthreshold=input('D threshold immobile to mobile? (enter: Dthres=0.05) ');

if isempty(x);
    x=1;
end

if isempty(Dthreshold)
    Dthreshold=0.05;
end


display('Loading set file ... select the Data folder')
opt=VB3_getOptions('General_set');
%Folder with data
display(sprintf('%s%s','JobID: ',opt.jobID))
if x==1;
D=load(sprintf('%s%s',opt.dataDirectory,'Results/D.mat'));
Traj=load(sprintf('%s%s',opt.dataDirectory,'Results/TrajChopped.mat'));
Traj=Traj.TrajC;
D=D.D;
else
TrajOri=load(sprintf('%s%s',opt.dataDirectory,'Results/Traj.mat'));
Traj=TrajOri.Traj;
D=Dcalc(Traj,opt.timestep);
end

%How does D look?
Nbins=100; bins=linspace(0, 2, Nbins); Dc=histc(D,bins); 
figure(1)
hold on
bar(bins,Dc/(sum(Dc)))
l=line([Dthreshold Dthreshold],[0 1])
set(l,'LineWidth',2,'color','r');
hold off
axis([0 1 0 0.1])
set(gca,'FontSize',16);
title('Distribution of Diffusion Constants')
xlabel('Diffusion Constant (um2/s)');
ylabel('Probability (-)');


Traj=Traj(~cellfun('isempty',Traj));

SizeT=size(Traj,2);
Tfast=cell(1,SizeT);
Tslow=cell(1,SizeT);

Ps=0.159; % PixelSize 0.159 micron per pixel

figure(2)
hold on
for i=1:SizeT   
    Traj{i}=Traj{i}./Ps;
    
    if D(i)>Dthreshold
        Tfast{i}=Traj{i};
    else
        Tslow{i}=Traj{i};
    end
    
end

for i=1:SizeT
    if ~isempty(Tfast{i});
    plot(Tfast{i}(:,1),Tfast{i}(:,2),'b','LineWidth',2);
    end
end

for i=1:SizeT   
    if ~isempty(Tslow{i});
    plot(Tslow{i}(:,1),Tslow{i}(:,2),'r','LineWidth',2);
    end
end
hold off
title(strcat({'N_{tracks}: ' },num2str(SizeT),'{ D_{thres}: }',num2str(Dthreshold)));
set(gca,'FontSize',16);
xlabel('X Position (-)');
ylabel('Y Position (-)');
end

