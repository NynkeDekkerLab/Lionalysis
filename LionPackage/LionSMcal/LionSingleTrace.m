clear all
close all

init.datapath = uigetdir(pwd,'Select Data Folder');
init.OSslash = '/';
init.datapath = strcat(init.datapath,init.OSslash);

[X]=xlsread(strcat(init.datapath,'Spot1_values.xls'));

X(:,2)=X(:,2)/100;

GolayFilt=sgolayfilt(X(1:100,2),30,99);

figure(1)
hold on
plot(X(1:40,2),'-k+')
plot(GolayFilt(1:40),'r','LineWidth',2);
hold off
set(gca,'FontSize',16)
xlabel('Time (s)')
ylabel('Average Integrated Intensity Counts (-)')

