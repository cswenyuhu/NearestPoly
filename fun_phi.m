clear
clc
close all

q=[0.5,1,1.5];   
t=0:0.2:2;
for i=1:length(q)
    p=q(i);
    x(i,:)=2/(2-p)*t-p/(2-p)*t.^(2/p)-1;
    figure(1); grid on; hold on; 
    if i==1
        plot(t,x(i,:),'rs-','LineWidth',1);
    elseif i==2
            plot(t,x(i,:),'kd-.','LineWidth',1);
    else i==3
        plot(t,x(i,:),'bo--','LineWidth',1);
    end
    hold off
end