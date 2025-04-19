clear all
close all



##RossbyRadius=gH*/((f*dx)^2);
##
##wfA=(1+(g*H/(f^2))*(sin(k*dx)/(dx^2)))^2


RossbyRadius=4;

kdx=[0:0.01:pi];

for ii=1:length(kdx)
wfAnalytic(ii)=(1+RossbyRadius*kdx(ii)^2)^(1/2);
wfA(ii)=(1+RossbyRadius*(sin(kdx(ii)))^2)^(1/2);
wfB(ii)=(1+4*RossbyRadius*(sin(kdx(ii)/(2)))^2)^(1/2);
wfC(ii)=((cos(kdx(ii)/2))^2+4*RossbyRadius*(sin(kdx(ii)/(2)))^2)^(1/2);
end

figure(1)
plot(kdx,wfA,'-b',kdx,wfB,'-r',kdx,wfAnalytic,'-k',kdx,wfC,'-g')
