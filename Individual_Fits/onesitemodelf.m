function ndh=onesitemodelf(par,xmt,mt,injv,V0,SyC,I0)
% Model Function for fitting ITC data using One set of Sites
n=par(1);Ka=par(2);dH=par(3);cac=par(4);

%Converting items to SI Units
injv=injv*1e-6;
mt=mt*1e-3;

b=1+xmt./(n)+1./(n*Ka.*mt);
Qv = (mt.*(n*dH*V0/2)).*(b-sqrt((b.*b)-(4/n).*(xmt)));

% Two lines below calculate Q for first injection, not used in the fit
b0=1+(I0(3)/n)+1/(n*Ka*I0(2));
Q0 = (I0(2)*(n*dH*V0/2))*(b0-sqrt((b0*b0)-(4/n)*(I0(3))));

qfit=Qv(1)+(I0(1)/V0)*0.5*(Qv(1)+Q0)-Q0;
for i=2:length(xmt)
    qfit=vertcat(qfit,Qv(i)+(injv(i)/V0)*((Qv(i)+Qv(i-1))/2)-Qv(i-1));
end

ndh=[];
for i=1:length(qfit)
    ndh=vertcat(ndh,qfit(i)/(injv(i)*SyC)); % Divide by number of moles of Xt injected
end

ndh=ndh+cac; % Added by deb to correct for shift

end