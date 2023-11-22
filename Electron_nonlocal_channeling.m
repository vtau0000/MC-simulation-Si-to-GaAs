amu=1.67e-27;
Z1=14;
Z2=64;
e=1.6e-19;
M1=28*amu;
M2=144*amu;
Z2Ga=31;
Z2As=33;
M2Ga=69*amu;
M2As=74*amu;

NN=5.32e-24*2/(M2*1000);%number density in (A^3)
N=100000;
lN=30;
E0=200000*e;

asa=rand(N,lN);
asa1=rand(N,lN);
idx=(asa<=0.9);%parameter to decide the ratio of dechanneling
idxatom=(asa1<0.5);%parameter to decide the Ga or As to be the substrate atom in collision
A=(idxatom.*M2Ga+~idxatom.*M2As)/M1;
Er=A.*E0./(1+A);%ZBL model relative energy
sela=0.9;
ll=1/((NN)^0.33);
pmax=1.35*(ll./((pi)^0.5));%maximum impact parameter
p=asa.^0.5.*pmax;
phii=2*pi*((rand(N,lN)+0.001)*0.99);
gamma=4*M1*M2/((M1+M2)^2);
gammaGa=4*M1*M2Ga/((M1+M2Ga)^2);
gammaAs=4*M1*M2As/((M1+M2As)^2);

ww=zeros(N,lN);
a0=0.529; %Bohr radius
a=0.8853*a0/((Z1^0.5+Z2^0.5)^0.666);
lnGa=(32.53*M2Ga*E0/(Z1*Z2Ga*(M1+M2Ga)*(Z1^0.23+Z2Ga^0.23)*e*1000));%
lnAs=(32.53*M2As*E0/(Z1*Z2As*(M1+M2As)*(Z1^0.23+Z2As^0.23)*e*1000));%
ln=(32.53*M2*E0/(Z1*Z2*(M1+M2)*(Z1^0.23+Z2^0.23)*e*1000));%reduce energy
f=1.45/(Z1^0.4);
amod=f*0.8853*a0*(idxatom.*1/((Z1^0.5+Z2Ga^0.5)^0.66)+~idxatom.*1/((Z1^0.5+Z2As^0.5)^0.66))./0.3;
ln11=amod.*Er.*(idxatom.*1/(Z1*(Z2Ga)*e)+(~idxatom.*1/(Z1*(Z2As)*e)));%reduce energy-1
lnmod=(lnGa+lnAs)/2;
r0=Er.*(1-p.^2)./(Z1*Z2*e);
preduce=p./amod;
hold on
%%%%%%
fid=fopen('D:\thz1\measuresample\SIMS\20200603\3-1-1.txt');
s=textscan(fid,'%f32 %f32');
fclose(fid);
x=s{1};
x1=s{2};
% fid1=fopen('D:\thz1\measuresample\ion implantation\SitoGaAs50kevTRIM_1E14.txt');
% s1=textscan(fid1,'%f32 %f32');
% fclose(fid1);
% x2=s1{1};
% x21=s1{2};
N1=300;
% [ff,xx]=hist(x1,N1);
% plot(x*10,x1,'g.');
% xlim([0 12000])
% ylim([1E16 1E21])
% set(gca, 'YScale', 'log')
% xlabel('depth(A)');
% %ylabel('lateral depth');
% ylabel('concentration(ion/cm^3)');
%%%%%%
C0=1e15;
% layercon1=x1(1,1)*((x(1,1))+(x(2,1)-x(1,1))*0.5);
% layercon2=x1(2,1)*((x(2,1)-x(1,1))*0.5+(x(3,1)-x(2,1))*0.5);
% layercon3=x1(3,1)*((x(3,1)-x(2,1))*0.5+(x(4,1)-x(3,1))*0.5);
% if x1(1,1)/x1(2,1)>1
%  C01=C0-(layercon1+layercon2+layercon3)*1e-7;
% else
%  C01=C0;
% end
% C02=(sum(x1.*(x(2,1)-x(1,1)))-(layercon1+layercon2+layercon3))/1e7;
%three order solution formula(Cardano)
aa=2.21./ln11;
bb=-8.5./ln11;
cc=1./ln11;
dd=-preduce.^2;
d3solutionp=(3*aa.*cc-bb.^2)./(3*aa.^2);
d3solutionq=(2*bb.^3-9*aa.*bb.*cc+27*a.^2*dd)./(27*aa.^3);
d3solutionln1=(-d3solutionq./2)+((d3solutionq./2)+(d3solutionp./3)).^0.5;
d3solutionln2=(-d3solutionq./2)-((d3solutionq./2)+(d3solutionp./3)).^0.5;
r0moliere=real(((1./ln11)-((1./ln11).^2-(((30./ln11)+4).*preduce.^2)).^0.5)./((15./ln11)+2));
r0moliere1=(d3solutionln1.^0.33+d3solutionln2.^0.33).*amod;
c=9e9*Z1*(Z2/2)*e^2./Er;
theta=asin((1+(2*ln11.*p).^2).^-0.5).*2;
theta1=pi-2*p.*real(asin((c+p.^2)./((c.^2+4*p.^2).^0.5).*r0)-asin(c./((c.^2+4*p.^2).^0.5)))./abs(p);
theta1moliere=pi-2*p.*real((asin((c+p.^2)./((c.^2+4*p.^2.*(0.35*exp(-0.3*(r0moliere1./amod))+0.55*exp(-1.2*(r0moliere1./amod))+0.1*exp(-6*(r0moliere1./amod)))).^0.5).*(r0moliere1./amod))-asin((c+p.^2)./((c.^2+4*p.^2).^0.5)))./abs(p));
phi=atan(sin(theta1moliere)./((M1./(idxatom.*M2Ga+~idxatom.*M2As))+cos(theta1moliere)));
re=1;
se=re*(1.21*Z1^(1/6)*(Z1*Z2/((Z1^(2/3)+Z2^(2/3))^1.5))*(1/((M1/amu)^0.5)))*(E0/e)^0.5;%Linherd electron stopping power
semod=re*sqrt(2)*(idxatom.*(1.21*Z1^(1/6)*(Z1*Z2Ga/((Z1^(2/3)+Z2Ga^(2/3))^1.5))*(1/((M1/amu)^0.5)))*(E0/e)^0.5+~idxatom.*(1.21*Z1^(1/6)*(Z1*Z2As/((Z1^(2/3)+Z2As^(2/3))^1.5))*(1/((M1/amu)^0.5)))*(E0/e)^0.5);


gammamod=(idxatom.*gammaGa+~idxatom.*gammaAs);
frac=idx.*(1-(gammamod.*max(sin(phi./2).^2)))+~idx.*1;%energy remaining 
% frac=(1-(gammamod.*max(sin(phi./2).^2)));
%after each collision(first term:collision,second term:channeling)
dE=zeros(N,lN);
u=zeros(N,lN);
v=zeros(N,lN);
xnl=zeros(N,lN);
xloc=zeros(N,lN);
dEnl=zeros(N,lN);
dEloc=zeros(N,lN);
snl=zeros(N,lN);
sn=zeros(N,lN);
dEnn=zeros(N,lN);
LL=zeros(N,lN);
C00=C0/1e16;
for j=1:N
for i=1:lN-1
    dE(j,1)=frac(j,1).*(E0/e);
    if dE(j,i)>5
      dE(j,i+1)=dE(j,i).*frac(j,i);
      
    else
      dE(j,i+1)=0;
    end
end
end
for j=1:N
for i=1:lN

if dE(j,i)>5
% xnl=1;%nonlocal stopping
%local stopping
xnl(j,i)=0.16*dE(j,i).^0.15;%nonlocal stopping
      xloc(j,i)=1-xnl(j,i);
      dEnl(j,i)=NN*ll.*semod(j,i).*(xnl(j,i)+xloc(j,i).*(1+pmax./amod(j,i)).*exp(-pmax./amod(j,i)));
    dEloc(j,i)=(xloc(j,i).*semod(j,i)./(2*pi*amod(j,i).^2)).*exp(-p(j,i)./amod(j,i));%Robinson model
    snl(j,i)=(log(1+1.1383*ln11(j,i))./(2*(ln11(j,i)+0.01321*ln11(j,i).^0.21226+0.19593*ln11(j,i).^0.5)));%analytical approximation to the numerical solution of the nuclear stopping power
    sn(j,i)=(84.62.*(Z1.*Z2Ga.*M1.*snl(j,i).*C00./((M1+M2Ga).*(Z1^0.66+Z2Ga^0.66))+Z1.*Z2As.*M1.*snl(j,i).*C00./((M1+M2As).*(Z1.^0.66+Z2As.^0.66))));%nuclear stopping power
    dEnn(j,i)=NN*ll.*sn(j,i);%energy due to nuclear stopping
    LL(j,i)=(dE(j,i)-dEloc(j,i)-dEnn(j,i))./(NN*semod(j,i).*(xnl(j,i)+xloc(j,i).*(1+pmax./amod(j,i)).*exp(-pmax./amod(j,i)))+NN*sn(j,i));%path length of scattering
else

      LL(j,i)=0;
end
end
end
lltotalRp=zeros(N,1);
lltotalRpx=zeros(N,1);
lltotalRpy=zeros(N,1);
for i=1:lN-1
    ww(:,1)=cos(phi(:,1).*gammamod(:,1));
    ww(:,i+1)=cos(phi(:,i+1)).*ww(:,i)-(1-ww(:,i+1).^2).^0.5.*(1-ww(:,i)).^0.5.*cos(phii(:,i));
    u(:,1)=abs((1-ww(:,1).^2)).^0.5.*cos(phii(:,1))+0.001;
    v(:,1)=abs((1-ww(:,1).^2)).^0.5.*sin(phii(:,1))+0.001;
    u(:,i+1)=cos(phi(:,i+1)).*u(:,i)+sin(phi(:,i+1)).*(u(:,i).*ww(:,i).*cos(phii(:,i))-v(:,i).*sin(phii(:,i)))./abs((1-ww(:,i).^2).^0.5);
    v(:,i+1)=cos(phi(:,i+1)).*v(:,i)+sin(phi(:,i+1)).*(v(:,i).*ww(:,i).*cos(phii(:,i))+u(:,i).*sin(phii(:,i)))./abs((1-ww(:,i).^2).^0.5);
end
for i=1:N
lltotalRp(i)=sum(LL(i,:).*real(ww(i,:)./((u(i,:).^2+v(i,:).^2+ww(i,:).^2).^0.5)));
lltotalRpx(i)=sum(LL(i,:).*real(u(i,:)./((u(i,:).^2+v(i,:).^2+ww(i,:).^2).^0.5)));
lltotalRpy(i)=sum(LL(i,:).*real(v(i,:)./((u(i,:).^2+v(i,:).^2+ww(i,:).^2).^0.5)));
end
rowsToDelete=lltotalRp>0;
rowsToDelete1=(lltotalRp>1000);
rowsToDelete2=(lltotalRp<1000)&(lltotalRp>0);
RRPP=lltotalRp(rowsToDelete);
RRPP1=lltotalRp(rowsToDelete1);
RRPP2=lltotalRp(rowsToDelete2);
[ff11,xx11]=hist(lltotalRp,N1);
[ff12,xx12]=hist(lltotalRpx,N1);
[ff13,xx13]=hist(lltotalRpy,N1);
rowsToDelete11=(xx11>3500);
rowsToDelete12=(xx11<=3500);
XX=[lltotalRp,lltotalRpy];
[fff21,xx22]=hist3(XX,'ctrs',{-2000:100:8000 -5000:100:5000});
aaa=-log10(C0/1e14)*1.15e-4+8.0e-4;%fitting parameter of exp(-x)
xx111=xx11-100;
xx1111=xx11./10000;
zz11origin=C0*ff11.*1e8./(trapz(xx11,ff11));
zz11origin1=zz11origin(rowsToDelete12);%MC profile shallow region
zz11=C0*ff11.*1e8./(trapz(xx11,ff11));%MC profile ion concentration
% zz12=C01*ff12.*1e8/((N*trapz(xx12,ff12)./size(xx12,2)));
% zz13=C01*ff13.*1e8/(trapz(xx13,ff13));
xmod=[xx11(rowsToDelete12),xx11(rowsToDelete11)];
% conmod=[zz11origin1,tailcon];
xx221=xx22{1,1};
xx222=xx22{1,2};
% zz3d=C01*fff21*5e11./trapz(xx222,trapz(xx221,fff21));
%mesh(xx221,xx222,zz3d);
figure(1);
plot(xx11,zz11origin,'b-',x*10,x1,'g.');%,x2*10,x21,'r-',x*10,x1,'g.');
xlim([0 12000])
ylim([1E16 1E21])
set(gca, 'YScale', 'log')
xlabel('depth(A)');
%ylabel('lateral depth');
ylabel('concentration(ion/cm^3)');
%scatter3(lltotalRpx,lltotalRpy,lltotalRp);
hold on
% fid1=fopen('D:\thz1\measuresample\ion implantation\ion implant range\SitoGaAs200keV1e15Rpmoliere-morepoints-density.txt','w');
% AA = [xx11;ff11./trapz(xx11,ff11)];
% fprintf(fid1,'%6f %6f\n',AA);
% fclose(fid1);
% hold on
% fid2=fopen('D:\thz1\measuresample\SIMS\20200603\SitoGaAs200keV1e15Rpmolieredechannelingprofile.txt','w');
% AA11 = [xmod*0.1;zz11origin];
% fprintf(fid2,'%6f %6f\n',AA11);
% fclose(fid2);

%fid2=fopen('D:\thz1\measuresample\ion implantation\SIMS\20200518\1-1-1.txt','w');
%AA1 = [x*1;x1*1];
%fprintf(fid2,'%6f %6f\n',AA1);
%fclose(fid2);