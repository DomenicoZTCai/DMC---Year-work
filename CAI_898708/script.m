%% Q2
clear all;
load('CAI_898708_mkr.mat')
MFF=M(1:73,1:73);
MFC=M(1:73,74:81);
MCF=M(74:81,1:73);
MCC=M(74:81,74:81);

RFF=R(1:73,1:73);
RFC=R(1:73,74:81);
RCF=R(74:81,1:73);
RCC=R(74:81,74:81);

KFF=K(1:73,1:73);
KFC=K(1:73,74:81);
KCF=K(74:81,1:73);
KCC=K(74:81,74:81);

n=length(MFF);
A=[MFF zeros(n);zeros(n) MFF];
B=[RFF KFF;-MFF zeros(n)];

[modes, eigenvalues]=eig(MFF\KFF);
freq=sqrt(diag(eigenvalues))/2/pi;
freq_sorted=sort(freq)

%% Q3

% [modes, eigenvalues]=eig(-(A)\B);
% freq=imag(diag(eigenvalues))/2/pi;

% freq_posi=zeros(length(freq)/2,1);
% for i=1:length(freq)/2
%     freq_posi(i)=freq(i*2-1);
% end
% freq_posi=sort(freq_posi)

%% Q4

i=sqrt(-1);
C0=1;
Q0=zeros(n,1);
Q0(12)=C0;
vett_f=0:0.01:25;
for k=1:length(vett_f)
    ome=2*pi*vett_f(k);
    A=-ome^2*MFF+i*ome*RFF+KFF;
    x0=A\Q0;
    yA=x0(9);
    yB=x0(12);
    yAdd=-ome^2*yA;
    yBdd=-ome^2*yB;
    mod1(k)=abs(yAdd);
    fas1(k)=angle(yAdd);
    mod2(k)=abs(yBdd);
    fas2(k)=angle(yBdd);
end

figure
subplot 211;plot(vett_f,mod1);grid;title('yAdd/C');xlabel('[Hz]');ylabel('[m/s^2]');
subplot 212;plot(vett_f,fas1*180/pi);grid;xlabel('[Hz]');ylabel('[rad]');

figure
subplot 211;plot(vett_f,mod2);grid;title('yBdd/C');xlabel('[Hz]');ylabel('[m/s^2]');
subplot 212;plot(vett_f,fas2*180/pi);grid;xlabel('[Hz]');ylabel('[rad]');

%% Q5

i=sqrt(-1);
Q0=zeros(n,1);
Q0(3)=-5/2;Q0(4)=-25/12;Q0(6)=-5;Q0(9)=-5/2;Q0(10)=25/12;
vett_f=0:0.01:25;
for k=1:length(vett_f)
    ome=2*pi*vett_f(k);
    A=-ome^2*MFF+i*ome*RFF+KFF;
    x0=A\Q0;
    yA=x0(9);
    yAdd=-ome^2*yA;
    x0d=i*ome*x0;
    x0dd=-ome^2*x0;
    R0=MCF*x0dd+RCF*x0d+KCF*x0;
    R_VC=R0(6);
    mod1(k)=abs(yAdd);
    fas1(k)=angle(yAdd);
    mod2(k)=abs(R_VC);
    fas2(k)=angle(R_VC);
end

figure
subplot 211;plot(vett_f,mod1);grid;title('yAdd/p0');xlabel('[Hz]');ylabel('[m/s^2]');
subplot 212;plot(vett_f,fas1*180/pi);grid;xlabel('[Hz]');ylabel('[rad]');

figure
subplot 211;plot(vett_f,mod2);grid;title('R_VC/p0');xlabel('[Hz]');ylabel('[N]');
subplot 212;plot(vett_f,fas2*180/pi);grid;xlabel('[Hz]');ylabel('[rad]');

%% Q6
i=sqrt(-1);
Q0=zeros(n,1);
Q0(3)=-5/2;Q0(4)=-25/12;Q0(6)=-5;Q0(9)=-5/2;Q0(10)=25/12;
vett_f=0:0.01:25;
for k=1:length(vett_f)
    ome=2*pi*vett_f(k);
    A=-ome^2*MFF+i*ome*RFF+KFF;
    x0=A\Q0;
    yA=x0(9);
    yAdd=-ome^2*yA;
    x0d=i*ome*x0;
    x0dd=-ome^2*x0;
    R0=MCF*x0dd+RCF*x0d+KCF*x0;
    R_VC=R0(6);
    mod1(k)=abs(yAdd);
    fas1(k)=angle(yAdd);
    mod2(k)=abs(R_VC);
    fas2(k)=angle(R_VC);
end

MFF_new=modes'*MFF*modes;
RFF_new=modes'*RFF*modes;
KFF_new=modes'*KFF*modes;
Q0_new=modes'*Q0;

a=find(freq==freq_sorted(1));
b=find(freq==freq_sorted(2));
c=find(freq==freq_sorted(3));

q0=zeros(n,1);
for k=1:length(vett_f)
    ome=2*pi*vett_f(k);
    A_new=-ome^2*MFF_new+i*ome*RFF_new+KFF_new;
    q0_full=A_new\Q0_new;
    q0(a)=q0_full(a);q0(b)=q0_full(b);q0(c)=q0_full(c);
    x0_new=modes*q0;
    yA_new=x0_new(9);
    yAdd_new=-ome^2*yA_new;
    %
    x0d_new=i*ome*x0_new;
    x0dd_new=-ome^2*x0_new;
    R0_new=MCF*x0dd_new+RCF*x0d_new+KCF*x0_new;
    R_VC_new=R0_new(6);
    mod3(k)=abs(yAdd_new);
    fas3(k)=angle(yAdd_new);
    mod4(k)=abs(R_VC_new);
    fas4(k)=angle(R_VC_new);
end

figure
subplot 211;plot(vett_f,mod1,'linewidth',0.5);hold on;...
    plot(vett_f,mod3,'--','linewidth',2);...
    grid;title('yAdd/p0');xlabel('[Hz]');ylabel('[m/s^2]');...
    legend('full modes','first three modes');
subplot 212;plot(vett_f,fas1*180/pi,'linewidth',0.5);hold on;...
    plot(vett_f,fas3*180/pi,'--','linewidth',2);...
    grid;xlabel('[Hz]');ylabel('[rad]');...
    legend('full modes','first three modes');

figure
subplot 211;plot(vett_f,mod2,'linewidth',0.5);hold on;...
    plot(vett_f,mod4,'--','linewidth',2);...
    grid;title('R_VC/p0');xlabel('[Hz]');ylabel('[N]');...
    legend('full modes','first three modes');
subplot 212;plot(vett_f,fas2*180/pi,'linewidth',0.5);hold on;...
    plot(vett_f,fas4*180/pi,'--','linewidth',2);...
    grid;xlabel('[Hz]');ylabel('[rad]');...
    legend('full modes','first three modes');