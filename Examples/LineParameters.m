function [Zrto,Yr]=LineParameters(s)
% Determine Frequency Depend parameters of overhead-Transmision Line
% Z=Zt+Zg+Zc and Y=Yg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s: Frequency or Laplace variable
% Nb: Number of Sub Current Carrying (Haz)
% NNhg: Number of Guard Cables 0 if there is no Guard Cables
% r: Sub Current Carrying Radius [meters] First phases, second Guard Cables
% R: Carrying Current haz Radius [meters] 1 if there is only one conductor by Haz
% Uo: Vacuum permeability [H/m]
% Eo: Medium permitivity [F/m]
% p: Ground resistivity [ohm.m]
% zd: Height vector [meters], Contains N+NNhg Datas
% yd: Horizontal Position vector , Contains N+NNhg Datas
% re: Phase Resistivity [ohm.m]
% rg: Guard Cables Resistivity [ohm.m]
% tr: Ripple factor 1 if it is not considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  nf = 3; ng = 2;
%         x = [-6   0   6   -3.932  3.932];
%         y = [15  24  15  30  30];
%         rad_fase=0.03105;
%         cond1 = [rad_fase  0.000  3.206E-8  1  1  rad_fase;
%                 rad_fase  0.000  3.206E-8  1  1  rad_fase;
%                 rad_fase  0.000  3.206E-8  1  1  rad_fase;
%                 0.0056    0.000  3.206E-8    1  1  0.0056;
%                 0.0056    0.000  3.206E-8    1  1  0.0056];
%         murt = 1.0;
%         rhot = 100.0; 
%         tr=1;
% zd=y;
% yd=x;
% p=rhot;
% NNhg=ng;
% r=cond1(:,1)'; R=cond1(:,6); Nb=cond1(:,5); re=cond1(:,3);
% rai=0;                              % Medium resistivity
% Ns=length(s);                       % Number of Frequencies
% Uo=4*pi*1e-7;                       % Vacuum permeability [T.m/A]
% Eo=8.85e-12;                        % Medium permeability [T.m/A]
% N=length(zd);                       % Number of conductors (phases + guard cables)
% Nf=nf;                              % Number of phases


%         nf = 3; ng = 2;
%         x = [-5.49   0   5.49   -3.05  3.05];
%         y = [14.4  15.62  14.4  18.21  18.21];
%         rad_fase=(3.09E-2)/2;
%         cond1 = [rad_fase  0.000  5.87E-8  1  1  rad_fase;
%                 rad_fase  0.000  5.87E-8  1  1  rad_fase;
%                 rad_fase  0.000  5.87E-8  1  1  rad_fase;
%                 (1.26E-2)/2    0.000  2.1E-8    1  1  (1.26E-2)/2;
%                 (1.26E-2)/2    0.000  2.1E-8    1  1  (1.26E-2)/2];
%         murt = 1.0;
%         rhot = 100.0; 
%         tr=1;


        nf = 3; ng = 2;
        x = [-12.8016   0   12.8016   -8.9916  8.9916];
        y = [20.7265  20.7265  20.7265  32.9185  32.9185];
        rad_fase=(3.556E-2)/2;
        cond1 = [rad_fase  0.000  0.0431E-8  1  1  rad_fase;
                rad_fase  0.000  0.0431E-8  1  1  rad_fase;
                rad_fase  0.000  0.0431E-8  1  1  rad_fase;
                (1.27E-2)/2    0.000  3.1069E-8    1  1  (1.27E-2)/2;
                (1.27E-2)/2    0.000  3.1069-8    1  1  (1.27E-2)/2];
        murt = 1.0;
        rhot = 100.0; 
        tr=3.5;


zd=y;
yd=x;
p=rhot;
NNhg=ng;
r=cond1(:,1)'; R=cond1(:,6); Nb=cond1(:,5); re=cond1(:,3);
rai=8.85E-12;                              % Medium resistivity
Ns=length(s);                       % Number of Frequencies
Uo=4*pi*1e-7;                       % Vacuum permeability [T.m/A]
Eo=8.85e-12;                        % Medium permeability [T.m/A]
N=length(zd);                       % Number of conductors (phases + guard cables)
Nf=nf;                              % Number of phases

% Equivalent Radius for lines with haz conductors
for j=1:N
    if Nb(j) ~= 1
        Req(j)=nthroot((Nb(j)*r(j)*R(j)^(Nb(j)-1)),Nb(j));        % Haz Equivalent Radius
    else
        Req(j)=r(j);                                              % Radius Without Haz
    end
end
Req=[Req r((Nf+1):N)];

% Distance Between Current carrier and image Current carrier
for i=1:N
    for j=1:N
        if i==j
            D(i,j)=2*zd(i);
        else
            D(i,j)=sqrt((zd(i)+zd(j))^2+(yd(j)-yd(i))^2);
        end
    end
end

% Distance Between Current carrier and Current carrier
for i=1:N
    for j=1:N
        if i==j
            d(i,j)=Req(i);
        else
            d(i,j)=sqrt((zd(i)-zd(j))^2+(yd(j)-yd(i))^2);
        end
    end
end

% Distance Between Current carrier and image Current carrier with complex penetration depth. This value depends on the frequency
Dc=zeros(N,N,Ns);
% Inverse Complex Penetration Depth neglecting Et and Ut aprox Uo
sigma=sqrt(s.*Uo/p);
for i=1:N
    for j=1:N
        if i==j
            Dc(i,j,:)=2*zd(i)+2./sigma;
        else
            Dc(i,j,:)=sqrt(((zd(i)+zd(j))*ones(1,Ns)+2./sigma).^2+ones(1,Ns).*(yd(j)-yd(i)).^2);
        end
    end
end

Rp=log(D./d);        % Potential Maxwell Matrix
Lg=Uo.*Rp/(2*pi);    % Geometric Inducatance Henry/meteres
C=2*pi*Eo*inv(Rp);   % Geometric Capacitance F/meters
G=rai*C/Eo;          % Geometric Conductance Siemens/meters
for i=1:N
    Rcd(i)=re(i)/(pi*(r(i))^2);                                     % DC Current resistance
    Zaf(i,:)=(tr/(2*pi*r(i))).*sqrt(s.*Uo*re(i));                   % High Frequency Impedance considering the sama radius in all phases
    zc(i,:)=sqrt(ones(1,Ns)*((Rcd(i))^2)+(Zaf(i,:)).^2)/Nb(i);      % Current carrier Impedance
    Zc(i,i,:)=zc(i,:);                                              % Conductor Impedance
end

for i=1:Ns
    Y(:,:,i)=G+s(i)*C;
    Zg(:,:,i)=s(i)*Lg;
    Zt(:,:,i)=(s(i)*Uo/(2*pi))*log(Dc(:,:,i)./D);
end

% Conductor Impedance
Ztot=Zg+Zc+Zt;
% Kron's Reduction for guard-cables
if NNhg~=0
    Yr=Y(1:Nf,1:Nf,:);
    G=G(1:Nf,1:Nf,:);
    for i=1:Ns
        Zrto(:,:,i)=Ztot(1:Nf,1:Nf,i)-Ztot(1:Nf,Nf+1:N,i)*inv(Ztot(Nf+1:N,Nf+1:N,i))*Ztot(Nf+1:N,1:Nf,i);
        Zgc(:,:,i)=Zg(1:Nf,1:Nf,i)-Zg(1:Nf,Nf+1:N,i)*inv(Zg(Nf+1:N,Nf+1:N,i))*Zg(Nf+1:N,1:Nf,i);
        Zcc(:,:,i)=Zc(1:Nf,1:Nf,i)-Zc(1:Nf,Nf+1:N,i)*inv(Zc(Nf+1:N,Nf+1:N,i))*Zc(Nf+1:N,1:Nf,i);
        Ztc(:,:,i)=Zt(1:Nf,1:Nf,i)-Zt(1:Nf,Nf+1:N,i)*inv(Zt(Nf+1:N,Nf+1:N,i))*Zt(Nf+1:N,1:Nf,i);
    end
else
    Zrto=Ztot;Yr=Y;
    Zgc=Zg;
end

end