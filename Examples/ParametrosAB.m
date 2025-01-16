function [A,B,Zabc,Yabc]=ParametrosAB(s,numfases,numlineas,numconducfase,Miu,SigmaSuelo,SigmaAl,SigmaAire,Epsilon,Ro,t,longitud,R,Mxy,w1)
 
%% Matriz de capacitancia 
P=MatrizP(numlineas,numconducfase,R,Mxy); % Se calcula la matriz P
C=2*pi*Epsilon*inv(P); % Se calcula la matriz C
A=zeros(numfases,numfases,length(s));
B=A;
for ciclos=1:length(s)
    
s1=s(ciclos);

 %% Matriz de impedancia externa geométrica (inductancia externa)

LG=(Miu/(2*pi))*P;
ZG=s1*LG; % Se calcula la matriz ZG

 %% Matriz de impedancia de retorno por tierra (método de las imágenes)

% Para frecuencias regulares
gammap=sqrt(s1*Miu*SigmaSuelo); % Constante de propagacion
p_suelo=1/gammap; % Profundidad de penetración suelo (imágenes)
ZT=MatrizZT(numlineas,p_suelo,Mxy,s1,Miu); % se obtiene la impedancia por retorno a tierra


 %% Impedancia propia o interna (Zc)

p_conduc=1/sqrt(s1*Miu*SigmaAl); % Profundidad de penetración conductor 
gammac=sqrt(s1*Miu*(SigmaAl+s1*Epsilon)); % Constante de propagación
r=Mxy(:,3); % Radios de los conductores
Zc=MatrizZcDS(numlineas,numconducfase,r,Ro,s1,Miu,t); % Aproximación por Deri-Semlyen

%% Matriz de admitancia y de impedancia total

G=(SigmaAire/Epsilon)*C; % Matriz de conductancia
Y=(G+s1*C); % Se construye Y
Ide=(LG*C)/(Miu*Epsilon); % Se comprueba el resultado mediante la identidad
Z=(Zc+ZG+ZT); % Se obtiene la impedancia total de la línea

%% Reducción de hilos de guarda o neutros

Yff=Y(1:numfases,1:numfases); % Se construye Y de fases

Zff=Z(1:numfases,1:numfases); % se genera Z de fases
Zgg=Z(numfases+1:end,numfases+1:end); % Se obtiene Z de guardas
Zfg=Z(1:numfases,numfases+1:end); % Se obtiene Z de fase-guardas
Zgf=Z(numfases+1:end,1:numfases); % se obtiene Z de guardas-fase
Zhg=Zfg*inv(Zgg)*Zgf; % Se genera la matriz de reducción de guardas
Zeq=Zff-Zhg; % Se hace la reducción de guardas

%% Análisis Modal

AM=Zeq*Yff; % Se construye la matriz A
[M,Lambda]=eig(AM); % Se obtienen los modos de A
V0=diag(Lambda); % Se extraen los eigenvalores solamente
gamma=sqrt(V0); % Se evalúa la cosntante de propagación (modos propagación)
% Vi=w./imag(gamma); % se obtiene las velocidades de propagación

%% Modelo de la linea en circuito PI

Psi=M*sqrt(Lambda)*inv(M); % Matriz modal                                    
Yc=Zeq\Psi; % Admitancia caracteristica
Z0=Psi\Zeq; % Impedancia caracteristica
X=gamma*longitud;

Yabc(:,:,ciclos)=Zeq;
Zabc(:,:,ciclos)=Yff;

% Yabc(:,:,ciclos)=Yc;
% Zabc(:,:,ciclos)=inv(Yc);

% A(ciclos,:,:)=Yc*coth(Psi*longitud); % Matriz A del equivalente PI
% B(ciclos,:,:)=-Yc*csch(Psi*longitud); % Matriz B del equivalente PI

A(:,:,ciclos)=Yc*cotanh(M,X,inv(M)); % Matriz A del equivalente PI
B(:,:,ciclos)=-Yc*cosech(M,X,inv(M)); % Matriz B del equivalente PI

end
end