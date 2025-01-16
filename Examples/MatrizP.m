%% Construccion de matriz de coeficientes de potencial de Maxwell
function [P]=MatrizP(numlineas,numconducfase,R,Mxy)
P=zeros(numlineas,numlineas);       
for k=1:numlineas                                                        
    for m=1:numlineas                                                      
        if k==m  
                if k<=numconducfase % Para los elementos en la diagonal de fases
                RMG=nthroot((numconducfase*(R)^(numconducfase-1)*Mxy(k,3)),numconducfase); % Se calcula el RMG.
                P(k,m)=P(k,m)+log(((2*Mxy(k,2)))/RMG); % se añaden los elementos en la diagonal trifásica
                else
                P(k,m)=P(k,m)+log(((2*Mxy(k,2)))/Mxy(k,3)); % se añaden los elementos en la diagonal de guardas.
                end                                                                 
        else
            if k<m % Para los elementos fuera de la diagonal
                d=sqrt((Mxy(k,1)-Mxy(m,1))^2+(Mxy(k,2)-Mxy(m,2))^2); % Se calcula Dik_p
                D=sqrt((Mxy(k,1)-Mxy(m,1))^2+(Mxy(k,2)+Mxy(m,2))^2); % Se calcula Dik
                P(k,m)=P(k,m)+log(D/d); % Se añaden los elementos fuera de la diagonal.
                P(m,k)=P(k,m);  
            end
        end
    end
end
end