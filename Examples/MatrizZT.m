%% Construccion de Matriz de impedancia de retorno por tierra
function [ZT]=MatrizZT(numlineas,p_suelo,Mxy,s,Miu)
LT=zeros(numlineas,numlineas);       
for k=1:numlineas                                                        
    for m=1:numlineas                                                      
         if k==m  
                LT(k,m)=LT(k,m)+log(1+(p_suelo/Mxy(k,2))); % se añaden los elementos en la diagonal trifásica                                                              
         else
                if k<m % Para los elementos fuera de la diagonal
                Dprima=sqrt((Mxy(k,1)-Mxy(m,1))^2+(Mxy(k,2)+Mxy(m,2)+2*p_suelo)^2); % Se calcula Dik_p
                D=sqrt((Mxy(k,1)-Mxy(m,1))^2+(Mxy(k,2)+Mxy(m,2))^2); % Se calcula Dik
                LT(k,m)=LT(k,m)+log(Dprima/D); % Se añaden los elementos fuera de la diagonal.
                LT(m,k)=LT(k,m);
                end
         end
    end
end
ZT=((s*Miu)/(2*pi))*LT;
end