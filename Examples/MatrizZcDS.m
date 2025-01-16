%% Matriz de impedancia del conductor o interna por Aproximación de Deri-Semlyen
function [Zc]=MatrizZcDS(numlineas,numconducfase,r,Ro,s,Miu,t)
Zc=zeros(numlineas,numlineas);
for k=1:numlineas
    for m=1:numlineas 
        if k==m 
            if k<=numconducfase % Para las fases:
                Rcd=Ro/(pi*r(k)^2);
                ZAF=t*(sqrt(s*Miu*Ro)/(2*pi*r(k)));
                Zc(k,m)=Zc(k,m)+sqrt(Rcd^2+ZAF^2)*(1/numconducfase); % se añaden los elementos en la diagonal.
            else % Para los hilos de guarda
                Rcd=Ro/(pi*r(k)^2);
                ZAF=t*(sqrt(s*Miu*Ro)/(2*pi*r(k)));
                Zc(k,m)=Zc(k,m)+sqrt(Rcd^2+ZAF^2); % se añaden los elementos en la diagonal.
            end   
        end
    end
end
end