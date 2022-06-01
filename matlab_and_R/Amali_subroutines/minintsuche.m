function [Von,Bis]=minintsuche(x,ll)
% sucht in x das Intervall mit Laenge ll, das die kleinsten Werte enth?lt

laenge=length(x);
m=laenge-ll;
werte=zeros(m,1);
for j=1:m
    werte(j)=sum(x(j:j+ll));
end
[w,Von]=min(werte);
Bis=Von+ll;