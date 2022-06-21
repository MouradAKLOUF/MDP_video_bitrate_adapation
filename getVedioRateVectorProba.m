function Pv = getVedioRateVectorProba(V,action,N)

if length(V)<= N, Pv=1; return; end   

L=length(V);
lastVector=V(L-N+1:L);
a=action(end-1);
c=0;

for k=0:1:L-N 
   CurentVector=V(L-N-k+1:L-k);
if isequal(lastVector,CurentVector) && a==action(L-k), c=c+1; end    
end

Pv=c/(L-N+1);
end