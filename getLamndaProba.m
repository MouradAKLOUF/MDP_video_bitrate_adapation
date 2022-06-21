
function Plmd = getLamndaProba(Lamnda,action)

 % this function ( & the other ) calculate the conditional probabilties
 % which form the transition probabilty
 % empirically from the current data , for instance if u have a vect (t) =
 % [ a b b a c b a b a  b] and u want to find Pr(a/b), u need to count how 
 % many time appear (b then a) so in the exemple u have  3a/b  1b/b 1a/c .
 
if length(Lamnda)==1 , Plmd=1; return; end  
 
L=length(Lamnda);
l1=Lamnda(end);
l0=Lamnda(end-1);
a=action(end-1);
c=0;

for k=L:-1:2  
if l1== Lamnda(k) && l0== Lamnda(k-1)&& a==action(k-1) , c=c+1; end    
end

Plmd=c/(L-1);
end