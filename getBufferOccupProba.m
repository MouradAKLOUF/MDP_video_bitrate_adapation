
function Pq = getBufferOccupProba(q,b,action,intrb,intrq)

if length(q)==1 , Pq=1; return; end  

L=length(q);
q1=q(end);
q0=q(end-1);
a=action(end-1);
b0=b(end-1);

q1_indx=ceil(q1/intrq);
q0_indx=ceil(q0/intrq);
b0_indx=ceil(b0/intrb);  

c=0;

for k=L:-1:2 
    
    qq1=ceil(q(k)/intrq);
    qq0=ceil(q(k-1)/intrq);
    bb0=ceil(b(k-1)/intrb);  
    
    if qq1== q1_indx && qq0== q0_indx && bb0==b0_indx && a==action(k-1) 
        c=c+1; 
    end    

end

Pq=c/(L-1);
end