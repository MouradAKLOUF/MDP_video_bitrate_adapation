
function Pq1 = getBufferChangRateProba(q_prime,b,action,intrb,min_q_prime,intr_qprime)

if length(q_prime)==1 , Pq1=1; return; end  

q_prime=q_prime-min_q_prime;

L=length(q_prime);

q1prime=q_prime(end);
a=action(end-1);
b0=b(end-1);

q1prime_indx=ceil(q1prime/intr_qprime);
b0_indx=ceil(b0/intrb);  

c=0;

for k=L:-1:2 
    
    qq1prime=ceil(q_prime(k)/intr_qprime);
    bb0=ceil(b(k-1)/intrb);  
    
    if qq1prime== q1prime_indx && bb0==b0_indx && a==action(k-1) 
        c=c+1; 
    end    

end

Pq1=c/(L);
end