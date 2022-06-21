
function [v,policy] = policy_iterr(Poly,T,R,gamma,theta)

  if (nargin<5)
      theta = 0.001;
  end

  [Ns,Na] = size(Poly);
  policy=zeros(Ns,1);
  v = zeros(Ns,1);
  va = zeros(Na,1);
  pi = Poly;
    
  stable=false;

 while (stable==false)
  
  initValue=v;
  v = policy_eval(pi,T,R,gamma,theta,initValue);
  b=pi;
      for s = 1:Ns
          for a=1:Na
              va(a) = 0;
              for sp=1:Ns
                  va(a) = va(a) + T(s,a,sp)*(R(sp,a) + gamma*v(sp));
              end
          end
          [v(s),idx(s)]=max(va);
      end
      pi = zeros(Ns,Na);  
      for s = 1:Ns
          pi(s,idx(s))=1;
      end

   if isequal(b,pi) 
     stable=true;
   end
 
 end

 for i=1:Ns
 policy(i,1)=find(pi(i,:));
 end
 
end
