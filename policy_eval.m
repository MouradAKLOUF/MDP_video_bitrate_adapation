function v = policy_eval(pi,T,R,gamma,theta,initValue)

if (nargin<5)
    theta = 0.001;
end

[Ns,Na] = size(pi);
v = initValue;
Delta = inf;

while (Delta>=theta)
    Delta = 0;
    v_old = v;
    for s = 1:Ns
        v(s) = 0;
        for a=1:Na
            for sp=1:Ns
                v(s) = v(s) + pi(s,a)*T(s,a,sp)*(R(sp,a) + gamma*v_old(sp));
            end
        end
        Delta = max(Delta,abs(v_old(s)-v(s)));
    end
end
