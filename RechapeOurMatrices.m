function [T,Rew] = RechapeOurMatrices(P,R)

S=size(P,1);
A=size(P,3);

    for s = 1:S
        for a=1:A
            for sp=1:S
                T(s,a,sp) = P(s,sp,a);
            end
        end
    end
sum(P(2,:,2))
Rew = zeros(S,A);

for a=1:A
   Rew(:,a) = sum(P(:,:,a).*R(:,:,a),2);
end;

end

