%%
clc; clear; close all;

%% define inputs

Qmax=30;  	% max level of buffer
Qlow=ceil(Qmax*16/100);	% 2 thresholds of the buffer
Qhigh=ceil(Qmax*83/100);

SegSize=2;    			%in sec
MovieSize=600;			%in sec
NbrFrames=MovieSize/SegSize; 	

b=zeros(1,NbrFrames);         % trace of estimated throughput
%%
b(1,1:50)=2000;
b(1,40:60)=4000; %%%%
b(1,61:100)=2000;
b(1,90:110)=500; %%%%
b(1,111:140)=2000;
b(1,141:160)=500; %%%%
b(1,161:190)=2000;
b(1,191:210)=4000; %%%%
b(1,211:240)=2000;
b(1,241:260)=500; %%%%
b(1,261:300)=2000;

%%
% b(1,1:100)=2000;
% b(1,101:150)=3000; %%%%
% b(1,151:200)=2500;
% b(1,201:230)=1500; %%%%
% b(1,231:300)=2000;

%%

Bmax=max(b);            % max measured throughput    
C=6;                    % segmenting bandwidth into C intervals ( i used it also with buffer occup )
intrb=Bmax/C;			% lenth of bandwidth segment --> we need it in proba calculation , please see line 123 and other functions


%Q=[500 1500 2500];		%all the available representation --> a(U(t)) --> v(t+1) 
Q=[300 700 1500 2500 3500];	
Act=size(Q,2);

a=0.4;         % try 035 035 015 015 then 04 04 01 01
bb=0.4;  
c=0.1;  
d=0.1;  
eps=1.5;	% smal cst so the rewards wont be zero

%%

intrq=Qmax/C;			% lenth of buffer uccup segment --> we need it in proba calculation
max_q_prime=(max(b)/min(Q))-1;              % max of buffer changing rate (q prime in the paper ))
min_q_prime=(min(b)/min(Q))-1;              % min of buffer changing rate --> should be negative
intr_qprime=(max_q_prime-min_q_prime)/C		% segment lenth  --> we need it in proba calculation

%% initialization

P= zeros(C,C,Act,'double'); 	% proba transition matrix from U(t) to U(t+1)
R= zeros(C,C,Act,'double');	% rewards matrix from U(t) to U(t+1)

P(:,:,:)=1/C;   	
R(:,:,:)=1;     	

discount=0.9; 		% discount factor
%%

V= [];          % video rate vector
Lambda= [];     % swithing rate
N= 5;           % how much you look in the past to find Lambda 
                % (sometimes you don't need to look to the beginning of the streaming N <=NbrFrames )

Vpd=zeros(1,NbrFrames,'double');  % vedio rate of the PD regulator 

Pchn=zeros(C,C,'double');   	% channel model : proba transition matrix from b(t) to b(t+1)
Pchn(:,:)=1/C;              	% equi-probable in the start
		
action=[];       %empty vector : action at state U(t) --> should be 1, 2 or 3 ( for 500, 1000 & 2500 resp)
q=[];            %empty vector : buffer occup
q_prime=[];      %empty vector : video changing rate
wait=[];         %empty vector : waiting time when buffer is full
t=[];
%% we can call it the "Main"
t(1)=1;
q(1)=0;
q_prime(1)=0;
V(1)=0;
Vpd(1)=0;
Lambda(1)=1;
wait(1)=0;
k=1;
                
Policy0 = (1/Act)*ones(C,Act);% initiation policy ; used in policy-iterr-algo to calculated the appropiate policy

NextState=ceil(b(1,1)/intrb);

while(k<NbrFrames) 

ThisState=NextState;                                 % the anticipated future state , is now the present state.                                       
%[Value, policy] = mdp_policy_iteration(P, R, discount);    

[T,Rew] = RechapeOurMatrices(P,R);                   % finding the right decision

% if q(k)>Qhigh || q(k)<Qlow
%     Rew0= zeros(C,C,Act,'double');
% end

[v,policy] = policy_iterr(Policy0,T,Rew,discount);
%
action(k)= policy(ThisState);

newV=Q(action(k));                                   % V(k+1)

delta_T=(newV/b(k))*SegSize;

if (q(k)+SegSize-(newV/b(k))*SegSize-wait(k))<=0
    q(k+1)=q(k)+SegSize;
    q_prime(k+1)=(b(k)/newV);
else
    q(k+1)=q(k)+SegSize-(newV/b(k))*SegSize-wait(k);
    q_prime(k+1)=(b(k)/newV)-1;
end

if q(k+1) > Qmax  
    wait(k+1)=(q(k+1)+SegSize-(newV/b(k))*SegSize)-Qmax;  %wait
else
    wait(k+1)=0;
end 

V = [V , newV];      % upgrading V, Lambda & Vpd

if length(V)<= N,
    f=2; L=length(V); % f (starting inx) and L (vector size) used in lines 113-***
else
    f=length(V)-N+1; L=N; end   
                    
if sum(V(f:end))/L == V(end),
    Lambda(k+1)=1; 
else
    Lambda(k+1)=0;
end


Kd=SegSize;		% PD regulator parameters ( have no idea how to get them )
Kp=1;
v0=b(k);
if q(k+1)>= Qhigh,
    q0=Qhigh;
elseif q(k+1)<= Qlow, 
    q0=Qlow;
else
    q0=(Qlow+Qhigh)/2;
end

% Kd=1;		% PD regulator parameters ( have no idea how to get them )
% Kp=1;
% v0=500;
% q0=Qlow;

if q(k+1)>Qhigh || q(k+1)<Qlow
    Vpd(k+1)= v0+(q(k+1)/SegSize)*(Kp*(q(k+1)-q0)+Kd*q_prime(k)); 
else
    Vpd(k+1)=V(end-1);
end

% finding the 4 new rewords

if length(V)<= N,
    f=2; L=length(V);
else
    f=length(V)-N+1; L=N-1; end    % need f & L in next line

m=sum(V(f:end))/(L);
sigma=(sqrt(sum(V(f:end)-m).^2)/(L));
r1=log(m+eps); %rmk

if (Lambda(k+1) ==0), r2=log(m+eps); else r2=0; end
if q(k+1)>Qhigh
    r3= -log(abs(newV-Vpd(k+1))+eps); %rfk
    r4= Qhigh-q(k+1);
elseif q(k+1)<Qlow
    r3= -log(abs(newV-Vpd(k+1))+eps); %rsighma
    r4= q(k+1)-Qlow;
else
    r3= 0;
    r4= -abs(q(k+1)-(Qlow+Qhigh)/2);
end
%
NextState=ceil(b(1,k+1)/intrb);
R(ThisState,NextState,action(k))= a*r1+bb*r2+c*r3+d*r4; 	%upgrading the reword to next state 

%		 % upgrading proba transition matrix of the channel
i= ceil(b(k)/intrb);           % find at wich bandwith interval you are; i=j=4000/200 
j= ceil(b(k+1)/intrb);

%%
temp=(C*Pchn(i,j)+1)/(C+1); 	% i'm using formula (15) from the paper , temp is just a temporary variable
Pchn(i,:)=(C*Pchn(i,:))/(C+1);
Pchn(i,j)=temp;
P(:,:,action(k))=Pchn; %upgrading the proba transition to next state
%%
% temp=(C*Pchn(i,j)+1)/(C+1); 	% i'm using formula (15) from the paper , temp is just a temporary variable
% Pchn(i,:)=(C*Pchn(i,j))/(C+1);
% Pchn(i,j)=temp;
% 
% Pq = getBufferOccupProba(q(2:end),b(2:k+1),action,intrb,intrq);	% Proba [q(t+1)/q(t),b(t),a(u(t))]	
% Pq1 = getBufferChangRateProba(q_prime(2:end),b(2:k+1),action,intrb,min_q_prime,intr_qprime);  % Proba [q_prime(t+1)/b(t),a(u(t))] 
% Pv = getVedioRateVectorProba(V(2:end),action,N);	% Proba [V(t+1)/V(t),a(u(t))] and V is vedio rate vector ( its vector not scalar)
% Plmd = getLamndaProba(Lambda(2:end),action);		% Proba [Lambda(t+1)/Lambda(t),(u(t))]
% 
% temp=Pq.*Pq1.*Pv.*Plmd.*Pchn(i,j);	%upgrading the proba transition to next state
% P(ThisState,:,action(k))=P(ThisState,:,action(k))*(1-temp)/(sum(P(ThisState,:,action(k)))-P(ThisState,NextState,action(k))); %the sum of proba should be 1 so we normalize all the other transition proba 
% P(ThisState,NextState,action(k))=temp;

%%
t(k)=k;
k=k+1; %go to next state pls;
end

while q(k)>0
    q(k+1)=q(k)-(V(end)/b(end))*SegSize;    % continue playing , not sure if this is correct
    k=k+1;
end

Waiting=ceil(wait./SegSize);   % waiting in nbr of segments

%% plese open and see the results
t;
V;
Vpd;
wait;
action;

%% figures 
maxx = max(q);minn= min(q);scal= max(V);scal2= min(V(2:end));
qplot=((scal*(q+minn))/maxx)+scal2;

figure;
subplot(221),plot(q,'r')
title('buff occup'); xlabel('nbr of segmt');ylabel('buff occup ( sec )');

%figure;
subplot(222),plot(Vpd,':b'),hold on, plot(b,'r');
title('Vpd'); xlabel('nbr of segmt');ylabel('Vpd( kbps ) '); legend('Vpd','bdw evolution');
axis([-1 301 -1 max(Q)+500])

%figure;
subplot(223), plot(qplot,'r'); hold on, plot(V(2:end),'b');
title('buff occup (scaled) & video bitrate');xlabel('nbr of segmt'); legend('buff occup','video bitrate');
axis([-1 301 -1 max(Q)+500])

%figure;
subplot(224), plot(b,'r'); hold on, plot(V(2:end),'b');
title('bndw & video bitrate');xlabel('nbr of segmt'); legend('bdw evolution','video bitrate');
axis([-1 301 -1 max(Q)+500])



