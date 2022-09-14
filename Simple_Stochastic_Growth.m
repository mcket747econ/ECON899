% PROGRAM NAME: Simple_Stochastic_Growth.M
% This program generates the value function and decision rules for
% a stochastic growth model. I wanted to get more loops involved but that
% ended up messing up the way the future value function was created so I
% had to separate out the two states of the world (I think my issue was
% ultimately linear algebra but I am not sure, I've done similar problems
% where everything worked out fine). A mod function may have made the loop
% work but I don't know
% Date: 9/13/22

% PARAMETERS
b=.99; %discount factor 
d=0.025; %depreciation rate
a=.36; %capital share (theta in the prompt) 

% ASSET VECTOR
klb=0.01; %lower bound of grid points
inc=0.025; %increments
kub=45;%upper bound of grid points
k=[klb:inc:kub];% asset (row) vector
N=size(k);
N=N(1,2);
c1=ones(N,1); % column vector of ones
K=c1*k; % rows are k, cols are k'

% SHOCKS STUFF
Q = [0.977, .023; .074, .926]; %transition matrix

Z = [1.25; .2]; %technology states

% TABULATE CURRENT RETURN (UTILITY) FUNCTION
cs=zeros(N,N,2); %rows are k, cols are k', needs 2 states
for s = 1:2
    cs(:,:,s)=max(Z(s)*K'.^a-(K-K'*(1-d)),0);%cons, no negatives
    us=log(cs);
    t=isinf(us);
    j=find(t==1);
    us(j)=-realmax;
end
%So I tried to get it work in a loop over states below but that didn't turn
%out well. So I split up my utility and consumption into two different
%matrices, one per state

us_good = us(:,:,1);
us_bad = us(:,:,2);
visr_good = us_good(:,1)';
visr_bad = us_bad(:,1)';



% TABULATE INITIAL VALUE FUNCTION GUESS
%visr=us(:,:,2)'; %choose utility associated with k'=0; initial value function
%%
tic %timing in here because that's when the real business starts 
pcntol=1; %tolerance for value function iteration
n=1; %if want to run vfi for a set number of iterations
while pcntol >.0001

   vis_good = c1*visr_good; %generates future value function matrix from above row vector; will include the transition matrix later one
   vis_bad = c1*visr_bad;


   %CONSTRUCT TOTAL RETURN FUNCTION
   wis_good = us_good+b*(vis_good*Q(1,1)+vis_bad*Q(1,2)); 
   wis_bad = us_bad+b*(vis_bad*Q(2,2)+vis_good*Q(2,1)); 

   %CHOOSE HIGHEST VALUE (ASSOCIATED WITH k' CHOICE)
   [vsr_good,I_good]=max(wis_good'); %since max gives highest element in each column of a matrix, I is the decision rule while vsr is the placeholder value function
   [vsr_bad,I_bad]=max(wis_bad');
   
   tol=max([abs(vsr_good-visr_good), abs(vsr_bad-visr_bad)]); %use sup norm for tolerance
   pcntol=tol/abs(vsr_bad(1,N));
   visr_good = vsr_good;%update value functions
   visr_bad = vsr_bad;
   n=n+1; 
end
save 899_PS1vdr vsr_bad vsr_good I_bad I_good k;
save 899_PS1parm b a 

toc %because who cares how fast the graphs take
%%
%plot(k,k([I])-k)% plot change in the decision rule
%% First figure
figure(1);
hold on;
plot(k, vsr_good, 'r'); 
plot(k, vsr_bad, 'b'); 
xlabel('k_t');
ylabel('V(k,z)'); 
legend({'Z=1.25', 'Z=.2'}, 'Location',"best");
title('Value Functions');
hold off;

%% Second Figure
figure(2);
hold on
plot(k, k([I_good]), 'r');
plot(k,k([I_bad]),'b')
xlabel('k_t');
legend({'Z=1.25', 'Z=.2'}, 'Location', "best")

title('Decision rule');
hold off;

%% Third figure
figure(3);
hold on;
plot(k, k([I_good])-k, 'r');
plot(k,k([I_bad])-k,'b')
xlabel('k_t');
legend({'Z=1.25', 'Z=.2'}, 'Location', "best")

title('Decision rule changes');
hold off;

