clc
clear

%% add paths
addpath(genpath("./PQN"));
addpath(genpath(".../gurobi910/linux64/matlab/")); % add in gruobi folder

%% Input parameters
k=80;
h=10;
pho = sqrt(nInstances);
gn = 10;
for i =1:gn
    Group{i}=(i-1)*100+1:i*100;
end

%% Load in data
load ./Xnd.txt
load ./yn.txt
load ./wd.txt

%% PQN method
X = Xnd;
[nInstances, nVars] = size(X);
y = yn;
w = wd;
[rid, cid] = find(wd);
utrue = w;
utrue(rid,1) = 1.;

% Initial guess of parameters
uSimplex = ones(nVars,1)*(1/nVars);%zeros(nVars,1);

% Set up Objective Function
funObj = @(w)L0Obj(w,X,y,pho);

% Set up Simplex Projection Function
funProj = @(w)ProjCSimplexGL_Gruobi(w,k,Group,h);

% Solve with PQN
options.maxIter = 50;
options.SPGiters = 5;
tStart = cputime;
[uout, err, Timespent, obj] = minConF_PQN_V3(funObj,uSimplex,funProj,options);
tEnd = cputime - tStart

[uout, zout] = ProjCSimplexGL_Gruobi_v2(uout, k, Group,h);
[B, Ranktmp] = sort(-uout);
Rank = sort(Ranktmp(1:k));
uout(Ranktmp(1:k));

%% Acc score
Indtrue = find(utrue);
C = intersect(Rank,Indtrue);
Acc = length(C)/k


