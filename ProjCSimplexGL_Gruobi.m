function [up] = ProjCSimplexGL_Gruobi(u,k,Group,h)

[n,m] = size(u);
g = length(Group);

H1 = eye(n);
f1 = -u;

H2 = zeros(g);
f2 = zeros(g,1);

H = blkdiag(H1, H2);
f = [f1; f2];

gcn = 0;
for i = 1 : g
    gcn = gcn + length(Group{i});
end

A1 = ones(1,n);
A2 = zeros(1,g);
for i = 1:g
    for j = 1:length(Group{i})
        A1t = zeros(1,n);
        A1t(Group{i}(j)) = 1;
        A1 = [A1; A1t];
        A2t = zeros(1,g);
        A2t(i) = -1;
        A2 = [A2; A2t];
    end
end
A1 = [A1; zeros(1,n)];
A2 = [A2; ones(1,g)];
A = [A1 A2];

b = zeros(gcn+2,1);
b(1) = k;
b(end) = h;

lb = zeros(n+g,1);
ub = ones(n+g,1);

%options = optimoptions('quadprog','Display','off');
%upt = quadprog(H,f,A,b,[],[],lb,ub,[],options);
model.A = sparse(A);
model.rhs = b;
model.ub = ub;
model.lb = lb;
model.obj = 2*f;
model.Q = sparse(H);
params.outputflag = 0;
output = gurobi(model,params);

up = output.x(1:n);
zp = output.x(n+1:end);


