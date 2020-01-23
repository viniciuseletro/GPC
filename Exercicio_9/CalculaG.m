function [AA,BB,DD,G,F,E,H] = CalculaG(N,Cd,C)
% ny = nº de saídas y
% nu = nº de entradas u
[ny,nu] = size(Cd);

%% Confere a existência de C;
if ~exist('C')
    cnex = 1;
else
    cnex = 0;
end

%% Calcula o Atil de cada saída
for j=1:ny
    aa=1;
    for i=1:nu
        aa = conv(aa,Cd(j,i).den{:});
    end
    Atil(j,1:length(aa)+1) = conv(aa,[1 -1]);
    if cnex ==1
        C(j,1:length(aa)) = aa;
    end
end

%% Cálcula o número de estados
[a,nest]=size(Atil);
nest = nest-1;

%% Calcula e monta a matriz Bbarra
p = 1;
Bbarra = zeros(ny*nu,length(Atil));
for k=1:ny
    for i=1:nu
        bb=1;
        for j=1:nu
            if j~=i
                bb=conv(bb,Cd(k,j).den{:});
            else
                bb=conv(bb,Cd(k,j).num{:});
            end
        end
        Bbarra(p,1:length(bb))=bb;
        p=p+1;
    end
end

Bbarra=Bbarra(:,2:end)';

bg = reshape(Bbarra,[nest,nu,ny]);

BB = zeros(nest*ny,nu);

j=1;
for i=1:nest:nest*ny
    BB(i:nest+i-1,:) = bg(:,:,j);
    j=j+1;
end

%% Monta a matriz Abarra
AA = zeros(nest*ny,nest*nu);

aa = 0;
i = 1;
for j=1:nest:nest*ny
    AA(j:nest+j-1,j:nest+j-1) = [-Atil(i,2:end)' [eye(nest-1);zeros(1,nest-1)]];
    i = i+1;
end

%% Monta a matriz D

if cnex==0
    C = repmat(C,ny,1);
end
C = [C,zeros(ny,(nest+1-length(C)))];
DD = zeros(nest*ny,ny);
d = C-Atil;
d=d(:,2:end);
i=1;
for j=1:nest:nest*ny
    DD(j:nest+j-1,i) = d(i,:)';
    i = i+1;
end

%% Cálculo H
H = zeros(ny,ny*nest);
j=1;
for i=1:nest:nest*ny
    H(j,i) = 1;
    j=j+1;
end

%% Cálculo G,F e E
G = zeros(N*ny,N*nu);
% F = zeros(N*ny,nest*nu);
% E = zeros(N*ny,ny);

k=0;
j=1;
for i=1:ny:(ny*N)
    if i==1
        gg = H*(AA^k)*BB;
        G(i:i+ny-1,k*nu+1:k*nu+nu) = gg;
        ff = H*(AA^(k+1));
        F(1:ny,:) = ff;
        ee = H*(AA^(k))*DD;
        E = ee;
    else
        gg = [H*(AA^k)*BB gg];
        G(i:i+ny-1,1:length(gg)) = gg;
        ff = H*(AA^(k+1));
        F = [F;ff];
        ee = H*(AA^(k))*DD;
        E = [E;ee];
    end
    
    k=k+1;
end


end

