%   UNIVERSIDADE FEDERAL DO CEARÁ 
%   PROGRAMA DE PÓS-GRADUAÇÃO EM ENGENHARIA ELÉTRICA
%   CONTROLE ADAPTATIVO PREDITIVO
%   Vinícius Alexandre de Mesquita

clear,close all, clc;
Ts = 0.03;
N = 3;
Nu = 2;
lambda = 0.05;     % Lambda (Ponderação)

Ct11 = tf(1,[0.7 1]);
Ct12 = tf(5,[0.3 1]);
Ct21 = tf(1,[0.5 1]);
Ct22 = tf(2,[0.4 1]);

Ct = [Ct11,Ct12;Ct21,Ct22];
Cd = c2d(Ct,Ts);
Cd.TimeUnit = 'minutes';

for ii=1:3
[ncd1,dcd1] = tfdata(Cd(1,1),'v');
[ncd2,dcd2] = tfdata(Cd(1,2),'v');
Atil1 = conv(conv(dcd1,dcd2),[1 -1]);

[ncd3,dcd3] = tfdata(Cd(2,1),'v');
[ncd4,dcd4] = tfdata(Cd(2,2),'v');
Atil2 = conv(conv(dcd3,dcd4),[1 -1]);

B1 = conv(ncd1,dcd2);
B1 = [B1,zeros(1,(length(Atil1)-length(B1)))];
B2 = conv(ncd2,dcd1);
B2 = [B2,zeros(1,(length(Atil1)-length(B2)))];
B3 = conv(ncd3,dcd4);
B3 = [B3,zeros(1,(length(Atil2)-length(B3)))];
B4 = conv(ncd4,dcd3);
B4 = [B4,zeros(1,(length(Atil2)-length(B4)))];

%Monta matriz Dbarra
Ce1(1,:) = [1 0 0 0];
Ce2(1,:) = [1 0 0 0];
Ce1(2,:) = conv(conv(dcd1,dcd2),[1 -0.8]);
Ce2(2,:) = conv(conv(dcd3,dcd4),[1 -0.8]);
Ce1(3,:) = conv([1 -0.5 0],[1 -0.5]);
Ce2(3,:) = conv([1 -0.5 0],[1 -0.5]);
D1= Ce1(ii,2:end)'-Atil1(2:end)';
D2= Ce2(ii,2:end)'-Atil2(2:end)';

%Monta Abarra
[mat,nat] = size(Atil1);
Ab1 = [-Atil1(2:end)' [eye(nat-2);zeros(1,nat-2)]];
Ab2 = [-Atil2(2:end)' [eye(nat-2);zeros(1,nat-2)]];
Bb1 = [B1(2:end)' B2(2:end)'];
Bb2 = [B3(2:end)' B4(2:end)'];

%Monta A
AA = [Ab1,zeros(length(Ab1),length(Ab2))...
    ;zeros(length(Ab2),length(Ab1)),Ab2];

% Monta B
BB=[Bb1;Bb2];

% Monta D
[md,nd]=size(D1);
DD = [D1,zeros(md,nd);zeros(md,nd),D2];

nest = length(Cd);
na = length(AA);


%Cálculo G

H1 = [1,zeros(1,length(Ab1)-1)];
H = [H1,zeros(1,length(Ab1));zeros(1,length(Ab1)),H1];

G = zeros(N*nest,N*nest);
F = zeros(N*nest,na);
E = zeros(N,1);

for i=1:N
    if i==1
        gg = H*BB;
        G(1:2,1:2) = gg;
        [mgg,ngg]=size(gg);
        ff = H*(AA^i);
        F = ff;
        ee = H*(AA^(i-1))*DD;
        E = ee;
    else
        gg = [H*(AA^(i-1))*BB gg];
        [mgg,ngg]=size(gg);
        G((i-1)*mgg+1:mgg*i,1:ngg) = gg;
        ff = H*(AA^i);
        [mff,nff] = size(ff);
        F = [F;ff];
        ee = H*(AA^(i-1))*DD;
        E = [E;ee];
    end
end

Gnu = G(:,1:nest*Nu);
Q = lambda*eye(Nu*nest);

K = ((Gnu'*Gnu)+Q)\Gnu';
K = K(1:nest,:);


k1 = K(1:2,1:2);
k2 = K(1:2,3:4);
k3 = K(1:2,5:6);

KR = (k1+k2+k3);

KF = K*F;

KE = K*E;

sim('Exercicio_4_2_Vinicius',3)

[yer,u] = simul.signals.values;
tempo = simul.time;
y1 = yer(:,1);
y2 = yer(:,2);
ref1 = yer(:,3);
ref2 = yer(:,4);
u1 = u(:,1);
u2 = u(:,2);

figure(ii)

subplot(2,1,1)
plot(tempo,y1,'LineWidth',1.5)
hold on
plot(tempo,y2,'LineWidth',1.5)
plot(tempo,ref1,'k--','LineWidth',0.1)
plot(tempo,ref2,'m--','LineWidth',0.1)
ylabel('Outputs and References')
legend('y_1','y_2','Referência 1','Referência 2',...
    'Location','southeast')
hold off
    title({...
    ['Incerteza na saída da planta = 1.2'],[...
    'c_{11} = ' num2str(Ce1(ii,2)) ...
    '   c_{12}= ' num2str(Ce1(ii,3)) ...
    '   c_{13}= ' num2str(Ce1(ii,4))]...
    ['c_{21} = ' num2str(Ce2(ii,2)) ...
    '   c_{22}= ' num2str(Ce2(ii,3)) ...
    '   c_{23}= ' num2str(Ce2(ii,4))]}...
    ,'FontSize',10)

subplot(2,1,2)
plot(tempo,u1,'LineWidth',1.5)
hold on
plot(tempo,u2,'LineWidth',1.5)
ylabel('Inputs')
legend('u_1','u_2')

end
