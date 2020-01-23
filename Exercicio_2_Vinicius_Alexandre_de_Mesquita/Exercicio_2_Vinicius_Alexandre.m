%   UNIVERSIDADE FEDERAL DO CEARÁ 
%   PROGRAMA DE PÓS-GRADUAÇÃO EM ENGENHARIA ELÉTRICA
%   CONTROLE ADAPTATIVO PREDITIVO
%   Vinícius Alexandre de Mesquita

clear;
close all;
Ts = 1;        % Tempo de amostragem
N = 3;        % Horizonte de predição
Nu = 3;       % Fator
d = 0;        % Atraso
A = [1 -0.8];     
B = [zeros(1,d) 0.4 0.6];

lambda = 0.8;     % Lambda (Ponderação)
lbd0 = 0;            % Índice do vetor lambda que será zerado
alfa = 0.8;

%Diferentes C
Ce(1,:) = [1 0 0];
Ce(2,:) = [1 -alfa 0];
Ce(3,:) = conv((Ce(2,1:2)),(Ce(2,1:2)));
Ce(4,:) = conv(A,Ce(2,1:2));
for i=1:4
    
    [kw,kf,ke,G,F,E,At,Bt,D,H] = calcula_K_ss_Vinicius(A,B,Ce(i,:),d,N,Nu,lambda,lbd0);
    
    sim('Exercicio_2_simulink_Vinicius',40)
    
    du_ant(:,1) = sld.signals.values(:,1);
    y_ant(:,1) = sld.signals.values(:,2);
    tempo = sld.time;
    y = sld.signals.values(:,4);
    du = sld.signals.values(:,3);
    u = sld.signals.values(:,5);
    if i==1
        y_novo = y;
        du_novo = du;
    end
    
    subplot(2,2,i)
    stairs(tempo,y,'Linewidth',1.5)
    hold on
    stairs(tempo,u,'Linewidth',1.5);
    stairs(tempo,du,'Linewidth',1.5);
    title({...
    ['c_1 = ' num2str(Ce(i,1)) '   c_2=' num2str(Ce(i,2)) '   c_3=' num2str(i,3)]}...
    ,'FontSize',14)
    legend('y','u','\Delta u','Location','east','FontSize',12)
    hold off
    
end


figure()
stairs(tempo(1:15,1),y_ant(1:15,1),'r--O','Linewidth',2)
hold on
stairs(tempo(1:15,1),y_novo(1:15,1),'b-','Linewidth',1)
stairs(tempo(1:15,1),du_ant(1:15,1),'m--O','Linewidth',2)
stairs(tempo(1:15,1),du_novo(1:15,1),'g-','Linewidth',1)
legend('y_{anterior}','y_{novo}','du_{anterior}','du_{novo}',...
    'Location','east','FontSize',14)
title('Comparação com o exercício anterior','FontSize',14)
