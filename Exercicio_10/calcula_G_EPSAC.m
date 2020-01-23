function [Gnu] = calcula_G_EPSAC(N,Nu,delta,qu,Ts,prmt_tanque)

qi = qu*ones(1,N);
qd = qi+delta;
imp = prmt_tanque.h0;
impd = prmt_tanque.h0;

for i = 1:N
    imp(i+1) = tanque(qi(i),imp(i),Ts,prmt_tanque);
end

for i = 1:N
    impd(i+1) = tanque(qd(i),impd(i),Ts,prmt_tanque);
end

g = (impd-imp)/delta;

for i=1:N
    hi(i)= g(i+1)-g(i);
end

G = zeros(N,N);
for i = 1:N
    G(i:end,i) = hi(1:end-i+1);
end

Gnu = G(:,1:Nu);
Gnu(Nu:end,end) = g(2:N-Nu+2);


end

