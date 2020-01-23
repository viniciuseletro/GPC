function [Hcl,Ur,dU] = GPCgen(B,A,d,Ts,N,Nu,lambda,lbd0)
W=ones(N+d,1);

Atio=conv([1 -1],A);
E=zeros(N+d,N+d);
% EB=zeros(N,N+1);
EB=zeros(N,N+1+2*d);
F=zeros(N,length(A));
G=zeros(N+d,N+d);
Gp=zeros(N,d+1);
K=zeros(Nu,Nu);

for i=1:N+d
    [Ej,Fj]=deconv([1 zeros(1,length(Atio)-2+i)],Atio);
    E(i,:)=[Ej zeros(1,length(E)-length(Ej))];
    F(i,:)=Fj(i+1:end);
    EB(i,:)=conv(E(i,:),B);
    G(i,:)=[flip(EB(i,1:i)) zeros(1,(N+d)-length(EB(i,1:i)))];
    Gp(i,:)=EB(i,i+1:i+1+d);
end

Gnu=G(:,1:Nu);
Q=lambda*eye(Nu);
if lbd0>0
    for i=1:lbd0
        Q(i,i)=0;
    end
end
K=(Gnu'*Gnu+Q)\Gnu';
k=K(1,:);

kr=k*W;
[mf,nf]=size(F);
for i=1:nf
    kf(i)=k*F(:,i);
end

kgp=k*Gp
% for i=1:d+1
%     kgp(1,i)=k*Gp(:,i);
% end
% ku=[1 0]-kgp*[1 -1]
kf
kr

integ=filt(1,[1 -1],Ts);
% H=filt([zeros(1,d+1) B],A,Ts);
H=filt([0 B],A,Ts);
Kr=filt(kr,1,Ts);
Kf=filt(kf,1,Ts);
Kgp=filt([0 kgp],1,Ts);
Hcl1=series(feedback(1,Kgp),integ);
HC=series(H,Hcl1);
Hcl2=feedback(HC,Kf);
Hcl=series(Kr,Hcl2);

Ur1=feedback(Hcl1,series(H,Kf));
Ur=series(Kr,Ur1);
dU=series(Ur,filt([1 -1],1,Ts));

% step(Hcl,dU)
% hold on
% step(Ur)


end

