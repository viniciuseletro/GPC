function [h] = tanque(qi,hant,Ts,prmt_tanque)
if hant<=0;
    hant=0;
end
h = hant + (Ts/prmt_tanque.A)*qi -...
    (Ts*prmt_tanque.k1/prmt_tanque.A)*sqrt(hant);
end

