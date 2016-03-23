function ep_new=new_extended_kalmaan(rcvd_sig,x,Qvn)
p(1)=1e-6;
ep(1)=.01;
% Qvn=1e-9;
N=length(rcvd_sig);
ph = i*2*pi/N;

for n=2:N
    z(n)=x(n).*exp(ph*(n-1)*ep(n-1));
    Hn=ph*(n-1)*exp(ph*(n-1)*ep(n-1)).*x(n);
    Kn = p(n-1)*conj(Hn)/(p(n-1)+Qvn);
    ep(n) = ep(n-1) + real(Kn*(rcvd_sig(n) -z(n)));
    p(n) = (1-Kn*Hn)*p(n-1);
end
    ep_new=ep(n);

end