function [T,Sa] = designSpectra(Sds,Sd1)

% Design Response Spectra 11.4.5
T0 = 0.2*Sd1/Sds;
Ts = Sd1/Sds;
Tl = 12;

% Range of T values
T = (0:0.001:Tl)';
n= numel(T);
Sa = NaN(n,1);
for i = 1:n
    if T(i) <= T0
        Sa(i) = Sds*(0.4 + 0.6*(T(i)/T0));
    elseif T(i) > T0 && T(i) <= Ts
        Sa(i) = Sds;
    elseif T(i) > Ts && T(i) <= Tl
        Sa(i) = Sd1/T(i);
    else
        Sa(i) = Sd1*Tl/T(i)^2;
    end
end