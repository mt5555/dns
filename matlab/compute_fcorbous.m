Lz=input('Aspect ratio Lz: ');
Fr=input('Fr: ');
Ro=input('Ro: ');
epsf = input('energy input rate eps_f: ');
k_f = input('forcing shell number k_f: ');

Lz_short = min(Lz,1)

k_w = k_f*2*pi/Lz_short;
U = (epsf/k_w)^(1/3);

if Lz <= 1
  H = 2*pi/k_w; %used to be 1/k_w for the older flows
  L = H/Lz;
end

if Lz > 1
  L = 2*pi/k_w
  H = L*Lz 
end

bous = U/H/Fr;
fcor = U/L/Ro;

disp(sprintf('bous = %f', bous));
disp(sprintf('fcor = %f', fcor));
disp(sprintf('delt = %f', 0.2*pi/(bous+fcor)));
