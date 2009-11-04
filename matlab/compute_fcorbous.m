Lz=input('Aspect ratio Lz: ');
Fr=input('Fr: ');
Ro=input('Ro: ');
epsf = input('energy input rate eps_f: ');
k_f = input('forcing shell number k_f: ');

k_f = k_f*2*pi/Lz;
U = (epsf/k_f)^(1/3);
H = 1/k_f;
L = H/Lz;

fcor = U/L/Ro;
bous = U/H/Fr;

disp(sprintf('fcor = %f', fcor));
disp(sprintf('bous = %f', bous));
