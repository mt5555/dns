function [Dlnorm,Dtnorm]=spher_harm_weight(name,ext,xx,ndir_use)
%
% weight the input structure function in each of the directions with the spherical harmonic function.
%
%
%


% x,y,z directions
  dir(1,:)=[1,0,0]
dir(2,:)=[0,1,0]
dir(3,:)=[0,0,1]

% face diagonals
dir(4,:)=[1,1,0]
dir(5,:)=[1,0,1]
dir(6,:)=[0,1,1]
dir(7,:)=[-1,1,0]
dir(8,:)=[-1,0,1]
dir(9,:)=[0,-1,1]

% body diagonals
dir(10,:)=[1,1,1]
dir(11,:)=[-1,1,1]
dir(12,:)=[1,-1,1]
dir(13,:)=[1,1,-1]


% face 1,2 directions
dir(14,:)=[0,1,2]
dir(15,:)=[0,2,1]
dir(16,:)=[0,-1,2]
dir(17,:)=[0,-2,1]

dir(18,:)=[1,2,0]
dir(19,:)=[1,0,2]
dir(20,:)=[-1,2,0]
dir(21,:)=[-1,0,2]

dir(22,:)=[2,1,0]
dir(23,:)=[2,0,1]
dir(24,:)=[-2,1,0]
dir(25,:)=[-2,0,1]

% body 1,1,2 directions
dir(26,:)=[1,1,2]
dir(27,:)=[1,2,1]
dir(28,:)=[2,1,1]

dir(29,:)=[-1,1,2]
dir(30,:)=[-1,2,1]
dir(31,:)=[-2,1,1]

dir(32,:)=[1,-1,2]
dir(33,:)=[1,-2,1]
dir(34,:)=[2,-1,1]

dir(35,:)=[1,1,-2]
dir(36,:)=[1,2,-1]
dir(37,:)=[2,1,-1]

% 2,2,1 directions
dir(38,:)=[1,2,2]
dir(39,:)=[2,1,2]
dir(40,:)=[2,2,1]

dir(41,:)=[-1,2,2]
dir(42,:)=[-2,1,2]
dir(43,:)=[-2,2,1]

dir(44,:)=[1,-2,2]
dir(45,:)=[2,-1,2]
dir(46,:)=[2,-2,1]

dir(47,:)=[1,2,-2]
dir(48,:)=[2,1,-2]
dir(49,:)=[2,2,-1]


% face 1,3 directions
dir(50,:)=[0,1,3]
dir(51,:)=[0,3,1]
dir(52,:)=[0,-1,3]
dir(53,:)=[0,-3,1]

dir(54,:)=[1,3,0]
dir(55,:)=[1,0,3]
dir(56,:)=[-1,3,0]
dir(57,:)=[-1,0,3]

dir(58,:)=[3,1,0]
dir(59,:)=[3,0,1]
dir(60,:)=[-3,1,0]
dir(61,:)=[-3,0,1]

% body 1,1,3 directions
dir(62,:)=[1,1,3]
dir(63,:)=[1,3,1]
dir(64,:)=[3,1,1]

dir(65,:)=[-1,1,3]
dir(66,:)=[-1,3,1]
dir(67,:)=[-3,1,1]

dir(68,:)=[1,-1,3]
dir(69,:)=[1,-3,1]
dir(70,:)=[3,-1,1]

dir(71,:)=[1,1,-3]
dir(72,:)=[1,3,-1]
dir(73,:)=[3,1,-1]

  ndir = length(dir);
ndir_vec=1:ndir;
%ndir_vec=[1,2,3, 4:3:ndir];

for i=ndir_vec
r = sqrt(sum(dir(i,:).^2));   % length of directional vector
rhat(1) = dir(i,1)/r;
rhat(2) = dir(i,2)/r;
rhat(3) = dir(i,3)/r;

if (sphere_harm == 2) then
% for j=2
%
% the five components m = -2, 1, 0, 1, 2 of the j=2 spherical harmonic
% projection of the mixed structure function $S_{ij}$, $i\neq j$
  
% the polar angle is t and the azimuthal angle  is p
  y2_2 = 2*rhat(1)*rhat(2)              !sint*sint*sin2p
  y22 = rhat(1)^2 - rhat(2)^2         !sint*sint*cos2p
  y2_1 = rhat(1)*rhat(3)                !sint*cost*cosp
  y21 = rhat(2)*rhat(3)                 !sint*cost*sinp
  y20 = (3*rhat(3)^2 - 1)              ! 3*cost*cost - 1  Dl(:,i,6)
  
  Dlnorm(:,i,1)=Dl(:,i,6)*y2_2
  Dlnorm(:,i,2)=Dl(:,i,6)*y22
  Dlnorm(:,i,3)=Dl(:,i,6)*y2_1
  Dlnorm(:,i,4)=Dl(:,i,6)*y21
  Dlnorm(:,i,5)=Dl(:,i,6)*y20


  Dtnorm(:,i,1)=Dt(:,i,6)*y2_2
  Dtnorm(:,i,2)=Dt(:,i,6)*y22
  Dtnorm(:,i,3)=Dt(:,i,6)*y2_1
  Dtnorm(:,i,4)=Dt(:,i,6)*y21
  Dtnorm(:,i,5)=Dt(:,i,6)*y20

end

end
