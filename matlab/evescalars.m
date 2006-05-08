%
%########################################################################
%#  plot of DNS structure function output file
%########################################################################
%
%
pack;
clear all;
apply_fix=0;
mu_hyper_value =0;
range=0:50;
%range = 1
%range=0:.5:2;
fid2=-1;

%fid=endianopen('/home/taylorm/ccs/dns/src/rot3d/rot3d_sto0000.0000.scalars','r');
%nx=128;
%***********************************************************************************
%fid=endianopen('/research/elunasin/ascitransfer_sept05/2D_1024a=1CFL.09v.8e-40000.0000.scalars','r')

%fid=endianopen('/superscratch/emanalo/ascitransfer_nov05/2D_1024a50dxCFL.09v.8e-40000.0000.scalars','r');

%most current for DNS run 10/31/05 good no hook
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r2_CFL.09v.8e-40004.0597.scalars','r');
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_1024r0_fw100_bw50CFL.2v1e-4hypo10_0000.0000.scalars','r'); %tturn = .025244

%most current run for alpha 11/06/05
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r910a20.48CFL.09to.05v.8e-4leta_a0007.6861.scalars','r');


%this for testing tturn 1 and 2
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r910a20.48CFL.09to.05v.8e-4leta_a0007.6861.scalars','r');
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048L_etaCFL.09v.8e-40000.0000.scalars','r');


%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r910a10.24CFL.09v.8e-4l_eta_a0008.1734.scalars','r');
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r910a2.048CFL.05v.8e-4l_eta_a0007.3972.scalars','r');

%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r9a20.48CFL.09to.05v.8e-4leta_a0007.0048.scalars','r');
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r9a10.24CFL.09v.8e-4l_eta_a0007.5000.scalars','r');
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_2048r9a2.048CFL.05v.8e-4l_eta_a0006.7345.scalars','r');

%1024 resolution
mu_hyper_value =0;
%local computer
%fid = fopen('C:\matlabR12\work\2D_Ekman_f100_mu_1e2_1024r0CFL.05v1e-40000.0000.scalars','r');% no hook l_eps = .001, l_eta = .003
%fid = fopen('C:\matlabR12\work\2D_Ekman_f0_mu_.1_1024r1CFL.01v1e-40001.9519.scalars','r');% no hook l_eps = 0.001225, l_eta = 0.005251

%fid = fopen('C:\matlabR12\work\2D_1024_4f_a3.25CFL.05v1e-40000.0000.scalars','r');% hook l_eps =0.000890 , l_eta=0.002450
%fid = fopen('C:\matlabR12\work\2D_1024_25f_a30.25CFL.05v1e-40000.0000.scalars','r');% hook l_eps =0.000521 , l_eta=.000732

 

%NOTES:

%D O N T    F O R G E T     T O    E D I T   L I N E    212  and line 251, line 246, defnition for epsilon and ke_diss respectively     

%TO TOGGLE BETWEEN ALPHA FINITE AND ALPHA  INFTY CASE
%           MAY
%*******************************************************
fid = fopen('C:\matlabR12\work\2D_1024r3_nohypoCFL.01v1e-40005.0771.scalars','r'); %001989713, 003589 , tf = 8.5, tturn = 1.37
%fid = fopen('C:\matlabR12\work\2D_a3.25_ffval_355.30scaled_Helm_r0_1024v1e-40000.0000.scalars','r'); %

%fid = fopen('C:\matlabR12\work\2D_a100_ffval_355.30scaled_Helm_r0_1024v1e-40000.0000.scalars','r'); %
fid = fopen('C:\matlabR12\work\2D_a1000_ffval_355.30scaled_Helm_r0_1024v1e-40000.0000.scalars','r'); %
%fid = fopen('C:\matlabR12\work\2D_a1e6_ffval_355.30scaled_Helm_r0_1024v1e-40000.0000.scalars','r'); %
%fid = fopen('C:\matlabR12\work\2D_a_inf_ffval_355.30scaled_Helm_r0_1024v1e-40000.0000.scalars','r'); %
%**********************************
%%%% april 
%**********************************
%case alpha = 0
 % of course the same
 
%case alpha = 0 with decaying force
%fid = fopen('C:\matlabR12\work\2D_1024decay_a00008.2209.scalars','r'); %
  
%case alpha = 3.25
%fid = fopen('C:\matlabR12\work\2D_1024f_1e-5r4_nohypo_a3.25CFL.09v1e-40036.2648.scalars','r'); %
   
%case alpha = 100
%fid = fopen('C:\matlabR12\work\2D_1024f_.009r1_nohypo_a100CFL.02v1e-40002.4980.scalars','r'); %

%case alpha = 100 decaying case
%fid = fopen('C:\matlabR12\work\2D_a100_fdecay_inftycode_nohypo_1024CFL.01v1e-40001.2164.scalars','r'); %
%fid = fopen('C:\matlabR12\work\2D_a100r1_decay_inftycode_nohypo_1024CFL.01v1e-40002.3539.scalars','r'); %
%fid = fopen('C:\matlabR12\work\2D_a100r0_xfac_decay_inftycode_nohypo_1024CFL.01v1e-40001.0000.scalars','r'); %

%case alpha = 1000 forcing scaled 
%fid = fopen('C:\matlabR12\work\2D_1024f_.9537r0_nohypo_a1000CFL.02v1e-40000.0000.scalars','r'); %

%case alpha = 1e3 old code old forcing
%fid = fopen('C:\matlabR12\work\2D_1024r1_a_1e3_rerunCFL.01v1e-40001.2508.scalars','r'); %
%fid = fopen('C:\matlabR12\work\2D_1024r2_a_1e3_rerunCFL.01v1e-40002.4922.scalars','r'); %
%case alpha = 2000 unscaled forcing
%fid = fopen('C:\matlabR12\work\2D_1024r0_nohypo_a2000CFL.01v1e-40000.0000.scalars','r'); %

%case alpha = 2000 scaled forcing
%fid = fopen('C:\matlabR12\work\2D_1024r0_nohypo_f3.8_a2000CFL.01v1e-40000.0000.scalars','r'); %
%fid = fopen('C:\matlabR12\work\2D_1024r1_a_2e3_rerunCFL.01v1e-40001.2641.scalars','r'); %

%case  alpha = 1e4 scaled forcing
 %fid = fopen('C:\matlabR12\work\2D_1024r0_nohypo_f95.36_a1e4CFL.01v1e-40000.0000.scalars','r'); %
   
%case alpha = 1e6 unscaled forcing
%fid = fopen('C:\matlabR12\work\2D_1024r0_nohypo_a1e6CFL.01v1e-40000.0000.scalars','r'); %

%case alpha = infinity from old code
%fid = fopen('C:\matlabR12\work\2D_1024r0_a_inf_old_code_rerunCFL.01v1e-40000.0000.scalars','r'); %
%fid = fopen('C:\matlabR12\work\2D_1024r0_a_inf_old_code_rerunCFL.01v1e-40001.2330.scalars','r'); %

%case alpha = infinity from old code old forcing decay run
%fid = fopen('C:\matlabR12\work\2D_1024dec_a_inf_old_cde_reruCFL.01v1e-40002.3450.scalars','r'); %


%case alpha = infinity 
%fid = fopen('C:\matlabR12\work\2D_a_infinity_r0_xfac_decay_inftycode_nohypo_1024CFL.01v1e-40008.6078.scalars','r'); %


%************************************************************************
%           JAN-FEBRUARY-MARCH runs  do nt forget to update mu if necessary
%************************************************************************
%case: alpha =0; no hypo
    %fid = fopen('C:\matlabR12\work\2D_1024r0_nohypoCFL.01v1e-40000.0000.scalars','r'); % 0.002084,   0.003717
    %fid = fopen('C:\matlabR12\work\2D_1024r1_nohypoCFL.01v1e-40002.0597.scalars','r'); %l_eps = 0.002100, l_eta =0.003690
    %fid = fopen('C:\matlabR12\work\feb06\2D_1024r2_nohypoCFL.01v1e-40003.0000.scalars','r'); %l_eps =.002053 , l_eta =.003594
    %fid = fopen('C:\matlabR12\work\2D_1024r3_nohypoCFL.01v1e-40005.0771.scalars','r'); %001989713, 003589 , tf = 8.5, tturn = 1.37
   %**********************************************************
   %                                  with hypo alpha = 0
   %***************************************************************
  %fid = fopen('C:\matlabR12\work\feb06\2D_1024_10hypoCFL.05v1e-40000.0000.scalars','r'); %tf = 10
  % note I rerun this case below because there is  Re = 0 above the same as below so I think use Euu to solve this err
  %fid = fopen('C:\matlabR12\work\2D_1024r_10hypoCFL.02v1e-40000.0000.scalars','r'); %tf = 2.0
 
  %alpha=3.25  no hypo
    %fid = fopen('C:\matlabR12\work\2D_1024r0_nohypo_a3.25CFL.01v1e-40000.0000.scalars','r');%0.002264,0.004342
    %fid = fopen('C:\matlabR12\work\2D_1024r1_nohypo_a3.25CFL.01v1e-40001.1881.scalars','r');%0.002082 , 0.003571 
    %fid = fopen('C:\matlabR12\work\2D_1024r2_nohypo_a3.25CFL.01v1e-40002.4596.scalars','r');% 0.002107, 0.003669 
    %fid = fopen('C:\matlabR12\work\2D_1024r3_nohypo_a3.25CFL02v1e-40003.5827.scalars','r');
    %fid = fopen('C:\matlabR12\work\2D_1024r4_nohypo_a3.25CFL02v1e-40006.0829.scalars','r'); 0.001166801, 0.004070  
    %****************************************************************
    %                 forcing scaled to \alpha^2 with alpha used = alpha numerics
    %********************************************************
    %fid = fopen('C:\matlabR12\work\2D_1024f_1e-5r0_nohypo_a3.25CFL.02v1e-40000.0000.scalars','r');%
    %fid = fopen('C:\matlabR12\work\2D_1024f_1e-5r1_nohypo_a3.25CFL.09v1e-40005.0025.scalars','r');
   %fid = fopen('C:\matlabR12\work\2D_1024f_1e-5r4_nohypo_a3.25CFL.09v1e-40036.2648.scalars','r');
   
    %v = 1e-6
    %fid = fopen('C:\matlabR12\work\2D_1024f_1e-5r2_nohypo_a3.25CFL.01v1e-60060.0000.scalars','r');
   
    
    %********************************************************************************
    %                          alpha = 3.25 with hypo
    %*********************************************************************************
    %fid = fopen('C:\matlabR12\work\2D_1024r1_10hypo_a3.25CFL.05v1e-40006.0000.scalars','r');%

    
%alpha = 10 (9.765625e-003)
    %fid = fopen('C:\matlabR12\work\2D_1024r0_a10_nohypo_CFL.01v1e-40000.0000.scalars','r');% 0.002261, 0.004210
    %fid = fopen('C:\matlabR12\work\2D_1024r1_a10_nohypo_CFL.01v1e-40001.2016.scalars','r'); %0.002063,0.003680 
    %fid = fopen('C:\matlabR12\work\2D_1024r3_a10_nohypo_CFL.02v1e-40003.7386.scalars','r'); 

%alpha = 100
%fid = fopen('C:\matlabR12\work\2D_1024r0_nohypo_a100CFL.01v1e-40000.0000.scalars','r');%0.001607905, 0.001991
%fid = fopen('C:\matlabR12\work\2D_1024r1_nohypo_a100CFL.02v1e-40001.1334.scalars','r')%

%case alpha = 200 but with forcing dvided by alpha^2
%fid = fopen('C:\matlabR12\work\2D_1024f_26.2144r1_nohypo_a200CFL.01v1e-40001.3968.scalars','r')%

%case: alpha = 1000
    %fid = fopen('C:\matlabR12\work\feb06\2D_1024r3_nohypo_a1000CFL.01v1e-40003.7717.scalars','r'); %0.000554760,0.000431
    %note that the energy diss length is small.  In the case of alpha =\infty the energy is not edited to handle alpha = inf 
    %but in the finite alpha the energy is handled correctly. 
    %fid = fopen('C:\matlabR12\work\2D_1024r3_nohypo_a1000CFL.01v1e-40003.7717.scalars','r'); %

%case alpha infinity
     %fid = fopen('C:\matlabR12\work\2D_a_infty_nohypo_1024CFL.01v1e-40000.0000.scalars','r');%  0.000000small, 0.000476
      %fid = fopen('C:\matlabR12\work\2D_a_infty_nohypo_1024r1CFL.01v1e-40001.2399.scalars','r'); %0.small, 0.000443
      %fid = fopen('C:\matlabR12\work\2D_a_infty_nohypo_1024r3CFL.01v1e-40003.7179.scalars','r');
      %fid = fopen('C:\matlabR12\work\2D_a_infty_nohypo_1024r4CFL.01v1e-40004.9627.scalars','r');
      %fid = fopen('C:\matlabR12\work\2D_a_infty_nohypo_1024r5CFL.01v1e-40006.2030.scalars','r');
      %fid = fopen('C:\matlabR12\work\2D_a_infty_nohypo_1024r6CFL.01v1e-40007.4556.scalars','r');
  %***************************************************************************************************
  %                                            DECEMBER
 %***************************************************************************************************
 
 %fid = fopen('C:\matlabR12\work\2D_1024_10hypoCFL.05v1e-40000.0000.scalars','r');

%******************************************************************************************************************
%fid = endianopen('/superscratch/emanalo/ascitransfer_jan06/2D_Ekman_mu_1e2_1024r0CFL.05v1e-40000.0000.scalars','r');% no hook l_eps = .001, l_eta = .003
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_Ekman_mu_0_a1.63_1024r0CFL.05v1e-40000.0000.scalars','r');% hook but l_eps =.001 ., l_eta =.003 consistent with the old code with hypo 10 
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_1024_10hypoCFL.05v1e-40000.0000.scalars','r'); % l_eps = .001, l_eta = .003 old code with hypo
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_1024_10hypo_alpha1.63CFL.05v1e-40000.0000.scalars','r');% hook but l_eps =.001 ., l_eta =.003 
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_1024_nohypo_alpha1.63CFL.05v1e-40000.0000.scalars','r');% hook but l_eps =.002 ., l_eta =.003
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_1024_nohypo_alpha5.2CFL.05v1e-40000.0000.scalars','r');% hook but l_eps =.0019 ., l_eta =.003
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_Ekman_mu_0_1024r0CFL.05v1e-40000.0000.scalars','r');% no hypo hook but l_eps =.001874 ., l_eta =.015820 
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_Ekman_mu_10_1024r0CFL.08v1e-40000.0000.scalars','r');%l_eps =.001034 ., l_eta = .003971
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_Ekman_mu_.1_1024r0CFL.08v1e-40000.0000.scalars','r');%hook but l_eps =.001100 ., l_eta =.005076 
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_Ekman_1024r0CFL.08v1e-40000.0000.scalars','r');%l_eps = .001, l_eta = .003
%fid = endianopen('/superscratch/emanalo/ascitransfer_dec05/2D_Ekman_1024r0CFL.08v.8e-40000.0000.scalars','r');% but l_eps =.000885 ., l_eta =.003646 
%fid=endianopen('/superscratch/emanalo/ascitransfer_oct05/2D_1024F2to4CFL.09v.8e-40000.0000.scalars','r'); % l_eps = .001, l_eta = .003

%fid=endianopen('/superscratch/emanalo/ascitransfer_nov05/2D_1024CFL.05v1e-40000.0000.scalars','r');
%fid=endianopen('/superscratch/emanalo/ascitransfer_nov05/2D_1024k_a10CFL.05v1e-40000.0000.scalars','r');
%fid=endianopen('/superscratch/emanalo/ascitransfer_nov05/2D_1024k_a100CFL.05v1e-40000.0000.scalars','r');

%nx=2048;
nx = 1024;
%nx=256;
%nx=512;
%***********************************************************************************

nscalars=0;
nscalars_e=0;
while (1) 
  [ni,count]=fread(fid,1,'float64');
  if (count~=1) 
     if (fid2<0) 
        break; 
     end
     fclose(fid);
     fid=fid2;
     fid2=-1;
     [ni,count]=fread(fid,1,'float64');
     if (count~=1) 
       break;
     end
  end;
  nints=ni;
  ns=fread(fid,1,'float64');
  mu=fread(fid,1,'float64');
  alpha=fread(fid,1,'float64');
  
[nints,ns,nscalars];
  data1=fread(fid,[nints,ns],'float64');
  data2=fread(fid,[nints,ns],'float64');

  if (nscalars==0) 
    ints=data1;
    maxs=data2;
  else
    ints=[ints,data1];
    maxs=[maxs,data2];
  end
  nscalars=nscalars+ns;
  

  % now read the "expensive" integrals, which are not computed every time step
  ns_e = fread(fid,1,'float64');
  time = fread(fid,1,'float64');
  data1 = fread(fid,[ns_e,1],'float64');
  data1=[time;data1]
  
end

  disp(sprintf('nints=%i  total scalars read=%i',nints,nscalars))
    
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the scalars computed every time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=size(ints);
l=l(2);

if (ns>=11) 
   h_diss=ints(11,:);   
end
%****************************************************************************************
%                     WARNING !!!!! DONT FORGET TO EDIT THIS PART
%****************************************************************************************
  
  ke_diss_d=ints(10,:);       % THIS IS THE ALPHA = FINITE CASE  
  ke_diss_d_inf=ints(1,:);         % THIS IS THE ALPHA = INFINITY CASE 

ke_diss_f=ints(3,:);        % < u,F >   F=f, except for alpha model F=f'
vor_z=ints(4,:);
hel=ints(5,:);
ens_diss = -(ints(12,:) + mu_hyper_value * ints(7,:));      %added for computing the enstrophy diss-rate with Ekman damping as well as alpha
grad_ke = ints(2,:);
ke=ints(6,:);   %  ke 
ke_inf=ints(2,:);
ens = ints(7,:);   %  enstrophy
% ints(8,:)   < u,div(tau)' >  (alpha model only)
% ints(9,:)   < u,f>           (alpha model only)
% ints(1,:)  < u_xx,u_xx> >   (used for E_alpha dissapation term)

maxU=maxs(1,:);  % at time_after
maxV=maxs(2,:);  % at time_after
maxW=maxs(3,:);  % at time_after
%  maxs(4,:)     % max used for CFL, at time_after
maxvor=maxs(5,:);
time_after=maxs(6,:);
time=maxs(7,:);

Ea = ints(6,:) + .5*alpha^2 *ints(2,:); % at time
E_uu_inf = ke_inf;
E_uv_inf = ints(2,:);
E_vv_inf = ints(1,:);
time_2 = .5*(time(2:l)+time(1:l-1));
ke_diss_tot=(ke(2:l)-ke(1:l-1))./(time(2:l)-time(1:l-1));
ke_diss_tot_inf=(ke_inf(2:l)-ke_inf(1:l-1))./(time(2:l)-time(1:l-1));
ens_diss_tot = (ens(2:l)-ens(1:l-1))./(time(2:l)-time(1:l-1));
Ea_diss_tot=(Ea(2:l)-Ea(1:l-1))./(time(2:l)-time(1:l-1)); %(alpha-nergy)

ke_v=ke + alpha^4*ints(1,:)/2 + 2*alpha^2*ints(2,:)/2;

%USING THE 3 + 2(for alpha = infinity) DIFFERENT TYPES OF ENERGY
Euse_uu=ke; %E_uu
Euse_uv=Ea; %E_uv
Euse_vv=ke_v; %E_vv
Euse_uu_inf = E_uu_inf;
Euse_uv_inf = E_uv_inf;
Euse_vv_inf = E_vv_inf;

mean_ens = mean(ens);  %dont forget to change ke_diss_d above for alpha = \infty
epsilon_inf = ke_diss_d_inf; % this is the alpha = inifnity case
epsilon=-( ke_diss_d-mu*alpha^2*ints(1,:)-mu_hyper_value*ints(6,:) );  % this is for the alpha = f inite case
%epsilon = -ke_diss_d + 1e-4 * alpha^2 *ints(1,:);                %new epsilon                           %added for Ekman pumping
%epsilon = -ke_diss_d + ints(8,:);  %a_diss = ints(8,:)

v_rms = sqrt(2*Euse_uu/3);

%this formula from lambda from Daryls formula Hamiltons Prin... remove division by 3 on Euse_uu etc
lambda_uu=sqrt( mu*(2*Euse_uu) / (epsilon/15)  );
lambda_uv=sqrt( mu*(2*Euse_uv) / (epsilon/15)  );
lambda_vv=sqrt( mu*(2*Euse_vv) / (epsilon/15)  );
lambda_uu_inf=sqrt( mu*(Euse_uu_inf) / (epsilon_inf/15)  );
lambda_uv_inf=sqrt( mu*(Euse_uv_inf) / (epsilon_inf/15)  );
lambda_vv_inf=sqrt( mu*(Euse_vv_inf) / (epsilon_inf/15)  );
%Re = lambda^(-2);
Re_uu = 1*sqrt(2*Euse_uu) / mu;
Re_uv = 1*sqrt(2*Euse_uv) / mu;
Re_vv = 1*sqrt(2*Euse_vv) / mu;
Re_uu_inf = 1*sqrt(2*Euse_uu_inf) / mu;
Re_uv_inf = 1*sqrt(2*Euse_uv_inf) / mu;
Re_vv_inf = 1*sqrt(2*Euse_vv_inf) / mu;

%remove division by 3 on Euse_uu
R_l_uu = lambda_uu.*sqrt(2*Euse_uu/3)/mu;
R_l_uv = lambda_uv.*sqrt(2*Euse_uv/3)/mu;
R_l_vv = lambda_vv.*sqrt(2*Euse_vv/3)/mu;
R_l_uu_inf = lambda_uu_inf.*sqrt(2*Euse_uu_inf/3)/mu;
R_l_uv_inf = lambda_uv_inf.*sqrt(2*Euse_uv_inf/3)/mu;
R_l_vv_inf = lambda_vv_inf.*sqrt(2*Euse_vv_inf/3)/mu;

R_large = 1 * sqrt(2*Euse_uu) / mu;
eta = (mu^3 ./ abs(epsilon)).^(.25); %(dissipation lengthscale)
eta_ens = (mu^3 ./ abs(ens_diss)).^(1/6) * 2^(1/6); %(ens dissipation lengthscale)
%lambda_ke=sqrt( mu*(2*ke/3) / (epsilon_ke/15)  );
%R_l_ke=lambda_ke.*sqrt(2*ke/3)/mu;


disp(sprintf('max vor_z = %e',max(vor_z)));
disp(sprintf('mean_ens = %e',mean_ens));
disp(sprintf('Re_uu = %e',Re_uu(1:1)));
disp(sprintf('Re_uv = %e',Re_uv(1:1)));
disp(sprintf('Re_vv = %e',Re_vv(1:1)));
disp(sprintf('Re_uv_inf = %e',Re_uv_inf(1:1)));
disp(sprintf('Re_vv_inf = %e',Re_vv_inf(1:1)));
disp(sprintf('alpha = %e',alpha));
%disp(sprintf('Ea_diss_tot = %e',Ea_diss_tot));
%figure(5)
%clf
%hold on

%plot(time,mean(ke_diss_f+ke_diss_d),'g') % d/dt ke shld bal to zero
%plot(time,ke,'k')
%%%%plot(time_2,ke_diss_tot,'g') %eve note ke_diss_tot = ke_diss_f+ke_diss_d
%%%%plot(time,ke_diss_f,'b') %eve     <f.u>
 %%%%%%plot(time,ke_diss_d,'b') % eve   <- visc * enstrophy> from the eqn conserv energy
%plot(time,hel,'r')
%%%%%%title('KE: black,   \epsilon: blue  d(KE)/dt: green,    hel: red');
%title('KE: black,    d(KE)/dt: green,    hel: red');
%hold off
%xlabel('time')

figure(6)
 clf;
 plot(time,ens_diss,'m'); hold on; 
%semilogy(time,ke,'r'); hold on
%semilogy(time, Ea_diss_tot,'b');hold on
%%%%plot(time_2,-ke_diss_tot,'b.')
%semilogy(time, grad_ke, 'm'); hold on
%semilogy(time,ke_diss_d-1e-4*alpha^2*ints(1,:),'b')
%***********************************
%eve plot(time,-ke_diss_d ,'b')
%******************************
%title('KE: red,    ke-diss-d: blue    KE-uv: greeg       |grad u|^2:magenta');
%%%%hold off
%xlabel('time')
%print -dpsc ke.ps
%print -djpeg -r72 ke.jpg

figure(7)
clf;
%t = size(time);
%time = time(1:t(2)-1);
ens_diss_tot = [ens_diss_tot 0];
plot(time,ens_diss_tot,'b'); hold on;
plot(time,ens_diss,'m'); hold on; 
  %   semilogy(time,(1+4*pi*pi*alpha^2)^2.*ens, 'b' ); hold on;
%plot(time,maxvor,'r'); hold on;
%plot(.4188,2500,'o')
%plot(.4328,2500,'o')
%plot(.4603,2500,'o')
%%%%plot(.4894,2500,'*')
%plot(.5551,2500,'o')
%plot(.6034,2500,'o')
%plot(.8149,2500,'o')
%%%%plot(time,50000*ke,'k');
%hold off;
%axis([0,1,0,5000]);
%title('red: \omega(t) blue: \omega_\alpha(t)')
%xlabel('time')
%print -djpeg -r72 vor.jpg


%figure(2);  hold on;
%plot(time,R_l_uu,'b'); hold on;
%%title('R_\lambda');
%legend('R_{\lambda}', 'R_{\lambda}as a function of(total KE)')
%xlabel('time')
%print -djpeg -r72 rl.jpg

%figure(3);
%plot(time,lambda)
%     title('\lambda')
%     xlabel('time')
%print -djpeg -r72 lambda.jpg


figure(4); subplot(1,1,1)
plot(time,eta_ens* nx*pi*2*sqrt(2)/3 )
title('k_{nmax} \eta')
xlabel('time')
%print -djpeg -r72 kmaxeta.jpg

% averge d/dt energy to a number
ke_diss_tot = ke_diss_tot(length(ke_diss_tot)/2:length(ke_diss_tot));
ke_diss_tot = sum(ke_diss_tot)/length(ke_diss_tot);
disp(sprintf('ke_diss_tot from numerics (average over last half of data) = %f ',ke_diss_tot));
% averge d/dt energy to a number
ke_diss_tot_inf = ke_diss_tot_inf(length(ke_diss_tot_inf)/2:length(ke_diss_tot_inf));
ke_diss_tot_inf = sum(ke_diss_tot_inf)/length(ke_diss_tot_inf);
disp(sprintf('ke_diss_tot_inf from numerics (average over last half of data) = %f ',ke_diss_tot_inf));

Ea_diss_tot = Ea_diss_tot(length(Ea_diss_tot)/2:length(Ea_diss_tot));
Ea_diss_tot = sum(Ea_diss_tot)/length(Ea_diss_tot);
disp(sprintf('Ea_diss_tot from numerics (average over last half of data) = %f ',Ea_diss_tot));



%compare with the averaage of ke_diss_d 
ke_diss_d = ke_diss_d(length(ke_diss_d)/2:length(ke_diss_d));
ke_diss_d = sum(ke_diss_d)/length(ke_diss_d);
disp(sprintf('ke_diss_d from myscalars formula (average over last half of data) = %f ',ke_diss_d));

grad_ke_diss_d = -ints(2,:)*1e-4;
grad_ke_diss_d = grad_ke_diss_d(length(grad_ke_diss_d)/2:length(grad_ke_diss_d));
grad_ke_diss_d = sum(grad_ke_diss_d)/length(grad_ke_diss_d);
disp(sprintf('grad_ke_diss_d = ke_diss_d_inf from myscalars formula (average over last half of data) = %f ',grad_ke_diss_d));

Laplace_alpha_ke_diss_d = -ints(1,:)*1e-4*alpha^2;
Laplace_alpha_ke_diss_d = Laplace_alpha_ke_diss_d(length(Laplace_alpha_ke_diss_d)/2:length(Laplace_alpha_ke_diss_d));
Laplace_alpha_ke_diss_d = sum(Laplace_alpha_ke_diss_d)/length(Laplace_alpha_ke_diss_d);
disp(sprintf('Laplace_alpha_ke_diss_d from myscalars formula (average over last half of data) = %f ',Laplace_alpha_ke_diss_d));

Laplace_alpha_ke_diss_d_inf = -ints(1,:)*1e-4;
Laplace_alpha_ke_diss_d_inf = Laplace_alpha_ke_diss_d_inf(length(Laplace_alpha_ke_diss_d_inf)/2:length(Laplace_alpha_ke_diss_d));
Laplace_alpha_ke_diss_d = sum(Laplace_alpha_ke_diss_d_inf)/length(Laplace_alpha_ke_diss_d_inf);
disp(sprintf('Laplace_alpha_ke_diss_d_inf from myscalars formula (average over last half of data) = %f ',Laplace_alpha_ke_diss_d_inf));



% averge d/dt enstrophy to a number
ens_diss_tot = ens_diss_tot(length(ens_diss_tot)/2:length(ens_diss_tot));
ens_diss_tot = sum(ens_diss_tot)/length(ens_diss_tot);
disp(sprintf('ens_diss_tot from numerics (average over last half of data) = %f ',ens_diss_tot));

%compare with the averaage of ens_diss 
ens_diss = ens_diss(length(ens_diss)/2:length(ens_diss));
ens_diss = sum(ens_diss)/length(ens_diss);
disp(sprintf('ens_diss from myscalars formula (average over last half of data) = %f ',ens_diss));

% averge eta to a number
eta = eta(length(eta)/2:length(eta));
eta = sum(eta)/length(eta);
disp(sprintf('eta_energy - dissipation lengthscale (average over last half of data) = %1.9f ',eta));

% averge eta_ens to a number
eta_ens = eta_ens(length(eta_ens)/2:length(eta_ens));
eta_ens = sum(eta_ens)/length(eta_ens);
disp(sprintf('eta_ens - ens dissipation lengthscale (average over last half of data) = %f ',eta_ens));

% averge R_l_uu to a number
R_l_uu = R_l_uu(length(R_l_uu)/2:length(R_l_uu));
R_l_uu = sum(R_l_uu)/length(R_l_uu);
disp(sprintf('R_l_uu (average over last half of data) = %f ',R_l_uu));
% averge R_l_uv to a number
R_l_uv = R_l_uv(length(R_l_uv)/2:length(R_l_uv));
R_l_uv = sum(R_l_uv)/length(R_l_uv);
disp(sprintf('R_l_uv (average over last half of data) = %f ',R_l_uv));
% averge R_l_vv to a number
R_l_vv = R_l_vv(length(R_l_vv)/2:length(R_l_vv));
R_l_vv = sum(R_l_vv)/length(R_l_vv);
disp(sprintf('R_l_vv (average over last half of data) = %f ',R_l_vv));

R_l_uu_inf = R_l_uu_inf(length(R_l_uu_inf)/2:length(R_l_uu_inf));
R_l_uu_inf = sum(R_l_uu_inf)/length(R_l_uu_inf);
disp(sprintf('R_l_uu_inf (average over last half of data) = %f ',R_l_uu_inf));
R_l_uv_inf = R_l_uv_inf(length(R_l_uv_inf)/2:length(R_l_uv_inf));
R_l_uv_inf = sum(R_l_uv_inf)/length(R_l_uv_inf);
disp(sprintf('R_l_uv_inf (average over last half of data) = %f ',R_l_uv_inf));
% averge R_l_vv_inf to a number
R_l_vv_inf = R_l_vv_inf(length(R_l_vv_inf)/2:length(R_l_vv_inf));
R_l_vv_inf = sum(R_l_vv_inf)/length(R_l_vv_inf);
disp(sprintf('R_l_vv_inf (average over last half of data) = %f ',R_l_vv_inf));

% averge R_l_ke to a number
%R_l_ke = R_l_ke(length(R_l_ke)/2:length(R_l_ke));
%R_l_ke = sum(R_l_ke)/length(R_l_ke);
%disp(sprintf('R_l_ke (average over last half of data) = %f ',R_l_ke));

disp(sprintf('1/250 in units of eta:  %f',(1/250)/eta));
disp(sprintf('1/500 in units of eta:  %f',(1/500)/eta));
disp(sprintf('1/512 in units of eta:            %f   %f', (1/512)/eta,eta*2*pi*512*sqrt(2)/3 )) ;
disp(sprintf('1/nx in units of eta:  %f',(1/nx)/eta));

tturn=-2*ke./(ke_diss_d);
tturn = tturn(length(tturn)/2:length(tturn));
tturn = sum(tturn)/length(tturn);
disp(sprintf('eddy turnover time (averaged over last haf of data) = %f ',tturn));

tturn2=ints(7,:)./(ints(12,:)); % enstrophy/enstrophy disipation rate
tturn2 = tturn2(length(tturn2)/2:length(tturn2));
tturn2 = sum(tturn2)/length(tturn2);
disp(sprintf('EDDY TURNOVER TIME (averaged over last haf of data) = %f ',tturn2));

tturn3=2*ke./epsilon;
tturn3 = tturn3(length(tturn2)/2:length(tturn3));
tturn3 = sum(tturn3)/length(tturn3);
disp(sprintf('eddy turnover time (averaged over last haf of data) = %f ',tturn3));

tturn4=sqrt(ones(size(ens))./ens);
tturn4 = tturn4(length(tturn4)/2:length(tturn4));
tturn4 = sum(tturn4)/length(tturn4);
disp(sprintf('eddy turnover time (averaged over last haf of data) = %f ',tturn4));

tturn5=sqrt(ens);
tturn5 = tturn5(length(tturn5)/2:length(tturn5));
tturn5 = sum(tturn5)/length(tturn5);
disp(sprintf('eddy turnover time (averaged over last haf of data) = %f ',tturn5));

epsilon = epsilon(length(epsilon)/2:length(epsilon));
epsilon = sum(epsilon)/length(epsilon);
disp(sprintf('epsilon_E (averaged over last haf of data) = %f ',epsilon));

%epsilon_ke = epsilon_ke(length(epsilon_ke)/2:length(epsilon_ke));
%epsilon_ke = sum(epsilon_ke)/length(epsilon_ke);
%disp(sprintf('epsilon_ke (averaged over last haf of data) = %f ',epsilon_ke));

%ens_diss = ens_diss(length(ens_diss)/2:length(ens_diss));
%ens_diss = sum(ens_diss)/length(ens_diss);
%disp(sprintf('ens_diss (averaged over last haf of data) = %f ',ens_diss));
%print -depsc scalars.ps

%figure(15)
%clf
%hold on

%plot(time_2,ens_diss_tot,'r') 
%title('d/dt enstrophy');
%hold off
%xlabel('time')

     
