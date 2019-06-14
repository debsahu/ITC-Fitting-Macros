function [  ] = Fcn_Fit_Cp_Ndil( ListFile  )
%function [ dQm_fit, dQm_exp, x, I0, TitrPar, p ] = Fcn_Fit_Cp_Ndil()
%function [  ] = Fcn_Fit_Cp_Ndil(   )
% fits ITC for datafiles in ListFile to Baldwin model & treats
% dilution/endpoint effects
%   

T_rm = 298.0;

%%%%%%% get data and titration info
%
%% Gets data for fitting itc calc.  
% TitrPar(i,:) = [ tau(i), V0(i), cLs(i), ndp(i), Nfit(i) ]; w/
%  tau(i):  tau = T / T_rm for titr i 
%   V0(i):  working vol for titr i (MUST be in same units as dV!!!)
%  cLs(i):  concen of Ligand in injection dV (MUST be same units as c_Mt!!!)
%  ndp(i):  number of initial datapts to be discarded for titr i
% Nfit(i):  number of total datapts in titr i for fitting 
%
% RawData - concatenated array with all data from all titrations
%          
%   I0(i) =  [ dV_dis(i)  c_Mt_dis(i)   f_Lt_dis(i) ] for titration i
%         -  corresponds to Deb's I0
%   dV_dis:  Vol of last discarded injection 
% c_Mt_dis:  Macromol concen (mt) in V0 after last discarded injection
% f_Lt_dis:  c_Ltot/c_Mtot (xmt) ratio of total ligand to total Macromol
%            concen in working vol V0 after last discarded injection
%
%    x    - concatenated array of data (x_1) for fitting from each titration
%   
%  	x_1 = [ dV(Nstart:Nend), c_Mt(Nstart:Nend), f_Lt(Nstart:Nend) ]
%   dV  : volume of last injection
%  c_Mt : Macromol concen (mt) in V0 after last injection
%  f_Lt : c_Ltot/c_Mtot (xmt) ratio of total ligand to total Macromol
%         concen in working vol V0 after last discarded injection
% 
% dQm_exp_fit - concatenated array of measured exptl heats for fitting 
%             - same number of rows as x

%Read Input paramters from list file
fidl = fopen(ListFile);
% Each line corresponds to:
% filename  Temp   V0   cLs     ndp
C = textscan (fidl, '%s %f %f %f %d' );
DataFiles = C{1};
Tempp = double( C{2} ); 
Vol0 = double( C{3} );
cLs  = double( C{4} ); 
ndp  = double( C{5} );
Temp = [ Tempp, Vol0, 0.001*cLs, ndp ];
fclose(fidl);
%End of reading Input paramters from list file

N_DataFiles = size(DataFiles,1);
TitrPar = [ Temp, ones(N_DataFiles,1) ];

% TitrPar(i,:) = [ Temp(i), V0(i), cLs(i), ndp(i), Nfit(i) ]; w/
% Temp(i):  temp for titr i (to be rescaled below)
%   V0(i):  working vol for titr i (MUST be in same units as dV!!!)
%  cLs(i):  concen of Ligand in injection dV (MUST be same units as c_Mt!!!)
%  ndp(i):  number of initial datapts to be discarded for titr i
% Nfit(i):  number of total datapts in titr i for fitting (DETERMINED BELOW !!!!) 

for i_DataSet = 1:N_DataFiles       
%     DataFiles{i_DataSet};   
    RawData1 = importdata( DataFiles{i_DataSet} );    
    RawData1 = RawData1.data;   
    RawData1(:,6) = 0.001 * RawData1(:,6); % convert to kcal   
    dV   = RawData1(:,2);		% injection vol  (injv)
								% nb i do NOT use col 3 (xt)
    c_Mt = 0.001 * RawData1(:,4); % converted from mM to M   
    f_Lt = RawData1(:,5);		% curr c_Ltot / c_Mtot (xmt)
    dQm  = RawData1(:,6);		% measured heat (ndh)  
    N_dis = int64( TitrPar(i_DataSet,4) ); % num of inital data pts NOT used in fit   
    % make sure that at least 1 data pt is discarded from titration 
    if N_dis < 1
         disp( 'Error in N_dis.' );
         exit; 
    end 
        
    I0_1  = [ dV(N_dis)  c_Mt(N_dis+1)   f_Lt(N_dis) ];  % I0: initial injection
    N_tot = size( f_Lt, 1 ); % total num data points
    Nstart = N_dis+1;
    Nend   = N_tot-1;        % last line only contains cLt and cMt, NO dQm  
    Nfit   = Nend - Nstart + 1;
    % nb shift in data to deal with 
    % modified.  dV, f_Lt, and dQm_q are all aligned in correct row.
    %            c_Mt and c_Lt are not
	x_1 = [ dV(Nstart:Nend), c_Mt(Nstart+1:Nend+1), f_Lt(Nstart:Nend) ];
    dQm_1 = dQm(Nstart:Nend); 
    TitrPar(i_DataSet,5) = double( Nfit );  % tot # data pts for fitting in i_DataSet
	% rescale temperature: tau = Temp / 298
	TitrPar(i_DataSet,1) = TitrPar(i_DataSet,1) / T_rm; 
    %%%%%%%% !! shift heats to zero!
    N_shift = 5;    % parameter determining no datapoints used in shift
    dQ_shift = mean( dQm_1(Nfit-N_shift+1:Nfit) );
    dQm_1 = dQm_1 - dQ_shift ;
    RawData1(:,6) = RawData1(:,6) - dQ_shift;
    %%%%%%%% !! end shift heats
    if ( i_DataSet == 1 ) 
         x   = [   x_1 ];
		 I0  = [  I0_1 ];
         dQm_exp = [ dQm_1 ];
	     RawData = [ RawData1 ];
    else 
         x  = vertcat(  x,  x_1 ); 
         I0 = vertcat( I0, I0_1 ); 
		 dQm_exp = vertcat( dQm_exp, dQm_1);
		 RawData = vertcat( RawData, RawData1 ); 
    end 
end
%
% 1. reads ListFile with list of data files and parameters for each titr 
% 2. reads each data files
%	a. vertcats data files into RawData array
%	b. truncates datafiles into x_1 and dQm to include only points used for fitting:
%	   x_1(i,:) = [ dV(i)  c_Mt(i)  f_Lt(i) ] for datapt i in given titr
%      dQm_1 = dQm(Nstart:Nend)
%      nb dQ is scaled by 0.001 to get results in kcal
%   c. stores last of the initial datapoints into I0_1 
%   d. vertcats dQm_1, x_1, and I0_1 over titrations into dQm_exp_fit, x, and I0
%   e. vertcats TitrPar over titrations:
%      TitrPar(i,:) = [ tau(i), V0(i), cLs(i), ndp(i), nt(i) ] for titr i 
%	   nb tau(i) = Temp(i) / T_rm where T_rm is input
%         T_rm = room temp, option allows for rescaling or not
%

%% setup parameters for nonlin fitting
% initialize parameters 
param = ones( 5, 1 ); 

% initial estimates
c_h = -75.0  ;   % kcal/mol by def R in Fcn_Calc_dQm_N
t_h = 295.0 / T_rm; % T_H=295K Baldwin PNAS'86; T_rm above
t_s = 386.0 / T_rm; % T_S=386K Baldwin PNAS'86; T_rm above
h_0 =   0.0 ;   % kcal/mol by def R in Fcn_Calc_dQm_N
s_0 = -10.0 ;   % kcal/mol by def R in Fcn_Calc_dQm_N

% fitting parameters - hydrophobic model
param(1) = c_h;       % c_h = T_rm Delta_r c_Hydrophobic
param(2) = t_h;       % t_h = T_H / T_rm    ! fixed by Baldwin model
param(3) = t_s;       % t_s = T_S / T_rm    ! fixed by Baldwin model

% fitting parameters  - secondary effects
param(4) = h_0;      % h_0 = Delta_r H_0 (nonhydrophobic H)
param(5) = s_0;      % s_0 = T_rm Delta_RT S
param0 = param;

% Id parameters that will be fixed (true) or vary (false) in fitting
fixed = [ false true  true  false false  ]';
%           1     2     3     4     5   

%%%%%%% solve for parameters 
[ p, err ] = nlinfit_Itc( param, fixed, TitrPar, I0, x, T_rm, dQm_exp );

%%%%%%% output
dQm_fit = Fcn_Calc_dQm_N( p, TitrPar, I0, x, T_rm ); 

%% Write the Fits
i_Data_N = 1;   % initialize index for total data array
                % nb one index loops over all data from all files
for i_DataSet = 1:size(TitrPar,1)       % loop over data files
    % name output file based on name of corresponding input file
    OutFile = strcat( 'Fit.', DataFiles{i_DataSet}, '.dat' );
    fid1 = fopen( OutFile, 'w' );    % open file
    N_DataPts_iSet = int64( TitrPar(i_DataSet,5) );
    
    for i_Data_1 = 1:N_DataPts_iSet % loop over pts in given ITC dataset       
        f_Lt_1 = x(i_Data_N,3);     % c_L_tot/c_M_tot
        dQm_exp_1 = dQm_exp(i_Data_N);
        dQm_fit_1 = dQm_fit(i_Data_N); 
        % this is for final formatted output
        fprintf(fid1, '%10.5f \t %10.5f \t %10.5f \n', ...
                      f_Lt_1, dQm_fit_1, dQm_exp_1 );
        i_Data_N = i_Data_N + 1;    % incr loop over all data points & files
    end
    fclose( fid1 );
end 
%% Write Output Parameters
%                Fcn_OutputParam( param0, p, err, T_rm, ListFile, 'OutputParSumm.txt' ); 
% function Out = Fcn_OutputParam( param0, param, err, T_rm, ListFile, OutFile )
% outputs parameters and error estimates for fitting
%  param0 - initial estimates for parameters
%  param  - final estimates for parameters
%         - [ c_h  t_h  t_s  h_0  s_0] in Baldwin model t_h, t_s fixed
%  err    - estimates of uncertainty in fit parameters
%  T_rm   - room temperature used as parameter for rescaling in fitting
%  ListFile - list of ITC data files for use in fitting
%  OutFile  - file name for output of function
% NB wgn06Sep2013 - converted cMt (Fcn_GetDataN.m) and cLs
%                 (Fcn_GetDataFilesPar.m) from mM to M
% nb wgn-28Sept2013 - corrected value for R

% Temporary Functions
%DH = D_r H(T); t=T/T_rm, c_h = T_rm * D_r c_Hphobic, t_h = T_H/T_rm,
% h_0=nonhydrobic contribution to enthalpy
Fcn_Trm_DSt = @(t,c_h,t_s,s_0) c_h * log( t/t_s ) + s_0;
%Trm_DSt = T_rm * D_r S(T); t=T/T_rm, c_h = T_rm * D_r c_Hphobic, t_s = T_S/T_rm,
% s_0 = T_rm * D_r S_RT (rot/trans contribution to entropy)
Fcn_DHt = @(t,c_h,t_h,h_0 ) c_h * ( t - t_h ) + h_0;

fid2 = fopen( 'OutputParSumm.txt', 'w' ); 
fprintf( fid2, 'Input file:\t %12s\n', ListFile ); 
fprintf( fid2, 'calculation performed on Date Year  Month  Date  Time:\n' );
fprintf( fid2, '%10.2f', clock ); 
fprintf( fid2, '\n'); 
fprintf( fid2, '\n --- Final transformed parameters --- \n' ); 

Dr_Chydro = p(1) / T_rm; 
T_H       = p(2) * T_rm; 
T_S       = p(3) * T_rm; 
H_0       = p(4);
S_0       = p(5) / T_rm; 

kT_rm = 0.592187 * T_rm / 298.0;   % in kcal/mol 

fprintf( fid2, 'Room Temp T_rm: %10.5f K\n', T_rm ); 
fprintf( fid2, 'Hydrophobic heat capacity.\n' ); 
fprintf( fid2, 'Delta_R C_hydrophobic: %10.5f ± %10.5f kcal/mol \n', Dr_Chydro, err(1)/T_rm ); 
fprintf( fid2, ' T_H                 : %10.5f K \n', T_H ) ;
fprintf( fid2, ' T_S                 : %10.5f K \n', T_S ) ; 
fprintf( fid2, 'NonHydrophobic enthalpy.\n' ); 
fprintf( fid2, ' H_0(T_rm)           : %10.5f ± %10.5f kcal/mol \n', H_0, err(4) ) ; 
fprintf( fid2, 'NonHydrophobic entropy (Rot + Trans + Conf?) at T_rm.\n' ); 
fprintf( fid2, ' Delta_R S_x(T_rm)   : %10.5f ± %10.5f kcal/mol K \n', S_0, err(5)/T_rm ) ; 
fprintf( fid2, 'Contrib of Delta_R S_x to Delta_R G at T_rm.\n' ); 
fprintf( fid2, ' -T_rm * Delta_R S_x : %10.5f kcal/mol \n', -1.* p(5) ) ; 
fprintf( fid2, '\n\n --- Initial / Final parameters computed --- \n' ); 
fprintf( fid2, 'i \t param0(i) \t\t param(i) \t\t Conf Int.\n' ); 

for i = 1:size(p,1)
    fprintf( fid2, '%d \t %10.5f \t %10.5f \t %10.5f \n', i, param0(i), p(i), err(i) ); 
end

fprintf( fid2, '\n\n --- Thermo properties for each file ---\n' );
fprintf( fid2, ' Temp in K.  dG, dH, -T dS in kcal/mol \n' ); 
fprintf( fid2, 'file \t\t\t name \t\t\t     Temp      dG ' );
fprintf( fid2, '     dH       dS      -TdS         K \n'); 
for i = 1:size(Temp,1)
    T     = Temp(i);
    Fname = DataFiles{i};
    fprintf( fid2, '%3d %30s %7.2f', i, Fname, T ); 
    tau = T / T_rm; 
    dH  =     Fcn_DHt( tau, p(1), p(2), p(4) ); 
    dS  = Fcn_Trm_DSt( tau, p(1), p(3), p(5) ) / T_rm; 
    TdS = T * dS; 
    dG = dH - TdS; 
    beta = 1.0 / ( tau * kT_rm ); 
    Kt = exp( -1.0 * beta * dG ); 
    fprintf( fid2, '%8.3f %8.3f %8.3f %8.3f %12.3e\n', dG, dH, dS, -1.0 * TdS, Kt ); 
end

fclose(fid2);

fprintf('\nFitting Done. \n Please check OutputParSumm.txt \n for fitted values of dCp, dH0 and dS0.\n\n');

return; 

end

