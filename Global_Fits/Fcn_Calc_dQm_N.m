function dQm_fit = Fcn_Calc_dQm_N( param, TitrPar, I0, x, T_rm )
%Fcn_Calc_dQm_N uses Baldwin model to estimate molar multiple ITC heats
%   uses hydrophobic model of Baldwin to estimate molar heats for ITC
%   titrations.  addresses dilution effects and end point effects.
%   requires concentrations should be expressed in M for comparison with K
%   which is implicitly calculated in M^-1
%   wgn-28Sept2013 corrects value for RT at T=298  

% Temporary functions
%Fcn_Trm_DSt = T_rm * D_r S(T); t=T/T_rm, c_h = T_rm * D_r c_Hphobic, t_s = T_S/T_rm,
% s_0 = T_rm * D_r S_RT (rot/trans contribution to entropy)
Fcn_Trm_DSt = @(t,c_h,t_s,s_0) c_h * log( t/t_s ) + s_0;
%Fcn_DHt = D_r H(T); t=T/T_rm, c_h = T_rm * D_r c_Hphobic, t_h = T_H/T_rm,
% h_0=nonhydrobic contribution to enthalpy
Fcn_DHt = @(t,c_h,t_h,h_0 ) c_h * ( t - t_h ) + h_0;
% Fcn_f_ML = c_LT/c_Mt, Mt_KT = Mt * KT, f_ML = c_ML / c_Mt (fraction ML bound)
%   NB c_Mt assumed fixed b/c vol assumed const
Fcn_f_ML = @(f_Lt,Mt_KT ) (0.5 * ( 1 + f_Lt ) + 0.5 ./ Mt_KT) ...
    - sqrt( ((0.5 * ( 1 + f_Lt ) + 0.5 ./ Mt_KT).* ...
    (0.5 * ( 1 + f_Lt ) + 0.5 ./ Mt_KT) - f_Lt) ); 

% global parameters
% fitting parameters - hydrophobic model
c_h   = param(1);     % c_h = T_rm Delta_r c_Hydrophobic
t_h   = param(2);     % t_h = T_H / T_rm   ----- fixed in Baldwin model
t_s   = param(3);     % t_s = T_S / T_rm   ----- fixed in Baldwin model

% fitting parameters  - secondary effects
h_0   = param(4);     % h_0 = Delta_r H_0 (nonhydrophobic H)
s_0   = param(5);     % s_0 = T_rm Delta_RT S

kT_rm = 0.592187 * T_rm / 298.0;   % input

N_fit_tot = 0;        % initialize no of heats calculated for counting

N_DataFiles = size( TitrPar, 1 );

% loop over titrations
for i_DataSet = 1:N_DataFiles
    
    % get titration-specific parameters
    t     = TitrPar(i_DataSet,1);   % rescaled temp.  t = T/T_rm
    V_0   = TitrPar(i_DataSet,2);   % working volume
    c_Ls  = TitrPar(i_DataSet,3);   % titrant concen in syringe

    % nb these are stored as doubles/reals and converted to ints
    N_dis = int64( TitrPar(i_DataSet,4) );  % # datapts discarded. NOT needed here
    N_fit = int64( TitrPar(i_DataSet,5) );  % # datapts for fitting 
    
    N_start = N_fit_tot+1;          % 1st  pt for fitting i_DataSet
    N_end   = N_start + N_fit - 1;  % last pt for fitting i_DataSet
    
    N_fit_tot = N_fit_tot + N_fit;  % tot # of fit pts including i_DataSet
    
    %i_DataSet
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % calc dH, dS, dG, K - fixed by temp (t) and model parameters
    % DHt = Delta_R H(T) 
    DHt = Fcn_DHt( t, c_h, t_h, h_0 );   

    % Trm_DSt = T_rm * Delta_R S(T)
    Trm_DSt = Fcn_Trm_DSt( t, c_h, t_s, s_0 );   

    % b_rm = beta_rm = 1 / k_B T_rm
    bt_rm = 1. / kT_rm; 

    % bt_DGt = beta Delta_R G(T) 
    bt_DGt = bt_rm * ( DHt / t - Trm_DSt );

    % Kt = K(T) = exp[ -beta Delta_R G(T) ]
    Kt = exp( -1. * bt_DGt );  
    
    %return
   
    %%%%%%%%%%%%%%%%%%%%%%%%
    % parse the data into titration-specific arrays!s
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % deal w/ last discarded data pt for i_DataSet
    I0_1 = I0( i_DataSet, : );
    
    dV_0   = I0_1(1);   % vol of last discarded inject.  NOT used
    c_Mt_0 = I0_1(2);   % Macromol concen for last discarded pt.  USED
    f_Lt_0 = I0_1(3);   % c_L_tot/c_M_tot for last discarded pt.  USED  

    %%%%%               % get heat evolved from last discarded data pt
    Mt_Kt = c_Mt_0 * Kt;                % c_Mt * K(T) for last disc datapt
    
    f_ML_0 = Fcn_f_ML( f_Lt_0, Mt_Kt ); % f_ML = = c_ML / c_M_tot
    
    Q_0 = V_0 * c_Mt_0 * DHt * f_ML_0;  % total heat evolved 
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % data for fitting from i_DataSet
    x_1 = x(N_start:N_end,:);
    
    dV   = x_1(:,1);    % vol injected at each pt
    c_Mt = x_1(:,2);    % total Macromol concen at each pt (c_M_tot)
    f_Lt = x_1(:,3);    % c_L_tot / c_M_tot at each pt
    
    
    f_ML = Fcn_f_ML( f_Lt, Kt .* c_Mt ); % = c_ML / c_M_tot
    
    Q = V_0 .* c_Mt .* DHt .* f_ML;      % total heat evolved
        
    %%%%%%%%%%%%%%%%%%%%%%%%
    % calc differential heats, dQ, and deal with dilution effects  
    
    dQ = zeros( N_fit, 1 );
    
    % treat first titration separately
    dQ(1) = Q(1)                                        ...
          + ( dV(1) / V_0 ) * 0.5 * ( Q(1) + Q_0 )      ...
          - Q_0; 

    % treat remaining titration points      
    for i = 2:N_fit

        dQ(i) = Q(i)                                        ...
 	          + ( dV(i) / V_0 ) * 0.5 * ( Q(i) + Q(i-1) )   ...
              - Q(i-1); 
          
    end 
          
    % differential heat per mol titrant added
    
    dQm = dQ ./ ( c_Ls .* dV );      % calc dQ_m = dQ_tot/dL

    %%%%%%%%%%%%%%%%%%%%%%%%    
    % store into total dQ_m array

    if ( i_DataSet == 1 )

       dQm_fit = [ dQm ];  

    else

       dQm_fit = vertcat( dQm_fit, dQm ); 

    end
    
end    

end

