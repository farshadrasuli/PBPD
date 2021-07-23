%{
Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.
%}
clear all;
close all;
clc;

%==================================================
% Input Data
%==================================================
    ProjectName = "Liao";
              T = 0.81; % Fundamental period in s
              C = 0.092; % Seismic response coefficient from code
              k = 1.00;  % Height distribution exponent
            R_u = 7.5;
         S_aDBE = 0.74;
         S_aMCE = 1.11;
     theta_yDBE = 0.005; % Yield Drif Ratio for DBE
     theta_yMCE = 0.005; % Yield Drif Ratio for MCE
     theta_uDBE = 0.02; % Target Drift Ratio for DBE
     theta_uMCE = 0.03; % Target Drift Ratio for MCE
            rho = 0.00; % Ratio of Stiffness Hardening after Yielding
              g = 32.1741; % Acceleration of gravity: 9.80665 m/s2 or 32.1741 ft/s2
WeightOfStories = { % Exclude seismic weight of base level
'Roof',	     54,	518
'Story 4',	 41,	519
'Story 3',	 28,	519
'Story 2',	 15,	519
};
%==================================================
% Now, You can run.
%==================================================






% Output Name
OutputName = ['PBPD_RC_' ProjectName];

%--------------------------------------------------
% Determining design base shear
%--------------------------------------------------

% Modification factor (C_2) to represent the effect of pinched hysteresis
if R_u < 3
  if T >= 0.2 && T < 0.4
    C_2 = 2.5 - 6.5 * (T - 0.2);
  endif
  if T >= 0.4
    C_2 = 1.1 - 0.077 * (T - 0.4);
  endif
elseif R_u >= 3
  if T >= 0.2 && T < 0.4
    C_2 = 3 - 7.5 * (T - 0.2);
  endif
  if T >= 0.4 && T < 0.8
    C_2 = 1.5 - T + 0.4;
  endif
  if T >= 0.8
    C_2 = 1.1 - 0.045 * (T - 0.8);
  endif
endif

% Modified theta_u, modifiedTheta_u
modifiedTheta_uDBE = theta_uDBE / C_2;
modifiedTheta_uMCE = theta_uMCE / C_2;

% Inelastic Drift Ratio, theta_p
theta_pDBE = modifiedTheta_uDBE - theta_yDBE; 
theta_pMCE = modifiedTheta_uMCE - theta_yMCE;

% Ductility factor, mu_s
mu_sDBE = modifiedTheta_uDBE / theta_yDBE;
mu_sMCE = modifiedTheta_uMCE / theta_yMCE;

% Ductility Reduction Facrot, R_mu
if T >= 0 & T < (0.57/10)
  R_muDBE = 1;
  R_muMCE = 1;
endif
if T >= (0.57/10) & T < (0.57/4)
  R_muDBE = sqrt(2 * mu_sDBE - 1) * ( ( 0.57 / (4 * T) ) ^ ( 2.513 * log10( 1 / sqrt(2 * mu_sDBE) ) ) );
  R_muMCE = sqrt(2 * mu_sMCE - 1) * ( ( 0.57 / (4 * T) ) ^ ( 2.513 * log10( 1 / sqrt(2 * mu_sMCE) ) ) );
  endif
if T >= (0.57/4) & T < (0.57 * sqrt(2 * mu_sDBE - 1) / mu_sDBE)
  R_muDBE = sqrt(2 * mu_sDBE - 1);
endif
if T >= (0.57/4) & T < (0.57 * sqrt(2 * mu_sMCE - 1) / mu_sMCE)
  R_muMCE = sqrt(2 * mu_sMCE - 1);
endif
if T >= (0.57 * sqrt(2 * mu_sDBE - 1) / mu_sDBE) & T < 0.57
  R_muDBE = T * mu_sDBE / 0.57;
endif
if T >= (0.57 * sqrt(2 * mu_sMCE - 1) / mu_sMCE) & T < 0.57
  R_muMCE = T * mu_sMCE / 0.57;
endif
if T >= 0.57
  R_muDBE = mu_sDBE;
  R_muMCE = mu_sMCE;
endif

% Energy Modification Factor, gamma
gamma_DBE = (2 * mu_sDBE - 1 + rho * (mu_sDBE - 1) ^ 2) / (R_muDBE ^ 2);
gamma_MCE = (2 * mu_sMCE - 1 + rho * (mu_sDBE - 1) ^ 2) / (R_muMCE ^ 2);

% Computing Number of stories
NumberOfStory = size(WeightOfStories,1);

%::::::::::::::::::::::::::::::::::::::::::::::::::
% Parameters for each elevation
%::::::::::::::::::::::::::::::::::::::::::::::::::

% Taking label of elevation, h_i in meter or foot, and w_i in kN or kips
Parameters = WeightOfStories;
% Computing Total Weight
Weight = 0;
for i = 1:NumberOfStory
  Weight = Weight + Parameters{i,3};
endfor
% Computing w_i * h_i in kN-m or kips-ft
for i = 1:NumberOfStory
  Parameters{i,4} = Parameters{i,2} * Parameters{i,3};
endfor
% Computing SIGMA (w_i * h_i) from i to top in kN-m or kips-ft
S = 0;
for i = 1:NumberOfStory
  S = S + Parameters{i,4};
  Parameters{i,5} = S;
endfor
clear S;
% Computing beta_i
for i = 1:NumberOfStory
  Parameters{i,6} = ( Parameters{i,5} / Parameters{1,4} ) ^ (0.75 * (T^-0.2));
endfor
% Computing beta_i - beta_(i+1)
Parameters{1,7} = Parameters{1,6};
for i = 2:NumberOfStory
  Parameters{i,7} = Parameters{i,6} - Parameters{i-1,6};
endfor
% Computing alpha
S1 = 0;
for i = 1:NumberOfStory
  S1 = S1 + Parameters{i,7} * Parameters{i,2};
endfor
S2 = (Parameters{1,4} / Parameters{NumberOfStory,5}) ^ (0.75 * (T^-0.2));
alpha_DBE = S1 * S2 * (theta_pDBE * 8 * (pi^2) / ( (T^2) * g ) );
alpha_MCE = S1 * S2 * (theta_pMCE * 8 * (pi^2) / ( (T^2) * g ) );
% Computing V/W, C_v
C_vDBE = (-alpha_DBE + sqrt( (alpha_DBE^2 + 4 * gamma_DBE * S_aDBE^2 ) )) / 2;
C_vMCE = (-alpha_MCE + sqrt( (alpha_MCE^2 + 4 * gamma_MCE * S_aMCE^2 ) )) / 2;
% Design Base Shear without P-Delta
V_DBE = C_vDBE * Weight;
V_MCE = C_vMCE * Weight;
% Design Base Shear with P-Delta
V_PDeltaDBE = V_DBE + Weight * theta_uDBE;
V_PDeltaMCE = V_MCE + Weight * theta_uMCE;
% Maximum Design Base Shear with P-Delta
V = V_PDeltaDBE; %max(V_PDeltaDBE,V_PDeltaMCE);

%..................................................
% Creating .xlsx file for Base Shear
%..................................................
BaseShear = {
'Parameters'               , 'DBE'              , 'MCE'
'--------'                 , '--------'         , '--------'
'T'                        , T                  , T
'S_a'                      , S_aDBE             , S_aMCE
'theta_y'                  , theta_yDBE         , theta_yMCE
'theta_u'                  , theta_uDBE         , theta_uMCE
'Modification factor, C_2' , C_2                , C_2
'Modified theta_u'         , modifiedTheta_uDBE , modifiedTheta_uMCE
'theta_p'                  , theta_pDBE         , theta_pMCE
'mu_s'                     , mu_sDBE            , mu_sMCE
'R_mu'                     , R_muDBE            , R_muMCE
'gamma'                    , gamma_DBE          , gamma_MCE
'alpha'                    , alpha_DBE          , alpha_MCE
'V/W'                      , C_vDBE             , C_vMCE
'V w/o P-delta'            , V_DBE              , V_MCE
'V w/ P-delta'             , V_PDeltaDBE        , V_PDeltaMCE
'--------'                 , '--------'         , '--------'
'Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.', '', ''
};
%xlswrite(OutputName, BaseShear,'Base Shear');

%..................................................
% Creating .xlsx file for Lateral Forces
%..................................................

for j = 1:9
  for i = 1:2
    LateralForces{i,j} = '';
  endfor
endfor
for i = 1:9
  LateralForces{3,i} = '-----';
endfor
LateralForces{1,4} = 'PBPD';
LateralForces{2,1} = 'Elevation';
LateralForces{2,2} = 'Height';
LateralForces{2,3} = 'Weight';
LateralForces{2,4} = 'beta_i';
LateralForces{2,5} = 'beta_i - beta_i+1';
LateralForces{2,6} = 'F_i';
LateralForces{2,7} = 'V_i';
LateralForces{1,8} = 'Code';
LateralForces{2,8} = 'F_ui';
LateralForces{2,9} = 'V_ui';
% Taking label of elevation, height in meter or foot, and weight in kN or kips
for j = 1:3
  for i = 1:NumberOfStory
    LateralForces{3+i,j} = Parameters{i,j};
  endfor
endfor
% Taking beta_i
for i = 1:NumberOfStory
  LateralForces{3+i,4} = Parameters{i,6};
endfor
% Computing Sigma beta_i
SigmaBeta_i = 0;
for i = 1:NumberOfStory
  SigmaBeta_i = SigmaBeta_i + Parameters{i,6};
endfor
% Taking beta_i - beta_i+1
for i = 1:NumberOfStory
  LateralForces{3+i,5} = Parameters{i,7};
endfor
% Computing Sigma (beta_i - beta_i+1)
SigmaBeta2 = 0;
for i = 1:NumberOfStory
  SigmaBeta2 = SigmaBeta2 + Parameters{i,7};
endfor
% Computing F_i for PBPD
for i = 1:NumberOfStory
  LateralForces{3+i,6} = Parameters{i,7} * S2 * V;
endfor
% Computing V_i for PBPD
for i = 1:NumberOfStory
  LateralForces{3+i,7} = ( Parameters{i,5} / Parameters{NumberOfStory,5} ) ^ (0.75 * (T^-0.2)) * V;
endfor
% Computing base shear for Code Design
V_u = C * Weight;
% Computing w_i * h_i ^ k
for i = 1:NumberOfStory
  Parameters{i,8} = Parameters{i,3} * Parameters{i,2}^k ;
endfor
% Computing F_ui for Code
S = 0;
for i = 1:NumberOfStory
  S = S + Parameters{i,8};
endfor
for i = 1:NumberOfStory
  LateralForces{3+i,8} = Parameters{i,8} / S * V_u;
endfor
clear S
% Computing V_ui for Code
S = 0;
for i = 1:NumberOfStory
  S = S + LateralForces{3+i,8};
  LateralForces{3+i,9} = S;
endfor
clear S
% Computing SIGMA
for i = 1:9
  LateralForces{3+NumberOfStory+1,i} = '';
endfor
LateralForces{3+NumberOfStory+1,1} = 'sum total';
LateralForces{3+NumberOfStory+1,3} = Weight;
LateralForces{3+NumberOfStory+1,4} = SigmaBeta_i;
LateralForces{3+NumberOfStory+1,5} = SigmaBeta2;
LateralForces{3+NumberOfStory+1,6} = V;
LateralForces{3+NumberOfStory+1,8} = V_u;
% Adding copyright statement
LateralForces{3+NumberOfStory+2,1} = 'Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.';
for j = 2:9
  LateralForces{3+NumberOfStory+2,j} = '';
endfor
%xlswrite(OutputName, LateralForces,'Lateral Forces');
clearvars -except Weight BaseShear LateralForces
%{
Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.
%}
