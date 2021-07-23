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
            Psi = 1.1;  % The over-strength factor for plastic hinge of the columns at the base
              x = 2.1;  % The ratio of the absolute value of M_pb_negative to M_pb_positive
            Ksi = 1.25; % The over-strength factor for plastic hinge of the beams
BayWidth = [ % Width of resisitant bays from left to right in meter or foot
30,30,30
];
w_tributary = [ % Factored gravity loads on beams in every bay at every elevation, kips/ft or kN/m
5.76, 5.76, 5.76
5.76, 5.76, 5.76
5.76, 5.76, 5.76
5.76, 5.76, 5.76
];
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

%--------------------------------------------------
% Design of DYMs
%--------------------------------------------------

% Computing number of resistant bays
B = size(BayWidth,2);
% Computing required plastic moment of the first-story columns
M_pc = (Psi * V * Parameters{NumberOfStory,2}) / (2 * B * (B+1));
% Determining force of each bay in every elevation
temp1 = 0;
for i = 1:B
  temp1 = temp1 + BayWidth(1,i) ^ 3;
endfor
for i = 1:B
  temp2(1,i) = BayWidth(1,i) ^ 3 / temp1;
endfor
for j = 1:B
  for i = 1:NumberOfStory
    ForceInBays{i,j} = temp2(1,j) * LateralForces{3+i,6};
  endfor
endfor
clear temp1 temp2
% Computing M_pb_positive and M_pb_negative
for j = 1:B
  for i = 1:NumberOfStory
    temp1(i,j) = ForceInBays{i,j} * Parameters{i,2};
  endfor
endfor
for j = 1:B
  S = 0;
  for i = 1:NumberOfStory
    S = S + temp1(i,j);
    temp2(1,j) = S;
  endfor
endfor
for j = 1:B
  for i = 1:NumberOfStory
    M_pb_positive(1,j) = (temp2(1,j) - 2 * M_pc) / ((1+x) / 0.9 * SigmaBeta_i);
  endfor
endfor
M_pb_negative = -1 * x .* M_pb_positive;
for j = 1:B
    M_pb(1,j) = max(abs(M_pb_positive(1,j)),abs(M_pb_negative(1,j)));
endfor
% Computing M_pb_positive for each bay
for j = 1:B
  for i = 1:NumberOfStory
    M_pb_posInStories{i,j} = Parameters{i,6} * M_pb_positive(1,j);
  endfor
endfor
% Computing M_pb_negative for each bay
for j = 1:B
  for i = 1:NumberOfStory
    M_pb_negInStories{i,j} = Parameters{i,6} * M_pb_negative(1,j);
  endfor
endfor
% Design parameters of beams
BeamDesign{1,1} = 'Elevation';
BeamDesign{1,2} = 'Height';
BeamDesign{1,3} = 'F*_i';
BeamDesign{1,4} = 'beta_i';
BeamDesign{1,5} = 'M_pb_positive';
BeamDesign{1,6} = 'M_pb_negative';
for j = 1:6
  BeamDesign{2,j} = '-----';
endfor
for b = 1:B
  for j = 1:2
    for i = 1:NumberOfStory
    BeamDesign{2+(b-1)*(NumberOfStory+1)+i,j} = Parameters{i,j};
    endfor
  endfor
  for j = 1:B
    for i = 1:NumberOfStory
      BeamDesign{2+(b-1)*(NumberOfStory+1)+i,3} = ForceInBays{i,j};
    endfor
    for i = 1:NumberOfStory
      BeamDesign{2+(b-1)*(NumberOfStory+1)+i,4} = Parameters{i,6};
    endfor
    for i = 1:NumberOfStory
      BeamDesign{2+(b-1)*(NumberOfStory+1)+i,5} = M_pb_posInStories{i,j};
    endfor
    for i = 1:NumberOfStory
      BeamDesign{2+(b-1)*(NumberOfStory+1)+i,6} = M_pb_negInStories{i,j};
    endfor
  endfor
endfor
BeamDesign{size(BeamDesign,1)+1,1} = 'Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.';
%xlswrite(OutputName, BeamDesign, 'Design of beams');

%
%--------------------------------------------------
% Design of non-DYMs
%--------------------------------------------------

% Computing expected shear force at the desired beam plastic hinge locations
for j = 1:B
  for i = 1:NumberOfStory
    ShearBeamDownward{i,j} = ((abs(Ksi * M_pb_negInStories{i,j}) + abs(M_pb_posInStories{i,j})) / (0.9 * BayWidth(1,j))) + 0.5 * w_tributary(i,j) * 0.9 * BayWidth(1,j);
  endfor
  for i = 1:NumberOfStory
    ShearBeamUpward{i,j} = ((abs(Ksi * M_pb_negInStories{i,j}) + abs(M_pb_posInStories{i,j})) / (0.9 * BayWidth(1,j))) - 0.5 * w_tributary(i,j) * 0.9 * BayWidth(1,j);
  endfor
endfor
% Computing total summation of expected shear force of the beams
for j = 1:B
  SigmaShearBeamDownward{1,j} = 0;
  for i = 1:NumberOfStory
    SigmaShearBeamDownward{1,j} = SigmaShearBeamDownward{1,j} + ShearBeamDownward{i,j};
  endfor
endfor
for j = 1:B
  SigmaShearBeamUpward{1,j} = 0;
  for i = 1:NumberOfStory
    SigmaShearBeamUpward{1,j} = SigmaShearBeamUpward{1,j} + ShearBeamUpward{i,j};
  endfor
endfor
% Computing axial forces of columns
ColumnAxialForce{1,1} = 'Elevation';
ColumnAxialForce{1,2} = 'Height';
ColumnAxialForce{1,3} = 'P_c';
for j = 1:3
  ColumnAxialForce{2,j} = '-----';
endfor
for j = 1:2
  for i = 1:NumberOfStory
    ColumnAxialForce{2+i,j} = Parameters{i,j};
  endfor
endfor
for j = 1:B+1
  if j == 1
    for i = 1:NumberOfStory
      ColumnAxialForce{2+i,2+j} = ShearBeamDownward{i,j};
    endfor
  endif
  if (j > 1) && (j < B+1)
    ColumnAxialForce{1,2+j} = '';
    ColumnAxialForce{2,2+j} = '-----';
    for i = 1:NumberOfStory
      ColumnAxialForce{2+i,2+j} = ShearBeamDownward{i,j} - ShearBeamUpward{i,j-1};
    endfor
  endif
  if j == B+1
    ColumnAxialForce{1,2+j} = '';
    ColumnAxialForce{2,2+j} = '-----';
    for i = 1:NumberOfStory
      ColumnAxialForce{2+i,2+j} = ShearBeamDownward{i,j-1};
    endfor
  endif
endfor
ColumnAxialForce{2+NumberOfStory+1,1} = 'Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.';
for j = 2:size(ColumnAxialForce,2)
  ColumnAxialForce{2+NumberOfStory+1,j} = '';
endfor
%xlswrite(OutputName, ColumnAxialForce,'Axial force of columns');

% Computing alpha for exterior and interior columns tree
for i = 1:NumberOfStory
  Parameters{i,9} = Parameters{i,7} / SigmaBeta2;
endfor
% Computing lateral forces of columns
ColumnLataralForce{1,1} = 'Elevation';
ColumnLataralForce{1,2} = 'Height';
ColumnLataralForce{1,3} = 'F_L';
for j = 1:3
  ColumnLataralForce{2,j} = '-----';
endfor
for j = 1:B
  S = 0;
  for i = 1:NumberOfStory
    S = S + Ksi * M_pb_negInStories{i,j};
    endfor
  tempsigma1{1,j} = S;
endfor
for j = 1:B
  S = 0;
  for i = 1:NumberOfStory
    S = S + ShearBeamDownward{i,j} * 0.5 * (0.9 * BayWidth(1,j));
  endfor
  tempsigma2{1,j} = S;
endfor
tempsigma3 = 0;
for i = 1:NumberOfStory
  tempsigma3 = tempsigma3 + Parameters{i,9} * Parameters{i,2};
endfor
for j = 1:B
  S = 0;
  for i = 1:NumberOfStory
    S = S + abs(Ksi * M_pb_posInStories{i,j});
  endfor
  tempsigma4{1,j} = S;
endfor
for j = 1:B
  S = 0;
  for i = 1:NumberOfStory
    S = S + abs(Ksi * M_pb_negInStories{i,j});
  endfor
  tempsigma5{1,j} = S;
endfor
for j = 1:B
  S = 0;
  for i = 1:NumberOfStory
    S = S + ShearBeamUpward{i,j} * 0.5 * (0.9 * BayWidth(1,j));
  endfor
  tempsigma6{1,j} = S;
endfor
for j = 1:2
  for i = 1:NumberOfStory
    ColumnLataralForce{2+i,j} = Parameters{i,j};
  endfor
endfor
for j = 1:B+1
  if j == 1
    for i = 1:NumberOfStory
      ColumnLataralForce{2+i,2+j} = Parameters{i,9} * (tempsigma1{1,j} + tempsigma2{1,j} + M_pc) / tempsigma3;
    endfor
  endif
  if (j > 1) && (j < B+1)
    ColumnLataralForce{1,2+j} = '';
    ColumnLataralForce{2,2+j} = '-----';
    for i = 1:NumberOfStory
      ColumnLataralForce{2+i,2+j} = Parameters{i,9} * (tempsigma4{1,j-1} + tempsigma5{1,j} + tempsigma2{1,j} + tempsigma6{1,j-1} + 2 * M_pc) / tempsigma3;
    endfor
  endif
  if j == B+1
    ColumnLataralForce{1,2+j} = '';
    ColumnLataralForce{2,2+j} = '-----';
    for i = 1:NumberOfStory
      ColumnLataralForce{2+i,2+j} = Parameters{i,9} * (tempsigma1{1,j-1} + tempsigma2{1,j-1} + M_pc) / tempsigma3;
    endfor
  endif
endfor
ColumnLataralForce{2+NumberOfStory+1,1} = 'Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.';
for j = 2:size(ColumnLataralForce,2)
  ColumnLataralForce{2+NumberOfStory+1,j} = '';
endfor
%xlswrite(OutputName, ColumnLataralForce,'Lateral force of columns');

% Computing shear forces of columns
ColumnShearForce = ColumnLataralForce;
ColumnShearForce{1,3} = 'V_c';
for j = 3:size(ColumnLataralForce,2)
  for i = 4:size(ColumnLataralForce,1)
    ColumnShearForce{i,j} = ColumnShearForce{i-1,j} + ColumnLataralForce{i,j};
  endfor
endfor
ColumnShearForce{2+NumberOfStory+1,1} = 'Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.';
for j = 2:size(ColumnShearForce,2)
  ColumnShearForce{2+NumberOfStory+1,j} = '';
endfor
%xlswrite(OutputName, ColumnShearForce,'Shear force of columns');
clearvars -except Weight BaseShear LateralForces BeamDesign ColumnAxialForce ColumnLataralForce ColumnShearForce
%{
Copyright (c) 2021 Farshad Rasuli <farshad.rasuli@gmail.com> under GNU GPL v3.0. All rights reserved.
%}
