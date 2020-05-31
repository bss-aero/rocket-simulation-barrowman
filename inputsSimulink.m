
clear
close all
clc
 
%% FOGUETE

% TODAS AS DIMENSÕES ABAIXO ESTÃO EM METROS
% 1)  GEOMETRIA:
% 1.1) COIFA
L_coifa = 0.5;         % comprimento da coifa
r = 0.155/2;           % raio de base

% Geometria da coifa
% 1 - cone
% 2 - ogiva
% 3 - parabola
C = 2;

% 1.2) CORPO
L_corpo = 2;          % comprimento do corpo do foguete

% 1.3) ALETA
%        |\        |
%  F     | \       |
%  O     |  \      | m
%  G     |   \     |
%  U  cr |\   \    |
%  E     | \   \   |
%  T     |  \l  |
%  E     |   \  |
%         \   \ |         
%          \   \|           
%           \   | ct       
%            \  |          
%             \ |            
%              \|                
%        --------
%            s

L_aleta = 1.85; % distância da ponta do foguete até o comepo das aletas
ct = 0.25;
cr = 0.25;
m = 0.20;
s = 0.20;
l = sqrt(s^2 + (0.5*(ct - cr) + m)^2);
    
N = 4;                        % número de aletas
tr = 0.0025;                  % espessura máxima do perfil na raíz do perfil

% 2) MOTOR
% 2.2) Motot sólido
curvas = [
   0     0
   0.017 55.478
   0.059 157.11
   0.101 167.832
   0.143 174.359
   0.185 179.953
   0.228 184.149
   0.272 187.879
   0.312 188.345
   0.353 189.744
   0.399 190.676
   0.439 187.413
   0.484 189.744
   0.522 186.48
   0.566 185.548
   0.606 181.818
   0.651 177.156
   0.694 170.629
   0.736 167.366
   0.777 158.974
   0.82 156.177
   0.863 154.312
   0.905 136.597
   0.947 74.592
   0.993 20.979
   1.036 4.662
   1.076 0.0
];

Empuxo_boster = 4*curvas(:,2);             
tempo_boster = curvas(:,1);
tempo_queima_bosters = max(tempo_boster);

Massa_motor_solido = [3.5 1]; % valores antes e depois da queima
tempo_motor_solido = [0 tempo_queima_bosters];

% 3) QUEDA
% 3.1) PARAQUEDAS
% Coeficiente de arrasto
Cd_paraquedas1 = 1.55;
Cd_paraquedas2 = 2.2;

Cd_lateral1 = 1.5;
Cd_lateral2 = 1.8;

% Areas
A_paraquedas1 = 0.3; %[m^2]
A_paraquedas2 = 4.15; %[m^2]
A_lateral1 = 1.5;
A_lateral2 = 2.5;

% ABERTURA DOS PARAQUEDAS
% O primeiro paraquedas será aberto quando a velocidade vertical for
% negativa e o foguete chegar a 0.95% do apogeu;
% O segundo paraquedas será aberto quando o foguete chegar a uma
% determinada altura.
tempo_abertura1 = 5; % tempo gasto para abrir o primeiro paraquedas
tempo_abertura2 = 7; % tempo gasto para abrir o segundo paraquedas

%% Massas e momentos de inércia
CG = 1.7; 	% Posição do cg em relação ao nariz do foguete
Massa_foguete = 1; % [kg]

%% Vento
Vento_eixo_x = 1;
Vento_eixo_y = 0;

%% LANÇAMENTO
phi = 0;	% Ângulo entre o eixo longitudinal e a vertical
theta = 0;	% Ângulo entre a projeção do foguete no plano horizontal e o eixo x

% Comprimento da guia do foguete
dist_guia = 12; % [m]

% Altura do lançamento acima do nível do mar
h0 = 1200;

%% AMBIENTE
g0 = 9.82;		
RE = 6371000;

%% Vari�veis 
% Modelo_arrasto_Barrowman
L_total = L_coifa + L_corpo;
Rs = 1e-5;
At = (cr + ct)*s/2;
A_ref = pi*r^2;
Awb = 2*r*L_corpo;
fn = L_coifa/(2*r);
tau_L = atan(m/s);
tau_C = acos(s/l);
AR = 2*(s^2)/At;

Geometria = [L_total,Rs,N,At,A_ref,Awb,tr,tau_L,cr,AR,tau_C,fn,r,s,C];

% Estabilidade
Geometria_est = [r,C,L_coifa,cr,ct,N,tau_C,AR,m,s,A_ref,L_aleta];

% Queda
% Correçãoo do coeficiente de arrasto
Cd_paraquedas1 = Cd_paraquedas1*A_paraquedas1/A_ref;
Cd_paraquedas2 = Cd_paraquedas2*A_paraquedas2/A_ref;
Cd_lateral1 = Cd_lateral1*A_lateral1/A_ref;
Cd_lateral2 = Cd_lateral2*A_lateral2/A_ref;

% Estabilidade
Mach_ref = 0.4;
[Cp_Barrowman_SB,Cn_SB] = stability(Mach_ref,Geometria_est,[],0); % sem os bosters
Cn = Cn_SB;
Cp = Cp_Barrowman_SB;

% Margem est�tica
ME = (Cp - CG)/(2*r);
