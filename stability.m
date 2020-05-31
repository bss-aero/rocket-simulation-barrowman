% C�LCULO DO CENTRO DE PRESS�O UTILIZANDO O M�TODO DE BARROWMAN E CUT OUT

function [Cp_Barrowman,Cn] = Estabilidade(Mach,Geometria_est,Geometria_est_bosters,Bosters)

if Bosters == 1
    [CP_bosters,CN_bosters] = Est(Mach,Geometria_est_bosters,1);
    [CP,CN] = Est(Mach,Geometria_est,0);
    Cp_Barrowman = (CP_bosters*CN_bosters + CP*CN)/(CN_bosters + CN);
    Cn = CN_bosters + CN;
else
    [Cp_Barrowman,Cn] = Est(Mach,Geometria_est,0);
end
    
    
    function [CP,CN] = Est(Mach,Geometria,Bosters)
    %% DADOS GEOM�TRICOS
    r = Geometria(1);
    C = Geometria(2);
    L_coifa = Geometria(3);
    cr = Geometria(4);
    ct = Geometria(5);
    N = Geometria(6);
    tau_c = Geometria(7);
    AR = Geometria(8);
    m = Geometria(9);
    s = Geometria(10);
    A_ref = Geometria(11);
    L_aleta = Geometria(12);

    ... COIFA

    % -> o formato da coifa pode ser padronizado ou � poss�vel entrar com a
    % equa��o do perfil para se calcular os coeficientes.

    r_base = r; % raio na base da coifa
    % C�nica - C = 1
    % Oviva - C = 2
    % Parab�lica - C = 3
    if C == 1    
        X_coifa = L_coifa*(2/3);
    elseif C == 2   % ogiva
        X_coifa = L_coifa*0.466;
    elseif C == 3
        X_coifa = L_coifa/2;
    end

    % -> A derivada do coeficiente de for�a normal em rela��o ao �ngulo de
    % ataque da coifa � considerado constante independente da geometria.
    C_coifa = 2;


    ... CONICAL SHOULDER AND CONICAL BOATTAIL
        % Conical shoulder: r1 < r2; Conical boattail: r1 > r2

    %    |          |
    %    |<---r1--->|   
    %    \          /   | L_cs ou L_cb
    %     \        /    |
    %      |<-r2->|
    %      |      |
    %      |      |

    % -> conical shoulder
    r1s = 0;
    r2s = 0;
    L_cs = 0;
    H_cs = 0; % dist�ncia da ponta do foguete at� o come�o do CS

    if r1s == r2s
        X_cs = 1;
        C_cs = 0;
    else
        X_cs = H_cs + L_cs*( 2*r2s^2 - r1s^2 - r1s*r2s )/(3*(r2s^2 - r1s^2));
        C_cs = (2/(r_base^2))*(r2s^2 - r1s^2);
    end

    % -> conical boattail
    r1b = 0;
    r2b = 0;
    L_cb = 0;
    H_cb = 0; % dist�ncia da ponta do foguete at� o come�o do CB

    if r1b == r2b
        X_cb = 1;
        C_cb = 0;
    else
        X_cb = H_cb + L_cb*( 2*r2b^2 - r1b^2 - r1b*r2b )/(3*(r2b^2 - r1b^2));
        C_cb = (2/(r_base^2))*(r2b^2 - r1b^2);
    end

    ... ALETAS


    %        
    %        |\        |
    %        | \       |
    %        |  \      | m
    %        |   \     |
    %      a |\   \    |
    %        | \   \   |
    %        |  \l  |
    %        |   \  |
    %         \   \ |         
    %          \   \|           
    %           \   | b         
    %            \  |          
    %             \ |            
    %              \|                
    %        --------
    %            s

    a = cr;
    b = ct;
    r_aletas = r_base; % raio do foguete onde as aletas est�o
    lambda = b/a;
    cm = (2/3)*a*(1 + lambda + lambda^2)/(1 + lambda);
    Area_t = (a + b)*s/2;

    H_aletas = L_aleta; % dist�ncia do ponto mais alto do foguete at� a o come�o das aletas
    X_aletas = H_aletas + m*(a + 2*b)/(3*(a + b)) + (1/6)*(a + b - a*b/(a + b));
    
     % Coeficiente que multiplica o coeficiente da derivada da aleta
     if Bosters == 0
         if (N == 3) || (N == 4)
             coef_aleta = 1;
         elseif N == 5
             coef_aleta = 0.948;
         elseif N == 6
             coef_aleta = 0.913;
         elseif N == 7
             coef_aleta = 0.854;
         elseif N == 8
             coef_aleta = 0.810;
         elseif N > 8
             coef_aleta = 0.750;
         end
         coef_aleta = coef_aleta*N/2;
     else
         if N == 3
             coef_aleta = sqrt(3);
         elseif N == 4
             coef_aleta = 3;
         elseif N == 2;
             coef_aleta = 2;
         elseif N == 1;
             coef_aleta = 1;
         end
     end

     if Bosters == 1
         coef_boster = coef_aleta;
         coef_aleta = 0;
     else
         coef_boster = 1;
     end

    % -> O c�lculo do coeficiente da derivada da for�a em rela��o ao �ngulo de
    % ataque depende do n�mero de aletas e de coeficientes de interfer�ncia que
    % ser�o calulados posteriormente.

    % C�LCULO DA POSI��O DO CENTRO DE PRESS�O PELO M�TODO DO BARROWMAN
    d = 2*r_base;

    beta = sqrt(1 - Mach^2);
    C_aleta = pi*AR*(Area_t/A_ref)/( 1 + sqrt(1 + ( 0.5*beta*AR/cos(tau_c) )^2 ) );
    kfi = 1 + r_aletas/(s + r_aletas); % fator de interfer�ncia entre as aletas e o corpo
    C_aletas = coef_aleta*kfi*C_aleta;
    CP = ( C_coifa*X_coifa + X_aletas*C_aletas + X_cs*C_cs + X_cb*C_cb )/( C_coifa + C_aletas + C_cs + C_cb );
    CN = coef_boster*(C_coifa + C_cs + C_cb + C_aletas);
    end

end