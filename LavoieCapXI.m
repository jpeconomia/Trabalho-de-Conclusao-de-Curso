%Resolução do Modelo GROWTH do capítulo 11 do livro Monetary Economics de Godley e Lavoie (2006), pelo MATLAB para o meu Trabalho de Conclusão de Curso

%GODLEY, W.; LAVOIE, M. Monetary Economics: An Integrated Approach to Credit, Money, Income, Production and Wealth. United Kingdom: Palgrave Macmillan. 2006.

clear;
clc;

T = 1000;
Z = 500;


%Parâmetros
addBL           = 0.02;              %mark-up dos bills
alpha1          = 0.75;              %propensão marginal a consumir a partir de renda disponível esperada e empréstimos
alpha2          = 0.064;             %propensão marginal a consumir da renda passada
bandB           = 0.01;              %poder de barganha
bandT           = 0.01;              %poder de barganha
beta            = 0.5;               %vendas esperadas entre 0 e 1
beta_b          = 0.4;               %capital próprio dos bancos
delta           = 0.10667;           %depreciação ???
delta_rep       = 0.1;               %repagamento de empréstimos
epsilon         = 0.8;               %ajuste do mark-up
epsilon_2       = 0.5;               %renda disponível esperada 
epsilon_b       = 0.25;              %calotes esperados    
eta_0           = 0.07416;           %empréstimos
eta_r           = 0.4;               %empréstimos
meta            = 0.6;               %ajuste para nível de emprego desejado
NCAR            = 0.1;               %regulação de alavancagem
gama            = 0.15;              %inventários esperados entre 0 e 1, de curto prazo
gama_u          = 0.05;              %impacto de utilização de capacidade na taxa de crescimento do capital
gama_r          = 0.1;               %impacto do custo de empréstimos (juros) na taxa de crescimento do capital
grg             = 0.03;              %Tx de crescimento do gasto do governo
grp             = 0.03;              %produtividade
gr0             = 0.00122;           %exógeno na taxa de crescimento do capital
ho              = 0.05;              %requisito de depósitos
lambda_b        = 0.0153;            %dividendos dos bancos
lambda_c        = 0.05;              %demanda transacional por moeda
lambda_10       = 0.11909;           %alocação de portfólio
lambda_11       = 6;
lambda_12       = 2;
lambda_13       = 2;
lambda_14       = 2;
lambda_15       = 0.3;
lambda_20       = 0.25;              
lambda_21       = 2;
lambda_22       = 6;
lambda_23       = 2;
lambda_24       = 2;
lambda_25       = 0.1;
lambda_30       = -0.04341;
lambda_31       = 2;
lambda_32       = 2;
lambda_33       = 6;
lambda_34       = 2;
lambda_35       = 0.1;
lambda_40       = 0.67432;
lambda_41       = 2;
lambda_42       = 2;
lambda_43       = 2;
lambda_44       = 6;
lambda_45       = 0.1;
omega_0         = -0.23322;          %intercepto do salário real desejado
omega_1         = 1;                 %reação do salário real desejado à produtividade
omega_2         = 2;                 %demanda salarial por poder de barganha
omega_3         = 0.45621;           %salário nominal
psi_d           = 0.15255;           %lucros distribuídos
psi_u           = 0.92;              %lucros retidos
sigma_N         = 0.1666;            %custo normal unitário histórico
sigma_T         = 0.2;               %target dos inventários no longo prazo
theta           = 0.22844;           %alíquota de imposto 
zeta_m          = 0.0008;            %variação do retorno dos depósitos



%Variáveis Exógenas 

rb          = 0.035.*ones(T,1);
npl         = 0.02.*ones(T,1);
N_fe        = 86.*ones(T,1); 
top         = 0.12.*ones(T,1);
bot         = 0.05.*ones(T,1);
z1          = ones(T,1);
z2          = ones(T,1);
z3          = ones(T,1);
z4          = ones(T,1);
z5          = ones(T,1);
%Pré-alocação de Variáveis Endógenas


y                 = 83867230.*ones(T,Z);                                    %decisão de produção real
p                 = 7.1537.*ones(T,Z);                                      %nível de preços
addl              = 0.04592.*ones(T,Z);                                     %spread entre empréstimos e depósitos
B_bd              = 42500560.*ones(T,Z);                                    %bills demandados pelos bancos
B_bs              = 42500560.*ones(T,Z);                                    %bills ofertados aos bancos
B_cbd             = 45083930.*ones(T,Z);                                    %bills demandados banco central
B_cbs             = 45083930.*ones(T,Z);                                    %bills ofertados banco central
B_hd              = 32340300.*ones(T,Z);                                    %bills demandados famílias
B_hs              = 32340300.*ones(T,Z);                                    %bills ofertados famílias
B_s               = (B_bs(1,1)+B_cbs(1,1)+B_hs(1,1)).*ones(T,Z);            %oferta total de bills
BL_d              = 8141406.*ones(T,Z);                                     %bonds demandados
BL_s              = 8141406.*ones(T,Z);                                     %bonds ofertados
BLR               = 0.10834.*ones(T,Z);                                     %liquidez dos bancos
BUR               = 0.06324.*ones(T,Z);                                     %peso da dívida das famílias
C                 = 7120643.*ones(T,Z);                                     %consumo nominal
c                 = 7120643.*ones(T,Z);                                     %consumo real
CAR               = 0.09245.*ones(T,Z);                                     %capital adequacy ratio
CG                = ones(T,Z);                                              %ganhos de capital dos bonds
ER                = ones(T,Z);                                              %taxa de emprego
e_d               = 5039.6.*ones(T,Z);                                      %ações demandadas
e_s               = 5039.6.*ones(T,Z);                                      %ações ofertadas
eta               = 0.0498.*ones(T,Z);                                      %razão entre empréstimos novos e renda
F_b               = 1688936.*ones(T,Z);                                     %lucros bancários
F_bt              = 1688947.*ones(T,Z);                                     %lucros bancários desejados
F_cb              = (rb(1,1).*B_cbd(1,1)).*ones(T,Z);                       %lucros do banco central
F_f               = 17508820.*ones(T,Z);                                    %lucros das firmas 
F_ft              = 17443490.*ones(T,Z);                                    %lucros desejados das firmas
FD_b              = 1283158.*ones(T,Z);                                     %dividendos dos bancos
FD_f              = 2586436.*ones(T,Z);                                     %dividendos das firmas
FU_b              = 405778.*ones(T,Z);                                      %lucros retidos dos bancos
FU_bt             = ones(T,Z);                                              %lucros retidos desejados dos bancos
FU_f              = 14674200.*ones(T,Z);                                    %lucros retidos das firmas
FU_ft             = 14589340.*ones(T,Z);                                    %lucros retidos desejados das firmas
G                 = (22681210.*7.1537).*ones(T,Z);                          %gasto nominal do governo
GL                = ones(T,Z);                                              %novos empréstimos desejados
g                 = 22681210.*ones(T,Z);                                    %gasto real do governo
gr_k              = 0.0300.*ones(T,Z);                                      %TX do crescimento de capital
H_bd              = 1961455.*ones(T,Z);                                     %moeda demandada pelos bancos
H_bs              = 1961455.*ones(T,Z);                                     %moeda ofertada aos bancos
H_hd              = 2546938.*ones(T,Z);                                     %moeda desejada pelas famílias
H_hs              = 2546938.*ones(T,Z);                                     %moeda ofertada às famílias
H_s               = (H_hs(1,1)+H_bs(1,1)).*ones(T,Z);                       %base monetária
HC_e              = ones(T,Z);                                              %custos esperados históricos
IN                = 11218840.*ones(T,Z);                                    %inventários nominais
i                 = 16376300.*ones(T,Z);                                    %investimento real bruto
I                 = (i(1,1).*p(1,1)).*ones(T,Z);                            %investimento bruto nominal
in                = 2004754.*ones(T,Z);                                     %inventários
in_e              = 2004754.*ones(T,Z);                                     %inventários esperados
in_T              = 2004754.*ones(T,Z);                                     %inventários desejados
k                 = 1725114.*ones(T,Z);                                     %estoque de capital real
K                 = (k(1,1).*p(1,1)).*ones(T,Z);                            %estoque nominal de capital
L_fd              = 15457450.*ones(T,Z);                                    %empréstimos desejados pelas firmas
L_fs              = 15457450.*ones(T,Z);                                    %empréstimos ofertados às firmas
L_hd              = 20923050.*ones(T,Z);                                    %empréstimos desejados pelas famílias
L_hs              = 20923050.*ones(T,Z);                                    %empréstimos ofertados pelas famílias
M_d               = 39229100.*ones(T,Z);                                    %depósitos demandados
M_h               = 39229100.*ones(T,Z);                                    %depósitos das famílias
M_s               = 39229100.*ones(T,Z);                                    %depósitos ofertados
N                 = 87.181.*ones(T,Z);                                      %nível de emprego
N_T               = 87.181.*ones(T,Z);                                      %nível de emprego desejado
NPL               = 299377.*ones(T,Z);                                      %non performing loans
NUC               = 5.5961.*ones(T,Z);                                      %custo unitário normal
NHUC              = ones(T,Z);                                              %custo unitário normal histórico
nl                = 661949.*ones(T,Z);                                      %empréstimos novos reais
NL                = (nl(1,1).*p(1,1)).*ones(T,Z);                           %empréstimos novos para as famílias
npl_e             = 0.02.*ones(T,Z);                                        %empréstimos não pagos esperados
omega_T           = 109566.*ones(T,Z);                                      %salário real desejado
OF_b              = 3363394.*ones(T,Z);                                     %fundos próprios dos bancos
OF_be             = 3363394.*ones(T,Z);                                     %fundos próprios esperados
OF_bt             = 3363394.*ones(T,Z);                                     %fundos próprios desejados                                                                                          
pe                = 17621.*ones(T,Z);                                       %preços das ações
PE                = ((pe(1,1).*e_s(1,1))./F_f(1,1))./ones(T,Z);             %razão preços e rendimentos
phi               = 0.26416.*ones(T,Z);                                     %mark-up
phi_T             = 0.26417.*ones(T,Z);                                     %mark-up desejado
pi                = 0.0026.*ones(T,Z);                                      %taxa de inflação
pr                = 134621.*ones(T,Z);                                      %produtividade
q                 = ones(T,Z);                                              %q de tobin
REP               = 2026110.*ones(T,Z);                                     %empréstimos pagos
rbl               = 0.055.*ones(T,Z);                                       %taxa de juros de longo prazo
pbl               = (1./rbl(1,1)).*ones(T,Z);                               %preço dos bonds                                              
rl                = 0.06522.*ones(T,Z);                                     %taxa de juros dos empréstimos
rm                = 0.02.*ones(T,Z);                                        %taxa de juros dos depósitos
rrl               = (((1+rl(1,1))./(1+pi(1,1))) - 1).*ones(T,Z);            %taxa de juros real dos empréstimos
rk                = 0.03008.*ones(T,Z);                                     %retorno sobre dividendos
s                 = 11677980.*ones(T,Z);                                    %vendas reais
S                 = (s(1,1).*p(1,1)).*ones(T,Z);                            %vendas nominais
s_e               = 11677990.*ones(T,Z);                                    %vendas esperadas
sigma_se          = (in(1,1)./s_e(1,1)).*ones(T,Z);                         %razão inventários/vendas esperadas
V                 = 156798400.*ones(T,Z);                                   %riqueza
V_fma             = (V(1,1)-H_hd(1,1)-OF_b(1,1)).*ones(T,Z);                %riqueza disponível para investir
v                 = (V(1,1)./p(1,1)).*ones(T,Z);                            %riqueza real
W                 = 753353.*ones(T,Z);                                      %salário nominal
WB                = (W(1,1).*N(1,1)).*ones(T,Z);                            %folha de pagamentos
Y                 = (y(1,1).*p(1,1)).*ones(T,Z);                            %PIB nominal
YD_r              = 54660390.*ones(T,Z);                                    %renda disponível regular
YD_hs             = ones(T,Z);                                              %renda disponível haig-simons
YP                = 70844000.*ones(T,Z);                                    %renda pessoal nominal
Tx                = (theta.*YP(1,1)).*ones(T,Z);                            %tributação
yd_r              = (YD_r(1,1)./p(1,1)).*ones(T,Z);                         %renda disponível regular real
yd_er             = 7585732.*ones(T,Z);                                     %renda disponível regular real esperada
GD                = (B_bs(1,1)+B_hs(1,1)+...
    BL_s(1,1).*pbl(1,1)+H_s(1,1)).*ones(T,Z);                               %dívida nominal do governo
u                 = (k(1,1)./y(1,1)).*ones(T,Z);                            %utilização de capacidade
UC                = (WB(1,1)./y(1,1)).*ones(T,Z);                           %custo unitário
PSBR              = (G(1,1) + rb(1,1).*(B_hs(1,1) + B_bs(1,1))...
    + BL_s(1,1) - Tx(1,1)).*ones(T,Z);                                      %déficit nominal do governo


%Solução do Modelo :

for t=2:T
    for z=2:Z
        if t==810
    %   omega_1 = 0.9;
    %   psi_u   = 0.91;
    %   psi_d   = 0.16;
    %   eta_0   = 0.08;
    %    rb(t,1) = 0.045;
        end
        if z==2 
addl(t,z-1)       = addl(t-1,end);
B_bd(t,z-1)       = B_bd(t-1,end);
B_bs(t,z-1)       = B_bs(t-1,end);
B_cbd(t,z-1)      = B_cbd(t-1,end);
B_cbs(t,z-1)      = B_cbs(t-1,end);
B_hd(t,z-1)       = B_hd(t-1,end);
B_hs(t,z-1)       = B_hs(t-1,end);
BL_d(t,z-1)       = BL_d(t-1,end);
BL_s(t,z-1)       = BL_s(t-1,end);
BLR(t,z-1)        = BLR(t-1,end);
B_s(t,z-1)        = B_s(t-1,end);
BUR(t,z-1)        = BUR(t-1,end);
C(t,z-1)          = C(t-1,end);
c(t,z-1)          = c(t-1,end);
CAR(t,z-1)        = CAR(t-1,end);
CG(t,z-1)         = CG(t-1,end);
ER(t,z-1)         = ER(t-1,end);
e_d(t,z-1)        = e_d(t-1,end);
e_s(t,z-1)        = e_s(t-1,end);
eta(t,z-1)        = eta(t-1,end);
F_b(t,z-1)        = F_b(t-1,end);
F_bt(t,z-1)       = F_bt(t-1,end);
F_cb(t,z-1)       = F_cb(t-1,end);
F_f(t,z-1)        = F_f(t-1,end);
F_ft(t,z-1)       = F_ft(t-1,end);
FD_b(t,z-1)       = FD_b(t-1,end);
FD_f(t,z-1)       = FD_f(t-1,end);
FU_b(t,z-1)       = FU_b(t-1,end);
FU_bt(t,z-1)      = FU_bt(t-1,end);
FU_f(t,z-1)       = FU_f(t-1,end);
FU_ft(t,z-1)      = FU_ft(t-1,end);
G(t,z-1)          = G(t-1,end);
GD(t,z-1)         = GD(t-1,end);
GL(t,z-1)         = GL(t-1,end);
g(t,z-1)          = g(t-1,end);
gr_k(t,z-1)       = gr_k(t-1,end);
H_s(t,z-1)        = H_s(t-1,end);
H_bd(t,z-1)       = H_bd(t-1,end);
H_bs(t,z-1)       = H_bs(t-1,end);
H_hd(t,z-1)       = H_hd(t-1,end);
H_hs(t,z-1)       = H_hs(t-1,end);
HC_e(t,z-1)       = HC_e(t-1,end);
I(t,z-1)          = I(t-1,end);
IN(t,z-1)         = IN(t-1,end);
in(t,z-1)         = in(t-1,end);
in_e(t,z-1)       = in_e(t-1,end);
in_T(t,z-1)       = in_T(t-1,end);
K(t,z-1)          = K(t-1,end);
k(t,z-1)          = k(t-1,end);
L_fd(t,z-1)       = L_fd(t-1,end);
L_fs(t,z-1)       = L_fs(t-1,end);
L_hd(t,z-1)       = L_hd(t-1,end);
L_hs(t,z-1)       = L_hs(t-1,end);
M_d(t,z-1)        = M_d(t-1,end);
M_h(t,z-1)        = M_h(t-1,end);
M_s(t,z-1)        = M_s(t-1,end);
N(t,z-1)          = N(t-1,end);
N_T(t,z-1)        = N_T(t-1,end);
NL(t,z-1)         = NL(t-1,end);
NPL(t,z-1)        = NPL(t-1,end);
NUC(t,z-1)        = NUC(t-1,end);
NHUC(t,z-1)       = NHUC(t-1,end);
nl(t,z-1)         = nl(t-1,end);
npl_e(t,z-1)      = npl_e(t-1,end); 
omega_T(t,z-1)    = omega_T(t-1,end);
OF_b(t,z-1)       = OF_b(t-1,end);
OF_be(t,z-1)      = OF_be(t-1,end);
OF_bt(t,z-1)      = OF_bt(t-1,end);
PE(t,z-1)         = PE(t-1,end);
PSBR(t,z-1)       = PSBR(t-1,end);
p(t,z-1)          = p(t-1,end);
pbl(t,z-1)        = pbl(t-1,end);
pe(t,z-1)         = pe(t-1,end);
phi(t,z-1)        = phi(t-1,end);
phi_T(t,z-1)      = phi_T(t-1,end);
pi(t,z-1)         = pi(t-1,end);
pr(t,z-1)         = pr(t-1,end);
q(t,z-1)          = q(t-1,end);
REP(t,z-1)        = REP(t-1,end);
rbl(t,z-1)        = rbl(t-1,end);
rl(t,z-1)         = rl(t-1,end);
rm(t,z-1)         = rm(t-1,end);
rrl(t,z-1)        = rrl(t-1,end);
rk(t,z-1)         = rk(t-1,end);
S(t,z-1)          = S(t-1,end);
s(t,z-1)          = s(t-1,end);
s_e(t,z-1)        = s_e(t-1,end);
sigma_se(t,z-1)   = sigma_se(t-1,end);
Tx(t,z-1)         = Tx(t-1,end);
u(t,z-1)          = u(t-1,end);
UC(t,z-1)         = UC(t-1,end);
V(t,z-1)          = V(t-1,end);
V_fma(t,z-1)      = V_fma(t-1,end);
v(t,z-1)          = v(t-1,end);
W(t,z-1)          = W(t-1,end);
WB(t,z-1)         = WB(t-1,end);
Y(t,z-1)          = Y(t-1,end);
YD_r(t,z-1)       = YD_r(t-1,end);
YD_hs(t,z-1)      = YD_hs(t-1,end);
YP(t,z-1)         = YP(t-1,end);
yd_r(t,z-1)       = yd_r(t-1,end);
yd_er(t,z-1)      = yd_er(t-1,end);
y(t,z-1)          = y(t-1,end);

        end

% 1) FIRMAS:
 
y(t,z)      = s_e(t,z-1) + in_e(t,z-1) - in(t-1,end);                       %11.1 poderia ser as vendas de fato que ainda seria um modelo keynesiano sem problemas
s_e(t,z)    = beta.*s(t,z-1) + (1-beta).*s(t-1,end).*(1+grp);               %11.2 colocar vendas presentes pra formular vendas esperadas
in_T(t,z)   = sigma_T.*s_e(t,z-1);                                          %11.3 
in_e(t,z)   = in(t-1,end) + gama.*(in_T(t,z-1) - in(t-1,end));              %11.4
in(t,z)     = in(t-1,end) + y(t,z-1) - s(t,z-1);                            %11.5
k(t,z)      = k(t-1,end).*(1 + gr_k(t,z-1));                                %11.6 
gr_k(t,z)   = gr0 + gama_u.*u(t,z-1) - gama_r.*rrl(t,z-1);                  %11.7
u(t,z)      = y(t,z-1)./k(t-1,end);                                         %11.8 utilização de capacidade vezes produtividade do capital (por que os dois não são do mesmo período?)
rrl(t,z)    = ((1+rl(t,z-1))./(1+pi(t,z-1))) - 1;                           %11.9
pi(t,z)     = (p(t,z-1)./p(t-1,end))-1;                                     %11.10
i(t,z)      = (gr_k(t,z-1) + delta).*k(t-1,end);                            %11.11
s(t,z)      = c(t,z-1) + g(t,z-1) + i(t,z-1);                               %11.12
S(t,z)      = s(t,z-1).*p(t,z-1);                                           %11.13
IN(t,z)     = in(t,z-1).*p(t,z-1);                                          %11.14
I(t,z)      = i(t,z-1).*p(t,z-1);                                           %11.15
K(t,z)      = k(t,z-1).*p(t,z-1);                                           %11.16
Y(t,z)      = s_e(t,z-1).*p(t,z-1) + (in_e(t,z-1) - in(t-1,end)).*p(t,z-1); %11.17 não tem equivalência entre essa definição de produto nominal e a do produto real

% 1.2) Decisões dos custos 

if ge(ER(t,z),1 - bandB) && ge(1 + bandT,ER(t,z))
    z3(t,1)      = 1;
else 
    z3(t,1)      = 0;
end
if ER(t,z)       > 1 + bandT 
    z4(t,1)      = 1; 
else
    z4(t,1)      = 0; 
end
if ER(t,z)       < 1 - bandB 
    z5(t,1)      = 1; 
else
    z5(t,1)      = 0; 
end


omega_T(t,z)= omega_0 + omega_1.*pr(t,z-1) -...
    omega_2.*(ER(t,z-1) + z3(t,1).*(1-ER(t,z-1))...
    - z4(t,1).*bandT + z5(t,1).*bandB);                                     %11.18
ER(t,z)     = N(t,z-1)./N_fe(t-1,end);                                      %11.19
W(t,z)      = W(t-1,end) + omega_3.*(omega_T(t,z-1).*p(t-1,end)...
    - W(t-1,end));                                                          %11.21
pr(t,z)     = pr(t-1,end).*(1 + grp);                                       %11.22
N_T(t,z)    = y(t,z-1)./pr(t,z-1);                                          %11.23
N(t,z)      = N(t-1,end) + meta.*(N_T(t,z-1) - N(t-1,end));                 %11.24 
WB(t,z)     = W(t,z-1).*N(t,z-1);                                           %11.25
UC(t,z)     = WB(t,z-1)./y(t,z-1);                                          %11.26 
NUC(t,z)    = W(t,z-1)./pr(t,z-1);                                          %11.27
NHUC(t,z)   = (1 - sigma_N).*NUC(t,z-1) +...
    sigma_N.*NUC(t-1,end).*(1 + rl(t,z-1));                                 %11.28

% 1.3) Decisões de preços 

p(t,z)      = (1 + phi(t,z-1)).*NHUC(t,z-1);                                %11.29
phi(t,z)    = phi(t-1,end) + epsilon.*(phi_T(t-1,end) - phi(t-1,end));      %11.30 deveria ser o alvo no presente
phi_T(t,z)  = F_ft(t,z-1)./HC_e(t,z-1);                                     %11.31
HC_e(t,z)   = (1 - sigma_se(t,z-1)).*s_e(t,z-1).*UC(t,z-1) + ...
    sigma_se(t,z-1).*s_e(t,z-1).*UC(t-1,end);                               %11.32
sigma_se(t,z)= in(t-1,end)./s_e(t,z-1);                                     %11.33
F_ft(t,z)   = FU_ft(t,z-1) + FD_f(t,z-1) + ...
    rl(t-1,end).*(L_fd(t,z-1) - IN(t-1,end));                               %11.34
FU_ft(t,z)  = psi_u.*I(t-1,end);                                            %11.35
FD_f(t,z)   = psi_d.*F_f(t-1,end);                                          %11.36
F_f(t,z)    = S(t,z-1) - WB(t,z-1) + (IN(t,z-1) - IN(t-1,end))...
    - rl(t,z-1).*IN(t-1,end);                                               %11.37 porque variação de inventário entra como positivo
FU_f(t,z)   = F_f(t,z-1) - FD_f(t,z-1) -...
    rl(t,z-1).*(L_fd(t-1,end) - IN(t-1,end) + NPL(t,z-1));                  %11.38
L_fd(t,z)   = L_fd(t-1,end) + I(t,z-1) + (IN(t,z-1) - IN(t-1,end)) ...
    - FU_f(t,z-1) - (e_s(t,z-1) - e_s(t-1,end)).*pe(t,z-1) - NPL(t,z-1);    %11.39
NPL(t,z)    = npl(t,1).*L_fd(t-1,end);                                      %11.40
e_s(t,z)    = e_s(t-1,end) + (1 - psi_u).*I(t-1,end)./pe(t,z-1);            %11.41
rk(t,z)     = FD_f(t,z-1)./(e_s(t-1,end).*pe(t-1,end));                     %11.42
PE(t,z)     = pe(t,z-1)./(F_f(t,z-1)./e_s(t-1,end));                        %11.43 obs: o preço das ações fica parado o tempo todo
q(t,z)      = (e_s(t,z-1).*pe(t,z-1) + L_fd(t,z-1))./(K(t,z-1) + IN(t,z-1));%11.44 não deveriam ser somente os empréstimos das firmas?

% 2) FAMÍLIAS
    
% 2.1) Renda e decisões de consumo


YP(t,z)     = WB(t,z-1) + FD_f(t,z-1) + FD_b(t,z-1)...
    + rm(t-1,end).*M_h(t-1,end) + rb(t-1,1).*B_hd(t-1,end)...
    + BL_d(t-1,end);                                                        %11.45
Tx(t,z)     = theta.*YP(t,z-1);                                             %11.46
YD_r(t,z)   = YP(t,z-1) - Tx(t,z-1) - rl(t-1,end).*L_hd(t-1,end);           %11.47
YD_hs(t,z)  = YD_r(t,z-1) + CG(t,z-1);                                      %11.48
CG(t,z)     = (pbl(t,z-1) - pbl(t-1,end)).*BL_d(t-1,end)...
    + (pe(t,z-1) - pe(t-1,end)).*e_d(t-1,end) +...
    (OF_b(t,z-1) - OF_b(t-1,end));                                          %11.49 nao deveriam ter a mudança de quantidades também?
V(t,z)      = V(t-1,end) + YD_r(t,z-1) - C(t,z-1);                          %11.50 desse jeito os ganhos de capital não afetam a riqueza
v(t,z)      = V(t,z-1)./p(t,z-1);                                           %11.51
C(t,z)      = p(t,z-1).*c(t,z-1);                                           %11.52
c(t,z)      = alpha1.*(yd_r(t,z-1) + nl(t,z-1)) + alpha2.*v(t-1,end);       %11.53
yd_er(t,z)  = epsilon_2.*yd_r(t,z-1) +...
    (1 - epsilon_2).*yd_r(t-1,end).*(1 + grp);                              %11.54 aqui também, renda presente definindo a renda esperada
yd_r(t,z)   = YD_r(t,z)./p(t,z-1) - pi(t,z-1).*V(t-1,end)./p(t,z-1);        %11.55

% 2.2) Decisões de Empréstimos

GL(t,z)     = eta(t,z-1).*YP(t,z-1);                                        %11.56
eta(t,z)    = eta_0 - eta_r.*rrl(t,z-1);                                    %11.57
NL(t,z)     = GL(t,z-1) - REP(t,z-1);                                       %11.58
REP(t,z)    = delta_rep.*L_hd(t-1,end);                                     %11.59
L_hd(t,z)   = L_hd(t-1,end) + NL(t,z-1);                                    %11.60
nl(t,z)     = NL(t,z-1)./p(t,z-1);                                          %11.61
BUR(t,z)    = (REP(t,z-1) + rl(t-1,end).*L_hd(t-1,end))./YP(t,z-1);         %11.62

% 2.3) Alocação de Portfólio 

M_d(t,z)    = (lambda_10 + lambda_11.*rm(t,z-1) - lambda_12.*rb(t,1) ...
    - lambda_13.*rbl(t,z-1) - lambda_14.*rk(t,z-1) ...
    + lambda_15.*(YP(t,z-1)./V_fma(t-1,end))).*V_fma(t-1,end);              %11.63 qual o sentido de deixar essa equação aqui se ela é a variável buffer
B_hd(t,z)   = (lambda_20 - lambda_21.*rm(t,z-1) + lambda_22.*rb(t,1) ...    % quando roda o modelo, aliás, esse aqui fica negativo
    - lambda_23.*rbl(t,z-1) - lambda_24.*rk(t,z-1) ...
    - lambda_25.*(YP(t,z-1)./V_fma(t-1,end))).*V_fma(t-1,end);              %11.64
BL_d(t,z)   = (lambda_30 - lambda_31.*rm(t,z-1) - lambda_32.*rb(t,1) ...
    + lambda_33.*rbl(t,z-1) - lambda_34.*rk(t,z-1) ...
    - lambda_35.*(YP(t,z-1)./V_fma(t-1,end))).*V_fma(t-1,end)./pbl(t,z);    %11.65
e_d(t,z)    = (lambda_40 - lambda_41.*rm(t,z-1) - lambda_42.*rb(t,1) ...
    - lambda_43.*rbl(t,z-1) + lambda_44.*rk(t,z-1) ...
    - lambda_45.*(YP(t,z-1)./V_fma(t-1,end))).*V_fma(t-1,end)./pe(t,z-1);   %11.66
M_h(t,z)    = V_fma(t,z-1) - B_hd(t,z-1) - pbl(t,z-1).*BL_d(t,z-1)...
    - pe(t,z-1).*e_d(t,z-1);                                                %11.67
V_fma(t,z)  = V(t,z-1) + L_hd(t,z-1) - H_hd(t,z-1) - OF_b(t,z-1);           %11.68 não deveria se definir a riqueza disponível para empréstimos somente a partir da riqueza líquida?
H_hd(t,z)   = lambda_c.*C(t,z-1);                                           %11.69
e_d(t,z)    = e_s(t,z-1);                                                   %11.70

% 3) SETOR PÚBLICO

% 3.1) Governo 

G(t,z)      = p(t,z-1).*g(t,z-1);                                           %11.71
g(t,z)      = g(t-1,end).*(1 + grg);                                        %11.72
PSBR(t,z)   = G(t,z-1) + rb(t-1,1).*(B_hs(t-1,end) + B_bs(t-1,end))...
    + BL_s(t-1,end) - Tx(t,z-1);                                            %11.73 por que não entram os títulos do banco central nessa definição?
B_s(t,z)    = B_s(t-1,end) + PSBR(t,z-1) - (BL_s(t,z-1)...
    - BL_s(t-1,end)).*pbl(t,z-1);                                           %11.74
GD(t,z)     = B_hs(t,z-1) + BL_s(t,z-1) + H_s(t,z-1);                       %11.75

% 3.2) Banco Central 

F_cb(t,z)   = rb(t-1,1).*B_cbd(t-1,end);                                    %11.76
BL_s(t,z)   = BL_d(t,z-1);                                                  %11.77
B_hs(t,z)   = B_hd(t,z-1);                                                  %11.78
H_hd(t,z)   = H_hs(t,z-1);                                                  %11.79
H_bd(t,z)   = H_bs(t,z-1);                                                  %11.80
H_s(t,z)    = H_bs(t,z-1) + H_hs(t,z-1);                                    %11.81
B_cbd(t,z)  = H_s(t,z-1);                                                   %11.82
B_cbs(t,z)  = B_cbd(t,z-1);                                                 %11.83
rbl(t,z)    = rb(t,1) + addBL;                                              %11.85
pbl(t,z)    = 1./rbl(t,z-1);                                                %11.86

% 4) SETOR BANCÁRIO 

% 4.1) Agregados Monetários

M_s(t,z)    = M_d(t,z-1);                                                   %11.87 porque o Ms nessas equações se ele é diferente da quantidade de moeda em posse das famílias?
L_fs(t,z)   = L_fd(t,z-1);                                                  %11.88
L_hs(t,z)   = L_hd(t,z-1);                                                  %11.89
H_bd(t,z)   = ho.*M_h(t,z-1);                                               %11.90
B_bs(t,z)   = B_s(t,z-1) - B_hs(t,z-1) - B_cbs(t,z-1);                      %11.91
B_bd(t,z)   = M_d(t,z-1) - L_fs(t,z-1) - L_hs(t,z-1)...
    - H_bd(t,z-1) + OF_b(t,z-1);                                            %11.92 parece que o sentido tá trocado nisso aqui
BLR(t,z)    = B_bd(t,z-1)./M_d(t,z-1);                                      %11.93


if BLR(t,z) < bot(t,1) 
    z1(t,1)      = 1;
else 
    z1(t,1)      = 0;
end
if BLR(t,z) > top(t,1)  
       z2(t,1)   = 1;
else 
       z2(t,1)   = 0;
end

rm(t,z)     = -zeta_m.*(z1(t,1) - z2(t,1)) + rm(t-1,end);                   %11.95 sinal dessa tá trocado no livro

% 4.2) Retornos dos Empréstimos 

rl(t,z)     = rm(t,z-1) + addl(t,z-1);                                      %11.98
OF_bt(t,z)  = NCAR.*(L_fs(t-1,end) + L_hs(t-1,end));                        %11.99 se o banco for manter essa regra e só dar os empréstimos que todo mundo aceita é muito fácil o spread se tornar negativo, pois, se poucos empréstimos forem tomados, o lucro desejado pelas firmas pode ser muito baixo
OF_be(t,z)  = OF_b(t-1,end) + beta_b.*(OF_bt(t,z-1) - OF_b(t-1,end));       %11.100 o livro diz que um NCAR menor faz as taxas de juros serem maiores, não é isso que acontece
FU_bt(t,z)  = OF_be(t,z-1) - OF_b(t-1,end) + npl_e(t,z-1).*L_fs(t-1,end);   %11.101
npl_e(t,z)  = epsilon_b.*npl_e(t-1,end) + (1 - epsilon_b).*npl(t-1,1);      %11.102
FD_b(t,z)   = lambda_b.*Y(t-1,end);                                         %11.103 porque 
F_bt(t,z)   = FD_b(t,z-1) + FU_bt(t,z-1);                                   %11.104
F_b(t,z)    = rl(t-1,end).*(L_fs(t-1,end) + L_hs(t-1,end) - NPL(t,z-1))...
    + rb(t-1,1).*B_bd(t-1,end) - rm(t-1,end).*M_h(t-1,end);                 %11.105
addl(t,z)   = (F_bt(t,z-1) - rb(t-1,1).*B_bd(t-1,end)...
    + rm(t-1,end).*(M_h(t-1,end) - (1 - npl_e(t,z-1)).*L_fs(t-1,end)...
    - L_hs(t-1,end)))./((1 - npl_e(t,z-1)).*L_fs(t-1,end) + L_hs(t-1,end)); %11.106
FU_b(t,z)   = F_b(t,z-1) - FD_b(t,z-1);                                     %11.107
OF_b(t,z)   = OF_b(t-1,end) + FU_b(t,z-1) - NPL(t,z-1);                     %11.108
CAR(t,z)    = OF_b(t,z-1)./(L_fs(t,z-1) + L_hs(t,z-1));                     %11.109
B_bd(t,z)   = B_bs(t,z-1);                                                  %11.110 essa equaçao não é redundante

    end
end
%%

% Fazendo os Gráficos 

%Vetorizando as variáveis endógenas:

addl              = addl(:,end);
B_bd              = B_bd(:,end);
B_bs              = B_bs(:,end);
B_cbd             = B_cbd(:,end);
B_cbs             = B_cbs(:,end);
B_hd              = B_hd(:,end);
B_hs              = B_hs(:,end);
BL_d              = BL_d(:,end);
BL_s              = BL_s(:,end);
BLR               = BLR(:,end);
B_s               = B_s(:,end);
BUR               = BUR(:,end);
C                 = C(:,end);
c                 = c(:,end);
CAR               = CAR(:,end);
CG                = CG(:,end);
ER                = ER(:,end);
e_d               = e_d(:,end);
e_s               = e_s(:,end);
eta               = eta(:,end);
F_b               = F_b(:,end);
F_bt              = F_bt(:,end);
F_cb              = F_cb(:,end);
F_f               = F_f(:,end);
F_ft              = F_ft(:,end);
FD_b              = FD_b(:,end);
FD_f              = FD_f(:,end);
FU_b              = FU_b(:,end);
FU_bt             = FU_bt(:,end);
FU_f              = FU_f(:,end);
FU_ft             = FU_ft(:,end);
G                 = G(:,end);
GD                = GD(:,end);
GL                = GL(:,end);
g                 = g(:,end);
gr_k              = gr_k(:,end);
H_s               = H_s(:,end);
H_bd              = H_bd(:,end);
H_bs              = H_bs(:,end);
H_hd              = H_hd(:,end);
H_hs              = H_hs(:,end);
HC_e              = HC_e(:,end);
I                 = I(:,end);
IN                = IN(:,end);
i                 = i(:,end); 
in                = in(:,end);
in_e              = in_e(:,end);
in_T              = in_T(:,end);
K                 = K(:,end);
k                 = k(:,end);
L_fd              = L_fd(:,end);
L_fs              = L_fs(:,end);
L_hd              = L_hd(:,end);
L_hs              = L_hs(:,end);
M_d               = M_d(:,end);
M_h               = M_h(:,end);
M_s               = M_s(:,end);
N                 = N(:,end);
N_T               = N_T(:,end);
NL                = NL(:,end);
NPL               = NPL(:,end);
NUC               = NUC(:,end);
NHUC              = NHUC(:,end);
nl                = nl(:,end);
npl_e             = npl_e(:,end); 
omega_T           = omega_T(:,end);
OF_b              = OF_b(:,end);
OF_be             = OF_be(:,end);
OF_bt             = OF_bt(:,end);
PE                = PE(:,end);
PSBR              = PSBR(:,end);
p                 = p(:,end);
pbl               = pbl(:,end);
pe                = pe(:,end);
phi               = phi(:,end);
phi_T             = phi_T(:,end);
pi                = pi(:,end);
pr                = pr(:,end);
q                 = q(:,end);
REP               = REP(:,end);
rbl               = rbl(:,end);
rl                = rl(:,end);
rm                = rm(:,end);
rrl               = rrl(:,end);
rk                = rk(:,end);
S                 = S(:,end);
s                 = s(:,end);
s_e               = s_e(:,end);
sigma_se          = sigma_se(:,end);
T                 = T(:,end);
u                 = u(:,end);
UC                = UC(:,end);
V                 = V(:,end);
V_fma             = V_fma(:,end);
v                 = v(:,end);
W                 = W(:,end);
WB                = WB(:,end);
Y                 = Y(:,end);
YD_r              = YD_r(:,end);
YD_hs             = YD_hs(:,end);
YP                = YP(:,end);
yd_r              = yd_r(:,end);
yd_er             = yd_er(:,end);
y                 = y(:,end);
TxN               = N(t)./N(t-1) - 1;
%%



figure
plot(T-200:T,WB(T-200:T)./Y(T-200:T),'-r','LineWidth',1.5)
grid on
legend('WB/Y')
title('Participação dos Salários')
xlabel('Tempo')
ylabel('Valor')


figure
plot(T-200:T,u(T-200:T),'-m','LineWidth',1.5)
grid on
legend('u')
title('Utilização de Capacidade')
xlabel('Tempo')
ylabel('Valor')


figure
plot(T-200:T,BUR(T-200:T),'-c','LineWidth',1.5)
grid on
legend('BUR')
title('Endividamento das Famílias')
xlabel('Tempo')
ylabel('Valor')


figure
plot(T-200:T,F_b(T-200:T)./F_f(T-200:T),'-c','LineWidth',1.5)
grid on
legend('Fb/Ff')
title('Razão Entre Lucros Bancários e Das Firmas')
xlabel('Tempo')
ylabel('Valor')


figure
plot(700:T,L_fd(700:T)./(K(700:T)+IN(700:T)),'-c','LineWidth',1.5)
grid on
legend('L/(K+IN)')
title('Razão Entre Passivos e Ativos Das Firmas')
xlabel('Tempo')
ylabel('Valor')



