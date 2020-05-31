%Filename   : bss_parameters.m
%Description: Parameters for BSS operation
%
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-06  1.0   Creation
%william         2017-02-21  1.4   Changed R_gl to negative
%william         2017-03-29  1.5   Added price_factor, hour_slice
%======================================================================

%PARAMETERS
%Independent Parameters 
T               = test_tot_time_in_hr * hour_slice;
I               = test_battery_number; 

annual_discount_rate = 0.05; 

R_s     = 20;        %Revenue from Swap
R_l     = 0.50;      %Revenue from Lease
R_gh    = 3;         %High Level Electricity Demand
R_gm    = 2;         %Mid Level Electricity Demand
R_gl    = 1;         %Low Level Electricity Demand

soc_max         = 0.85; 
soc_rdy_full    = 0.80; % (soc_rdy <= soc_max)
soc_rdy_partial = 0.65;
soc_min         = 0.25; 
soc_crt         = 0.40; %arbitrary 

soh_max     = 1.00; 
soh_best    = 0.80; 
soh_betr    = 0.70; 
soh_good    = 0.65;
soh_min     = 0.60; 

eta_c   = 0.95; 
eta_d   = 0.95; 
E_l     = 175; 
E_d     = 50; 

beta    = (0.20/2500);   %2500 cycles down to 20%

%Dependent Parameters 
price_factor = 0.0067/hour_slice;  %Assumes a 30 kWH charged in 4.5 hours
alpha       = (0.20/157680)/hour_slice; %18 years: 157680 hours down to 20%

h_cg    = 4*hour_slice; 
h_dg    = 4*hour_slice; 
h_dh    = 48*hour_slice; 
r_cg    = (soc_max - soc_min)/h_cg; %rate of charging to grid
r_dg    = (soc_max - soc_min)/h_dg; %rate of discharging to grid
r_de    = (soc_max - soc_min)/h_dh; %rate of discharging to ev


load        = test_load_type;
h_customer  = test_customer_data; 
l_customer  = []; 
T_YEARS     = ((T/HOURS_IN_A_DAY)/DAYS_IN_A_YEAR)/hour_slice; 
HOURS_IN_A_YEAR = HOURS_IN_A_DAY*DAYS_IN_A_YEAR; 

%Variables
%Continuous 
soc             = soc_max*ones(I,T);           %State of Charge
soh             = soh_max*ones(I,T);           %State of Health 
u               = zeros(I,T);                  %Used battery cycles 
R_g             = R_gm*ones(1, T);             %Revenue from Grid 
                                               % (default value:R_gm)   

%Integer 
v              = zeros(I,T);                  %Loyalty of customer 
dh             = zeros(1,T);

loyalty        = zeros (I,T);

%VARIABLES (with initial values)    
%flags / binary        
f       = ones(I,T); %(soh >= soh_best)
g       = ones(I,T); %(soh >= soh_betr)
k       = ones(I,T); %(soh >= soh_good)
p       = ones(I,T); %(soc <= soc_crt) 
q       = ones(I,T); %(soc >= soc_rdy_full)

%VALIDATION
transitions = 19;                
modes       = 5; 
states      = 3; 
state_space = zeros(T,30,I); 
current_mode            = zeros(I,T); 
current_transition      = zeros(I,T); 
objective               = zeros(I,T); 
dh_validation           = zeros (I,T);
sh_validation           = zeros (I,T);
m4bss_validation        = zeros (I,T);
loop_validation         = 3*ones(35, transitions);
validate_c25            = zeros (I,T);
validate_c26            = zeros (I,T);