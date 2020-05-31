%Filename   : bss_strategy_time.m
%Description: BSS operation in time per battery
%
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-06  1.0   Creation
%william         2016-12-29  1.1   Added Objectives for Scenarios
%william         2017-01-13  1.3   Added NSWD 2016 Electricity Data
%william         2017-03-29  1.5   Added SA 2016 Electricity Data
%                                  Added Net Present Value
%                                  Changed to 30-minute intervals
%======================================================================

% External Condition A: Generate Electricity Demand
if (hour_slice == 1) 
    if (strcmp(load, 'nswd')) 
        level_hour_year = importdata('ns_lev_365_hr.mat');
        price_hour_year = ...
            price_factor*importdata('ns_rrp_365_hr.mat');
    end
    if (strcmp(load, 'sa')) 
        level_hour_year = importdata('sa_lev_365_hr.mat');
        price_hour_year = ...
            price_factor*importdata('sa_rrp_365_hr.mat');
    end
else %hour_slice == 2 
    if (strcmp(load, 'nswd')) 
        level_hour_year = importdata('ns_lev_365_30.mat');
        price_hour_year = price_factor*importdata('ns_rrp_365_30.mat');
    end
    if (strcmp(load, 'sa')) 
        level_hour_year = importdata('sa_lev_365_30.mat');
        price_hour_year = price_factor*importdata('sa_rrp_365_30.mat');
    end    
end
       
electricity_level = zeros(T, 1);
price_hour_years = zeros(T, 1);  

for year = 1:T_YEARS
    electricity_level(hour_slice*HOURS_IN_A_YEAR*(year - 1) ...
        + 1: hour_slice*HOURS_IN_A_YEAR*(year)) ...
        = level_hour_year(1:hour_slice*HOURS_IN_A_YEAR);

    price_hour_years(hour_slice*HOURS_IN_A_YEAR*(year - 1) ...
        + 1: hour_slice*HOURS_IN_A_YEAR*(year)) ...
        = price_hour_year(1:hour_slice*HOURS_IN_A_YEAR);        
end

t_electricity_level = transpose(electricity_level); 
for t = (1 + 1):T 
    switch t_electricity_level(t)
        case 1
            R_g(t) = R_gl;
        case 2
            R_g(t) = R_gm; 
        case 3
            R_g(t) = R_gh; 
    end        
end

%External Condition B: Generate Customer Data
hct_mx  = zeros(length(h_customer), T); %initial battery demand 
hcb_mx  = zeros(length(h_customer), T); %customers with battery IDs 
hct_mx2 = zeros(length(h_customer), T); %battery demand customer matrix 
for customer_index = 1:length(h_customer)
    hct_mx(customer_index, h_customer(customer_index)) = 1; 
end
hct_mx2 = hct_mx; 
hct_demand = sum(hct_mx); 


%Time sequence operation
for t = (1 + 1):T  %initial value at t = 1 already given    
    dh(t)   = dh(t-1) + hct_demand(t);    
    hct_mx2(:, t) = hct_mx2(:, t-1) + hct_mx(:, t);
    hcb_mx(:, t) = hcb_mx(:, t-1);
    
    for i = 1:I            
        all_conditions = false;                        
        loop = 0;         
        %%generate permutation list ; 
        transition_priority = randperm(transitions); 
        while (all_conditions == false) 
            %%search until satisfied to be used for while condition             
            loop = loop + 1;            
            m0_to_m0(i,t) = false; 
            m1_to_m0(i,t) = false; 
            m2_to_m0(i,t) = false;            
            m4_to_m0(i,t) = false;
            m0_to_m1(i,t) = false;
            m1_to_m1(i,t) = false;
            m2_to_m1(i,t) = false;            
            m4_to_m1(i,t) = false;
            m0_to_m2(i,t) = false;
            m1_to_m2(i,t) = false;
            m2_to_m2(i,t) = false;                       
            m0_to_m4(i,t) = false;
            m1_to_m4(i,t) = false;
            m2_to_m4(i,t) = false;
            m4_to_m4(i,t) = false;
            m0_to_m5(i,t) = false; 
            m1_to_m5(i,t) = false; 
            m2_to_m5(i,t) = false; 
            m5_to_m5(i,t) = false; 

            current_transition(i,t) = transition_priority(loop);            
            switch current_transition(i,t) 
                case 1 
                    m0_to_m0(i,t) = true; 
                case 2 
                    m1_to_m0(i,t)  = true; 
                case 3 
                    m2_to_m0(i,t)  = true;   
                case 4 
                    m4_to_m0(i,t)  = true;  
                case 5 
                    m0_to_m1(i,t)  = true;  
                case 6 
                    m1_to_m1(i,t)  = true; 
                case 7 
                    m2_to_m1(i,t)  = true; 
                case 8 
                    m4_to_m1(i,t)  = true;                     
                case 9 
                    m0_to_m2(i,t)  = true; 
                case 10 
                    m1_to_m2(i,t)  = true; 
                case 11 
                    m2_to_m2(i,t)  = true; 
                case 12
                    m0_to_m4(i,t)  = true; 
                case 13
                    m1_to_m4(i,t)  = true; 
                case 14
                    m2_to_m4(i,t)  = true; 
                case 15
                    m4_to_m4(i,t)  = true; 
                case 16
                    m0_to_m5(i,t)  = true; 
                case 17
                    m1_to_m5(i,t)  = true;                       
                case 18
                    m2_to_m5(i,t)  = true;                       
                case 19
                    m5_to_m5(i,t)  = true;                 
                otherwise
                    warning('Unexpected value');
            end
            
            %DEFINITIONS 
            %Modes
            m0(i,t) = m0_to_m0(i,t) + m1_to_m0(i,t) + m2_to_m0(i,t) + ...
                      m4_to_m0(i,t); 
            m1(i,t) = m0_to_m1(i,t) + m1_to_m1(i,t) + m2_to_m1(i,t) + ...
                      m4_to_m1(i,t); 
            m2(i,t) = m0_to_m2(i,t) + m1_to_m2(i,t) + m2_to_m2(i,t);  
            m4(i,t) = m0_to_m4(i,t) + m1_to_m4(i,t) + m2_to_m4(i,t) + ...
                      m4_to_m4(i,t);
            m5(i,t) = m0_to_m5(i,t) + m1_to_m5(i,t) + m2_to_m5(i,t) + ...
                      m5_to_m5(i,t);

            if(m0(i,t) == 1)
                current_mode(i,t) = 0;
            end
            if(m1(i,t) == 1)
                current_mode(i,t) = 1;
            end              
            if(m2(i,t) == 1)
                current_mode(i,t) = 2;
            end  
            if(m4(i,t) == 1)
                current_mode(i,t) = 4;
            end                  
            if(m5(i,t) == 1)
                current_mode(i,t) = 5;
            end
            
            %States
            b0(i,t) = m0_to_m0(i,t) + m1_to_m0(i,t) + m2_to_m0(i,t) + ...
                      m4_to_m0(i,t) + m0_to_m5(i,t) + ...
                      m1_to_m5(i,t) + m2_to_m5(i,t) + m5_to_m5(i,t);
            bc(i,t) = m0_to_m1(i,t) + m1_to_m1(i,t) + m2_to_m1(i,t) + ...
                      m4_to_m1(i,t);
            bd(i,t) = m0_to_m2(i,t) + m1_to_m2(i,t) + m2_to_m2(i,t) + ...
                      m0_to_m4(i,t) + m1_to_m4(i,t) + ...
                      m2_to_m4(i,t) + m4_to_m4(i,t);

            %Functions
            m_w(i,t)     = m0_to_m0(i,t) + ...
                           m1_to_m0(i,t) + ...
                           m2_to_m0(i,t) + ...
                           m4_to_m0(i,t);
            m_c(i,t)     = m0_to_m1(i,t) + ...
                           m1_to_m1(i,t) + ...
                           m2_to_m1(i,t) + ...
                           m4_to_m1(i,t);
            m_d(i,t)     = m0_to_m2(i,t) + ...
                           m1_to_m2(i,t) + ...
                           m2_to_m2(i,t); 
            m_s(i,t)     = m0_to_m4(i,t) + ...
                           m1_to_m4(i,t) + ...
                           m2_to_m4(i,t); 
            m_e(i,t)     = m4_to_m4(i,t);  
            m_r(i,t)     = m0_to_m5(i,t) + ...
                           m1_to_m5(i,t) + ...
                           m2_to_m5(i,t); 
            m_f(i,t)     = m5_to_m5(i,t);

            %Extra Definitions
            m_b(i,t)     = m0(i,t) + m1(i,t) + m2(i,t);
            m_return(i,t)    = m4_to_m0(i,t) + m4_to_m1(i,t);
            
            
            c01_special_order_v1 = (1 == ...
                m0_to_m0(i,t) + m1_to_m0(i,t) + m2_to_m0(i,t) + ...
                m4_to_m0(i,t) + ...
                m0_to_m1(i,t) + m1_to_m1(i,t) + m2_to_m1(i,t) + ...
                m4_to_m1(i,t) + ...
                m0_to_m2(i,t) + m1_to_m2(i,t) + m2_to_m2(i,t) + ...
                m0_to_m4(i,t) + m1_to_m4(i,t) + m2_to_m4(i,t) + ...
                m4_to_m4(i,t) + ...
                m0_to_m5(i,t) + m1_to_m5(i,t) + m2_to_m5(i,t) + ...
                m5_to_m5(i,t));            
             
            %c02_soc_computation
            soc(i,t) = soc(i,t-1) + eta_c*r_cg*m1(i,t) ...
                - (1/eta_d)*(r_dg*m2(i,t) + r_de*m4(i,t));
                        
            c03_soc_upper_limit = (soc(i,t) <= soc_max + r_cg);             
            c04_soc_lower_limit = (soc_min <= soc(i,t)); 
            
            %c05_soh_use 
            u(i,t) = u(i,t-1) + (1/2)*(r_cg*m1(i,t) + r_dg*m2(i,t) + ...
                     r_de*m4(i,t)); 
            
            %c06_soh_computation
            soh(i,t) = (1 - alpha*t - beta*u(i,t))*k(i,t-1) + ...
                soh(i,t-1)*(1 - k(i,t-1));   
                        
            
            c07_soh_upper_limit = (soh(i,t) <= soh_max);             
            c08_soh_lower_limit = (soh_min <= soh(i,t));             
            c09_soh_sequence = (soh(i,t) <= soh(i,t-1)); 
                        
            f(i,t) = (soh(i,t) >= soh_best); 
            g(i,t) = (soh(i,t) >= soh_betr);                        
            k(i,t) = (soh(i,t) >= soh_good);
            p(i,t) = (soc(i,t) <= soc_crt); 
            q(i,t) = (soc(i,t-1) >= soc_rdy_full); 
            
            c10_soh_high_ev_use = (m_s(i,t) <= f(i,t));             
            c12_soh_grid_only = (m0(i,t) + m1(i,t) + m2(i,t) <= k(i,t-1));              
            c13_soh_disposal = (m5(i,t) <= (1 - k(i,t-1)) ); 
            
            c14_strategy_m_0 = (m0(i,t-1) <= m0_to_m0(i,t) + ...
                m0_to_m1(i,t) + m0_to_m2(i,t)  + ...
                m0_to_m4(i,t) + m0_to_m5(i,t));
            c15_strategy_m_1 = (m1(i,t-1) <= m1_to_m0(i,t) + ...
                m1_to_m1(i,t) + m1_to_m2(i,t)  + ...
                m1_to_m4(i,t) + m1_to_m5(i,t));            
            c16_strategy_m_2 = (m2(i,t-1) <= m2_to_m0(i,t) + ...
                m2_to_m1(i,t) + m2_to_m2(i,t)  + ...
                m2_to_m4(i,t) + m2_to_m5(i,t));                   
            c18_strategy_m_4 = (m4(i,t-1) <= m4_to_m0(i,t) + ...
                m4_to_m1(i,t) + m4_to_m4(i,t));    
            c19_strategy_m_5 = (m5(i,t-1) <= m5_to_m5(i,t));    
            
            c21_ev_high      = (m4(i,t-1) <= m_e(i,t) + p(i,t));  %soc_crt                                                
            c22_reliability_hev = ( (m_s(i,t)) == ...
                ((dh(t) > 0) & q(i,t) & m_b(i,t-1) & f(i,t)));                                                                  
            c24_discharge_rel = ( ((dh(t) > 0) & f(i,t)) <= ...
                (1 - m2(i,t) - m2_to_m0(i,t) - m1_to_m0(i,t)) );
            
           %Additional constraints for smart grid services 
            c25_elec_high = ( ( (dh(t) == 0) & q(i,t) & m_b(i,t-1) & ...
                k(i,t) & (R_g(t) == R_gh)) <=  m2(i,t) );              
            c26_elec_low = ( ( (dh(t) == 0) & (q(i,t) == 0) & ...
                m_b(i,t-1) & k(i,t) & (R_g(t) == R_gl) )  <=  m1(i,t) );
            c27_elec_mid = ( ( (dh(t) == 0) & q(i,t) & ...
                m_b(i,t-1) & k(i,t) & (R_g(t) == R_gm) )  <=  m0(i,t) );           
            validate_c26(i,t) =  ((dh(t) == 0) & (q(i,t) == 0) & ...
                m_b(i,t-1) & k(i,t) & (R_g(t) == R_gl)); 
            
            if (has_smart_grid_services == false)                
                    all_conditions = c01_special_order_v1   && ...
                                 c03_soc_upper_limit        && ...
                                 c04_soc_lower_limit        && ...
                                 c07_soh_upper_limit        && ...
                                 c08_soh_lower_limit        && ...
                                 c09_soh_sequence           && ...
                                 c10_soh_high_ev_use        && ...
                                 c12_soh_grid_only          && ...
                                 c13_soh_disposal           && ...
                                 c14_strategy_m_0           && ...
                                 c15_strategy_m_1           && ...
                                 c16_strategy_m_2           && ...
                                 c18_strategy_m_4           && ...
                                 c19_strategy_m_5           && ...
                                 c21_ev_high                && ...
                                 c22_reliability_hev        && ...
                                 c24_discharge_rel          ;                                     
            else 
                if (reserved_batt == 0) 
                    all_conditions = c01_special_order_v1   && ...
                                 c03_soc_upper_limit        && ...
                                 c04_soc_lower_limit        && ...
                                 c07_soh_upper_limit        && ...
                                 c08_soh_lower_limit        && ...
                                 c09_soh_sequence           && ...
                                 c10_soh_high_ev_use        && ...
                                 c12_soh_grid_only          && ...
                                 c13_soh_disposal           && ...
                                 c14_strategy_m_0           && ...
                                 c15_strategy_m_1           && ...
                                 c16_strategy_m_2           && ...
                                 c18_strategy_m_4           && ...
                                 c19_strategy_m_5           && ...
                                 c21_ev_high                && ...
                                 c22_reliability_hev        && ...
                                 c24_discharge_rel          && ... 
                                 c25_elec_high              && ...                
                                 c26_elec_low               && ...
                                 c27_elec_mid    ;       
                else 
                    all_conditions = c01_special_order_v1   && ...
                                 c03_soc_upper_limit        && ...
                                 c04_soc_lower_limit        && ...
                                 c07_soh_upper_limit        && ...
                                 c08_soh_lower_limit        && ...
                                 c09_soh_sequence           && ...
                                 c10_soh_high_ev_use        && ...
                                 c12_soh_grid_only          && ...
                                 c13_soh_disposal           && ...
                                 c14_strategy_m_0           && ...
                                 c15_strategy_m_1           && ...
                                 c16_strategy_m_2           && ...
                                 c18_strategy_m_4           && ...
                                 c19_strategy_m_5           && ...
                                 c21_ev_high                && ...
                                 c22_reliability_hev        && ...
                                 c24_discharge_rel          && ...                                              
                                 c26_elec_low               && ...
                                 c27_elec_mid               ;    
                end
            end
             
            if (all_conditions == true) 
                dh_validation(i,t) = dh(t);

                if (m_s(i,t))
                    dh(t) = dh(t) - 1;     
                    sh_validation(i,t) = 1;                     
                    first_index = find(hct_mx2(:, t)==1,1); 
                    hcb_mx(first_index, t) = i;
                    hct_mx2(first_index, t) = 0;                                           
                end
                if (m4_to_m0(i,t))
                    dh(t) = dh(t) + 1;    
                    first_index = find(hcb_mx(:, t)==i,1); 
                    hcb_mx(first_index, t) = 0;
                    hct_mx2(first_index, t) = 1;                                           
                    m4bss_validation(i,t) = 1; 
                end         
                if (m4_to_m1(i,t))
                    dh(t) = dh(t) + 1;                         
                    first_index = find(hcb_mx(:, t)==i,1); 
                    hcb_mx(first_index, t) = 0;
                    hct_mx2(first_index, t) = 1;                      
                    m4bss_validation(i,t) = 1; 
                end
            end
        end 
    end     
end 

real_m_d = bsxfun(@times, transpose(price_hour_years), m_d); 
real_m_c = bsxfun(@times, transpose(price_hour_years), m_c);

increment_objective =   R_s*(m_s) ...
                      + R_l*(m_s + m_e) - E_l*(v) ...
                      - E_d*(m_r) ...
                      + real_m_d ...
                      - real_m_c;
obj_year = zeros(T_YEARS,1); 
disc_obj_year = zeros(T_YEARS,1); 

for year = 1:T_YEARS
    obj_year(year) = sum(sum( ...
        increment_objective(:,hour_slice*HOURS_IN_A_YEAR*(year - 1) + 1 ...
        : hour_slice*HOURS_IN_A_YEAR*(year))) );           
    disc_obj_year(year) = obj_year(year)/(1+annual_discount_rate)^year;
end                 
                  
obj_s1_l1_g1 = sum(sum(R_s*(m_s) ...
                       + R_l*(m_s + m_e) - E_l*(v) ...
                       - E_d*(m_r) ...
                       + real_m_d ...
                       - real_m_c ...
                   )...
                );
                        
m_d_tot =  sum(sum(real_m_d));
m_c_tot =  sum(sum(real_m_c)); 
m_l_tot = sum(sum(R_l*(m_s + m_e)));
m_v_tot = sum(sum(E_l*(v))); 
m_s_tot = sum(sum(m_s)); 
    
defective_batt = zeros(I,1); 
for batteries = 1:I
    defective_batt(batteries) = find(m5(batteries,:) == 1, 1);     
end

soh_good_batt = zeros(I,1); 
for batteries = 1:I
    soh_good_batt(batteries) = find(soh(batteries,:) <= soh_good, 1);     
end

soh_betr_batt = zeros(I,1); 
for batteries = 1:I
    soh_betr_batt(batteries) = find(soh(batteries,:) <= soh_betr, 1);     
end

soh_best_batt = zeros(I,1); 
for batteries = 1:I
    soh_best_batt(batteries) = find(soh(batteries,:) <= soh_best, 1);     
end

yy = min(soh_best_batt); 
hct_mx2c = hct_mx2(:,1:yy); 
add_column = zeros(length(h_customer),1);

hct_mx2ch = (abs(diff([hct_mx2c add_column], 1, 2)));
hct_mx2ch3 = 1 - (abs(diff([hct_mx2ch add_column], 1, 2)));

hct_x = hct_mx2ch3 .* hct_mx2c; 
hct_x_sum = sum(sum(hct_x));

%Transfer to Excel
if(isempty(h_customer)) 
    scenario = 'pur';
else 
    scenario = 'mix';
end

FileName = ...
sprintf('%s_bss_%s_grid_%d_res_%d_%s_test_seed_%d.xlsx',  ...
datestr(now, 'yymmdd_HHMM'), scenario, ...
has_smart_grid_services, ...
reserved_batt, load, test_random_number);    

f_mx.filename               = FileName;      
f_mx.objective_type         = objective_type;
f_mx.objective              = obj_s1_l1_g1;
f_mx.objective_revised      = obj_s1_l1_g1 - E_l*hct_x_sum;
f_mx.npv_obj                = sum(disc_obj_year); 
f_mx.seed                   = test_random_number;
f_mx.has_grid_scheduling    = has_smart_grid_services;
f_mx.reservation            = reserved_batt;
f_mx.ev_customers           = length(h_customer); 
f_mx.battery_count          = I;
f_mx.m_s_tot                = m_s_tot;
f_mx.m_d_tot                = m_d_tot;
f_mx.m_l_tot                = m_l_tot;
f_mx.m_c_tot                = m_c_tot;
f_mx.m_v_tot                = m_v_tot;

for counting_c = 1:I
 batt_charging_count(counting_c,:)   = R_g .* m_c(counting_c,:);
end

for counting_d = 1:I
 batt_discharging_count(counting_d,:)   = R_g .* m_c(counting_d,:);
end

f_mx.sum_m_s                        = sum(sum(m_s));
f_mx.sum_m_e                        = sum(sum(m_e));
f_mx.sum_v                          = hct_x_sum;
f_mx.sum_m_r                        = sum(sum(m_r));
f_mx.sum_m_c                        = nnz(m_c);
f_mx.sum_m_c_high                   = nnz(batt_charging_count == R_gh);
f_mx.sum_m_c_mid                    = nnz(batt_charging_count == R_gm);
f_mx.sum_m_c_low                    = nnz(batt_charging_count == R_gl);
f_mx.sum_m_d                        = nnz(m_d);
f_mx.sum_m_d_high                   = nnz(batt_discharging_count == R_gl);
f_mx.sum_m_d_mid                    = nnz(batt_discharging_count == R_gm);
f_mx.sum_m_d_low                    = nnz(batt_discharging_count == R_gh);   
f_mx.load                           = load; 

table_f_mx            = struct2table(f_mx); 
cell_f_mx             = table2cell(table_f_mx); 
cell_f_mx_with_header = [table_f_mx.Properties.VariableNames;  cell_f_mx]; 

sheetnames = ...
    {'final_state', 'defective','soh_good', 'soh_betr', 'soh_best'}; 
xlsheets(sheetnames, FileName);   
xlswrite(FileName, defective_batt, 'defective');
xlswrite(FileName, soh_good_batt, 'soh_good');
xlswrite(FileName, soh_betr_batt, 'soh_betr');
xlswrite(FileName, soh_best_batt, 'soh_best');
xlswrite(FileName, transpose(cell_f_mx_with_header), 'final_state');