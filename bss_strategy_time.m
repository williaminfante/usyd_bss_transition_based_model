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
%======================================================================

%External Condition A: Generate Electricity Demand
if (strcmp(load, 'duck')) 
    heavy_interval = [19 20 21 22];
    light_interval = [12 13 14 15];
elseif (strcmp(load, 'nswd')) 
    heavy_interval = [17 18 19 20];
    light_interval = [2 3 4 5];    
elseif (strcmp(load, 'mojo'))
    heavy_interval = [15 16 17 18];
    light_interval = [0 1 2 3];
else
    error('no type of load existing'); 
end

if (stage == 0)
    for t = (1 + 1):T 
        if (ismember(mod(t, 24), heavy_interval))
            R_g(t) = R_gh;
        elseif (ismember(mod(t, 24), light_interval))
            R_g(t) = R_gl;
        else
            R_g(t) = R_gm; 
        end        
    end
else
    %Stage 1 & 2 : Based on 2016 NSWD Level Pricing
    level_hour_365 = importdata('level_hour_365.mat');
    level_hour_15_years = zeros(T, 1);
    for year = 1:15
        level_hour_15_years(8760*(year - 1) + 1: 8760*(year)) ...
            = level_hour_365(1:8760);
    end
    t_level_15_years = transpose(level_hour_15_years); 
    for t = (1 + 1):T 
        switch t_level_15_years(t)
            case 1
                R_g(t) = R_gl;
            case 2
                R_g(t) = R_gm; 
            case 3
                R_g(t) = R_gh; 
        end        
    end
end


%Stage 2: Based on PROPORTIONAL pricing
    
if (stage == 2) 
    price_hour_365 = importdata('2016_elec_price_365_days.mat');
    price_hour_15_years = zeros(T, 1);    
    for year = 1:15
        price_hour_15_years(8760*(year - 1) + 1: 8760*(year)) ...
            = price_hour_365(1:8760);
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

lct_mx = zeros(length(l_customer), T);
for customer_index = 1:length(l_customer)
    lct_mx(customer_index, l_customer(customer_index)) = 1; 
end
lct_demand = sum(lct_mx); 


%Time sequence operation
for t = (1 + 1):T  %initial value at t = 1 already given    
    dh(t)   = dh(t-1) + hct_demand(t);
    dl(t)   = dl(t-1) + lct_demand(t);
    
    hct_mx2(:, t) = hct_mx2(:, t-1) + hct_mx(:, t);
    hcb_mx(:, t) = hcb_mx(:, t-1);
    for i = 1:I            
        all_conditions = false;                        
        loop = 0;         
        %%generate permutation list ; 
        transition_priority = randperm(25); 
        while (all_conditions == false) 
            %%search until satisfied to be used for while condition             
            loop = loop + 1; 
            
            m0_to_m0(i,t) = false; 
            m1_to_m0(i,t) = false; 
            m2_to_m0(i,t) = false;
            m3_to_m0(i,t) = false; 
            m4_to_m0(i,t) = false;
            m0_to_m1(i,t) = false;
            m1_to_m1(i,t) = false;
            m2_to_m1(i,t) = false;
            m3_to_m1(i,t) = false;
            m4_to_m1(i,t) = false;
            m0_to_m2(i,t) = false;
            m1_to_m2(i,t) = false;
            m2_to_m2(i,t) = false;
            m0_to_m3(i,t) = false;
            m1_to_m3(i,t) = false;
            m2_to_m3(i,t) = false;
            m3_to_m3(i,t) = false;
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
                    current_mode(i,t) = 0;
                case 2 
                    m1_to_m0(i,t)  = true; 
                    current_mode(i,t) = 0;
                case 3 
                    m2_to_m0(i,t)  = true;   
                    current_mode(i,t) = 0;
                case 4 
                    m3_to_m0(i,t)  = true;  
                    current_mode(i,t) = 0';
                case 5 
                    m4_to_m0(i,t)  = true;  
                    current_mode(i,t) = 0;
                case 6 
                    m0_to_m1(i,t)  = true;  
                    current_mode(i,t) = 1;
                case 7 
                    m1_to_m1(i,t)  = true; 
                    current_mode(i,t) = 1;
                case 8 
                    m2_to_m1(i,t)  = true; 
                    current_mode(i,t) = 1;
                case 9 
                    m3_to_m1(i,t)  = true; 
                    current_mode(i,t) = 1;
                case 10 
                    m4_to_m1(i,t)  = true; 
                    current_mode(i,t) = 1;
                case 11 
                    m0_to_m2(i,t)  = true; 
                    current_mode(i,t) = 2;
                case 12 
                    m1_to_m2(i,t)  = true; 
                    current_mode(i,t) = 2;
                case 13 
                    m2_to_m2(i,t)  = true; 
                    current_mode(i,t) = 2;
                case 14 
                    m0_to_m3(i,t)  = true; 
                    current_mode(i,t) = 3;
                case 15 
                    m1_to_m3(i,t)  = true; 
                    current_mode(i,t) = 3;
                case 16
                    m2_to_m3(i,t)  = true; 
                    current_mode(i,t) = 3;
                case 17
                    m3_to_m3(i,t)  = true; 
                    current_mode(i,t) = 3;
                case 18
                    m0_to_m4(i,t)  = true; 
                    current_mode(i,t) = 4;
                case 19
                    m1_to_m4(i,t)  = true; 
                    current_mode(i,t) = 4;
                case 20
                    m2_to_m4(i,t)  = true; 
                    current_mode(i,t) = 4;
                case 21
                    m4_to_m4(i,t)  = true; 
                    current_mode(i,t) = 4;
                case 22
                    m0_to_m5(i,t)  = true; 
                    current_mode(i,t) = 5;
                case 23
                    m1_to_m5(i,t)  = true; 
                    current_mode(i,t) = 5;                       
                case 24
                    m2_to_m5(i,t)  = true; 
                    current_mode(i,t) = 5;                       
                case 25
                    m5_to_m5(i,t)  = true; 
                    current_mode(i,t) = 5;                    
                otherwise
                    warning('Unexpected value');
            end
            
            %DEFINITIONS 
            %Modes
            m0(i,t) = m0_to_m0(i,t) + m1_to_m0(i,t) + m2_to_m0(i,t) + ...
                      m3_to_m0(i,t) + m4_to_m0(i,t); 
            m1(i,t) = m0_to_m1(i,t) + m1_to_m1(i,t) + m2_to_m1(i,t) + ...
                      m3_to_m1(i,t) + m4_to_m1(i,t); 
            m2(i,t) = m0_to_m2(i,t) + m1_to_m2(i,t) + m2_to_m2(i,t);  
            m3(i,t) = m0_to_m3(i,t) + m1_to_m3(i,t) + m2_to_m3(i,t) + ...
                      m3_to_m3(i,t);
            m4(i,t) = m0_to_m4(i,t) + m1_to_m4(i,t) + m2_to_m4(i,t) + ...
                      m4_to_m4(i,t);
            m5(i,t) = m0_to_m5(i,t) + m1_to_m5(i,t) + m2_to_m5(i,t) + ...
                      m5_to_m5(i,t);

            %States
            b0(i,t) = m0_to_m0(i,t) + m1_to_m0(i,t) + m2_to_m0(i,t) + ...
                      m3_to_m0(i,t) + m4_to_m0(i,t) + m0_to_m5(i,t) + ...
                      m1_to_m5(i,t) + m2_to_m5(i,t) + m5_to_m5(i,t);
            bc(i,t) = m0_to_m1(i,t) + m1_to_m1(i,t) + m2_to_m1(i,t) + ...
                      m3_to_m1(i,t) + m4_to_m1(i,t);
            bd(i,t) = m0_to_m2(i,t) + m1_to_m2(i,t) + m2_to_m2(i,t) + ...
                      m0_to_m3(i,t) + m1_to_m3(i,t) + m2_to_m3(i,t) + ...
                      m3_to_m3(i,t) + m0_to_m4(i,t) + m1_to_m4(i,t) + ...
                      m2_to_m4(i,t) + m4_to_m4(i,t);

            %Functions
            m_w(i,t)     = m0_to_m0(i,t) + ...
                           m1_to_m0(i,t) + ...
                           m2_to_m0(i,t) + ...
                           m3_to_m0(i,t) + ...
                           m4_to_m0(i,t);
            m_c(i,t)     = m0_to_m1(i,t) + ...
                           m1_to_m1(i,t) + ...
                           m2_to_m1(i,t) + ...
                           m3_to_m1(i,t) + ...
                           m4_to_m1(i,t);
            m_d(i,t)     = m0_to_m2(i,t) + ...
                           m1_to_m2(i,t) + ...
                           m2_to_m2(i,t); 
            m_s(i,t)     = m0_to_m3(i,t) + ...
                           m1_to_m3(i,t) + ...
                           m2_to_m3(i,t) + ...
                           m0_to_m4(i,t) + ...
                           m1_to_m4(i,t) + ...
                           m2_to_m4(i,t);
            m_sl(i,t)    = m0_to_m3(i,t) + ...
                           m1_to_m3(i,t) + ...
                           m2_to_m3(i,t); 
            m_sh(i,t)    = m0_to_m4(i,t) + ...
                           m1_to_m4(i,t) + ...
                           m2_to_m4(i,t); 
            m_e(i,t)     = m3_to_m3(i,t) + ...
                           m4_to_m4(i,t); 
            m_el(i,t)    = m3_to_m3(i,t); 
            m_eh(i,t)    = m4_to_m4(i,t); 
            m_r(i,t)     = m0_to_m5(i,t) + ...
                           m1_to_m5(i,t) + ...
                           m2_to_m5(i,t); 
            m_f(i,t)     = m5_to_m5(i,t);

            %Extra Definitions
            m_b(i,t)     = m0(i,t) + m1(i,t) + m2(i,t);
            m_rh(i,t)    = m4_to_m0(i,t) + m4_to_m1(i,t);
            
            
            c01_special_order_v1 = (1 == ...
                m0_to_m0(i,t) + m1_to_m0(i,t) + m2_to_m0(i,t) + ...
                m3_to_m0(i,t) + m4_to_m0(i,t) + ...
                m0_to_m1(i,t) + m1_to_m1(i,t) + m2_to_m1(i,t) + ...
                m3_to_m1(i,t) + m4_to_m1(i,t) + ...
                m0_to_m2(i,t) + m1_to_m2(i,t) + m2_to_m2(i,t) + ...
                m0_to_m3(i,t) + m1_to_m3(i,t) + m2_to_m3(i,t) + ...
                m3_to_m3(i,t) + ...
                m0_to_m4(i,t) + m1_to_m4(i,t) + m2_to_m4(i,t) + ...
                m4_to_m4(i,t) + ...
                m0_to_m5(i,t) + m1_to_m5(i,t) + m2_to_m5(i,t) + ...
                m5_to_m5(i,t));            
%             c01_special_order_v2 = (1 == b0(i,t) + bc(i,t) + bd(i,t)); 
%             c01_special_order_v3 = (1 == m0(i,t) + m1(i,t) + m2(i,t) +...
%                 m3(i,t) + m4(i,t) + m5(i,t)); 
%             c01_special_order_v4 = (1 == m_w(i,t) + m_c(i,t) + ...
%                 m_d(i,t) + m_r(i,t) + m_sl(i,t) + m_sh(i,t) + ...
%                 m_el(i,t) + m_eh(i,t) + m_f(i,t)); 
%             c01_special_order_v5 = (1 == m_w(i,t) + m_c(i,t) + ...
%                 m_d(i,t) + m_r(i,t) + m_s(i,t) + m_e(i,t) + m_f(i,t));  
            
            %c02_soc_computation
            soc(i,t) = soc(i,t-1) + eta_c*r_cg*m1(i,t) ...
                - (1/eta_d)*(r_dg*m2(i,t) + r_dl*m3(i,t) + r_dh*m4(i,t));
                        
            c03_soc_upper_limit = (soc(i,t) <= soc_max + r_cg);             
            c04_soc_lower_limit = (soc_min <= soc(i,t)); 
            
            %c05_soh_use 
            u(i,t) = u(i,t-1) + (1/2)*(r_cg*m1(i,t) + r_dg*m2(i,t) + ...
                     r_dl*m3(i,t) + r_dh*m4(i,t)); 
            
            %c06_soh_computation
            soh(i,t) = (1 - alpha*t - beta*u(i,t))*k(i,t-1) + ...
                soh(i,t-1)*(1 - k(i,t-1));
            
            c07_soh_upper_limit = (soh(i,t) <= soh_max); 
            
            c08_soh_lower_limit = (soh_min <= soh(i,t)); 
            
            c09_soh_sequence = (soh(i,t) <= soh(i,t-1)); 
                        
            f(i,t) = (soh(i,t) >= soh_best); 
            c10_soh_high_ev_use = (m_sh(i,t) <= f(i,t)); 
            
            g(i,t) = (soh(i,t) >= soh_betr);
            c11_soh_low_ev_use = (m_sl(i,t) <= g(i,t));
                        
            k(i,t) = (soh(i,t) >= soh_good);
            c12_soh_grid_only = (m0(i,t) + m1(i,t) + m2(i,t) <= k(i,t-1));              
            c13_soh_disposal = (m5(i,t) <= (1 - k(i,t-1)) ); 

            c14_strategy_m_0 = (m0(i,t-1) <= m0_to_m0(i,t) + ...
                m0_to_m1(i,t) + m0_to_m2(i,t) + m0_to_m3(i,t) + ...
                m0_to_m4(i,t) + m0_to_m5(i,t));
            c15_strategy_m_1 = (m1(i,t-1) <= m1_to_m0(i,t) + ...
                m1_to_m1(i,t) + m1_to_m2(i,t) + m1_to_m3(i,t) + ...
                m1_to_m4(i,t) + m1_to_m5(i,t));            
            c16_strategy_m_2 = (m2(i,t-1) <= m2_to_m0(i,t) + ...
                m2_to_m1(i,t) + m2_to_m2(i,t) + m2_to_m3(i,t) + ...
                m2_to_m4(i,t) + m2_to_m5(i,t));            
            c17_strategy_m_3 = (m3(i,t-1) <= m3_to_m0(i,t) + ...
                m3_to_m1(i,t) + m3_to_m3(i,t));            
            c18_strategy_m_4 = (m4(i,t-1) <= m4_to_m0(i,t) + ...
                m4_to_m1(i,t) + m4_to_m4(i,t));    
            c19_strategy_m_5 = (m5(i,t-1) <= m5_to_m5(i,t));                 
            
            p(i,t) = (soc(i,t) <= soc_crt); 
            c20_ev_low       = (m3(i,t-1) <= m_el(i,t) + p(i,t));  %soc_crt                    
            c21_ev_high      = (m4(i,t-1) <= m_eh(i,t) + p(i,t));  %soc_crt            
                        
            q1(i,t) = (soc(i,t-1) >= soc_rdy_full); 
            c22_reliability_hev = ( (m_sh(i,t)) == ...
                ((dh(t) > 0) & q1(i,t) & m_b(i,t-1) & f(i,t)));
                                  
            q2(i,t) = (soc(i,t-1) >= soc_rdy_partial);
            c23_reliability_lev = ( m_sl(i,t) == ...
                ( (dl(t) > 0) & q2(i,t) & m_b(i,t-1) & g(i,t)));            
                        
            c24_discharge_rel = ( ((dh(t) > 0) & f(i,t)) <= ...
                (1 - m2(i,t) - m2_to_m0(i,t) - m1_to_m0(i,t)) );
            c25_discharge_rel = ( ((dl(t) > 0) & g(i,t)) <= ...
                (1 - m2(i,t) - m2_to_m0(i,t) - m1_to_m0(i,t)) );
%             %Additional constraints for smart grid services 
            c25_elec_high = ( ( (dh(t) == 0) & q2(i,t) & m_b(i,t-1) & ...
                k(i,t) & (R_g(t) == R_gh)) <=  m2(i,t) );              
            c26_elec_low = ( ( (dh(t) == 0) & (q2(i,t) == 0) & ...
                m_b(i,t-1) & k(i,t) & (R_g(t) == R_gl) )  <=  m1(i,t) );
            c27_elec_mid = ( ( (dh(t) == 0) & q1(i,t) & ...
                m_b(i,t-1) & k(i,t) & (R_g(t) == R_gm) )  <=  m0(i,t) );
            
            validate_c25(i,t) = ( (dh(t) == 0) & q2(i,t) & ...
                m_b(i,t-1) & k(i,t) & (R_g(t) == R_gh) ); 
            validate_c26(i,t) =  ((dh(t) == 0) & (q2(i,t) == 0) & ...
                m_b(i,t-1) & k(i,t) & (R_g(t) == R_gl)); 
            
            if (has_smart_grid_services == false)                
                    all_conditions = c01_special_order_v1   && ...
                                 c03_soc_upper_limit        && ...
                                 c04_soc_lower_limit        && ...
                                 c07_soh_upper_limit        && ...
                                 c08_soh_lower_limit        && ...
                                 c09_soh_sequence           && ...
                                 c10_soh_high_ev_use        && ...
                                 c11_soh_low_ev_use         && ...
                                 c12_soh_grid_only          && ...
                                 c13_soh_disposal           && ...
                                 c14_strategy_m_0           && ...
                                 c15_strategy_m_1           && ...
                                 c16_strategy_m_2           && ...
                                 c17_strategy_m_3           && ...
                                 c18_strategy_m_4           && ...
                                 c19_strategy_m_5           && ...
                                 c20_ev_low                 && ...
                                 c21_ev_high                && ...
                                 c22_reliability_hev        && ...
                                 c23_reliability_lev        && ...
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
                                 c11_soh_low_ev_use         && ...
                                 c12_soh_grid_only          && ...
                                 c13_soh_disposal           && ...
                                 c14_strategy_m_0           && ...
                                 c15_strategy_m_1           && ...
                                 c16_strategy_m_2           && ...
                                 c17_strategy_m_3           && ...
                                 c18_strategy_m_4           && ...
                                 c19_strategy_m_5           && ...
                                 c20_ev_low                 && ...
                                 c21_ev_high                && ...
                                 c22_reliability_hev        && ...
                                 c23_reliability_lev        && ...
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
                                 c11_soh_low_ev_use         && ...
                                 c12_soh_grid_only          && ...
                                 c13_soh_disposal           && ...
                                 c14_strategy_m_0           && ...
                                 c15_strategy_m_1           && ...
                                 c16_strategy_m_2           && ...
                                 c17_strategy_m_3           && ...
                                 c18_strategy_m_4           && ...
                                 c19_strategy_m_5           && ...
                                 c20_ev_low                 && ...
                                 c21_ev_high                && ...
                                 c22_reliability_hev        && ...
                                 c23_reliability_lev        && ...
                                 c24_discharge_rel          && ...                                              
                                 c26_elec_low               && ...
                                 c27_elec_mid               ;    
                end
            end


             loop_validation(1, loop)  = c01_special_order_v1;
             loop_validation(2, loop)  = soc(i,t);
             loop_validation(3, loop)  = c03_soc_upper_limit;
             loop_validation(4, loop)  = c04_soc_lower_limit;
             loop_validation(5, loop)  = u(i,t);
             loop_validation(6, loop)  = soh(i,t);
             loop_validation(7, loop)  = c07_soh_upper_limit;
             loop_validation(8, loop)  = c08_soh_lower_limit;
             loop_validation(9, loop)  = c09_soh_sequence;
             loop_validation(10, loop) = c10_soh_high_ev_use;
             loop_validation(11, loop) = c11_soh_low_ev_use;
             loop_validation(12, loop) = c12_soh_grid_only;
             loop_validation(13, loop) = c13_soh_disposal;
             loop_validation(14, loop) = c14_strategy_m_0;
             loop_validation(15, loop) = c15_strategy_m_1;
             loop_validation(16, loop) = c16_strategy_m_2;
             loop_validation(17, loop) = c17_strategy_m_3;
             loop_validation(18, loop) = c18_strategy_m_4;
             loop_validation(19, loop) = c19_strategy_m_5;
             loop_validation(20, loop) = c20_ev_low;
             loop_validation(21, loop) = c21_ev_high;
             loop_validation(22, loop) = c22_reliability_hev;
             loop_validation(23, loop) = c23_reliability_lev;
             loop_validation(24, loop) = c24_discharge_rel;             
             loop_validation(25, loop) = c25_elec_high;
             loop_validation(26, loop) = c26_elec_low;
             loop_validation(27, loop) = current_mode(i,t);
             loop_validation(28, loop) = transition_priority(loop);
             
            if (all_conditions == true) 
                dh_validation(i,t) = dh(t);
                dl_validation(i,t) = dl(t); 
                
                if (m_sl(i,t))
                    dl(t) = dl(t) - 1;                                                          
                end                
                if (m3_to_m0(i,t))
                    dl(t) = dl(t) + 1;     
                end

                if (m3_to_m1(i,t))
                    dl(t) = dl(t) + 1;     
                end
                
                if (m_sh(i,t))
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
                %(strcmp(load, 'duck'))
                if (strcmp(objective_type,'swap_only'))                  
                    objective(i,t) = objective(i,t-1) + ...
                       R_s*(m_s(i,t)) + ...
                       - E_l*(v(i,t)) - E_d*(m_r(i,t));                          
                elseif (strcmp(objective_type, 'lease_only'))
                    objective(i,t) = objective(i,t-1) + ...
                       R_s*(m_s(i,t)) + ...
                       R_l*(m_s(i,t) + m_e(i,t)) + ...
                       R_g(t)*(m_d(i,t) - m_c(i,t)) ...
                       - E_l*(v(i,t)) - E_d*(m_r(i,t));                                         
                elseif (strcmp(objective_type, 'swap_and_grid'))
                    objective(i,t) = objective(i,t-1) + ...
                       R_s*(m_s(i,t)) + ...
                       R_g(t)*(m_d(i,t) - m_c(i,t)) ...
                       - E_l*(v(i,t)) - E_d*(m_r(i,t));   
               elseif (strcmp(objective_type, 'swap_and_lease'))
                    objective(i,t) = objective(i,t-1) + ...
                       R_s*(m_s(i,t)) + ...
                       R_l*(m_s(i,t) + m_e(i,t)) + ...
                       - E_l*(v(i,t)) - E_d*(m_r(i,t));   
                elseif (strcmp(objective_type, 'lease_and_grid'))
                    objective(i,t) = objective(i,t-1) + ...
                       R_l*(m_s(i,t) + m_e(i,t)) + ...
                       R_g(t)*(m_d(i,t) - m_c(i,t)) ...
                       - E_l*(v(i,t)) - E_d*(m_r(i,t));                       
                elseif (strcmp(objective_type, 'swap_lease_and_grid'))
                    if (stage == 2) 
                        objective(i,t) = objective(i,t-1) + ...
                           R_s*(m_s(i,t)) + ...
                           R_l*(m_s(i,t) + m_e(i,t)) + ...
                           price_hour_15_years(t)*(m_d(i,t) - m_c(i,t)) ...
                           - E_l*(v(i,t)) - E_d*(m_r(i,t));    
                    else 
                        objective(i,t) = objective(i,t-1) + ...
                           R_s*(m_s(i,t)) + ...
                           R_l*(m_s(i,t) + m_e(i,t)) + ...
                           R_g(t)*(m_d(i,t) - m_c(i,t)) ...
                           - E_l*(v(i,t)) - E_d*(m_r(i,t)); 
                    end 
                end
            end
        end 

        state_space(t,1,i)  = current_mode(i,t);
        state_space(t,2,i)  = objective(i,t);
        state_space(t,3,i)  = current_transition(i,t);
        state_space(t,4,i)  = soc(i,t);
        state_space(t,5,i)  = soh(i,t);
    end     
end 


m_c_01 = sum(sum(m_c(:,1:24:T)));
m_c_02 = sum(sum(m_c(:,2:24:T)));
m_c_03 = sum(sum(m_c(:,3:24:T)));
m_c_04 = sum(sum(m_c(:,4:24:T)));
m_c_05 = sum(sum(m_c(:,5:24:T)));
m_c_06 = sum(sum(m_c(:,6:24:T)));
m_c_07 = sum(sum(m_c(:,7:24:T)));
m_c_08 = sum(sum(m_c(:,8:24:T)));
m_c_09 = sum(sum(m_c(:,9:24:T)));
m_c_10 = sum(sum(m_c(:,10:24:T)));
m_c_11 = sum(sum(m_c(:,11:24:T)));
m_c_12 = sum(sum(m_c(:,12:24:T)));
m_c_13 = sum(sum(m_c(:,13:24:T)));
m_c_14 = sum(sum(m_c(:,14:24:T)));
m_c_15 = sum(sum(m_c(:,15:24:T)));
m_c_16 = sum(sum(m_c(:,16:24:T)));
m_c_17 = sum(sum(m_c(:,17:24:T)));
m_c_18 = sum(sum(m_c(:,18:24:T)));
m_c_19 = sum(sum(m_c(:,19:24:T)));
m_c_20 = sum(sum(m_c(:,20:24:T)));
m_c_21 = sum(sum(m_c(:,21:24:T)));
m_c_22 = sum(sum(m_c(:,22:24:T)));
m_c_23 = sum(sum(m_c(:,23:24:T)));
m_c_24 = sum(sum(m_c(:,24:24:T)));

m_d_01 = sum(sum(m_d(:,1:24:T)));
m_d_02 = sum(sum(m_d(:,2:24:T)));
m_d_03 = sum(sum(m_d(:,3:24:T)));
m_d_04 = sum(sum(m_d(:,4:24:T)));
m_d_05 = sum(sum(m_d(:,5:24:T)));
m_d_06 = sum(sum(m_d(:,6:24:T)));
m_d_07 = sum(sum(m_d(:,7:24:T)));
m_d_08 = sum(sum(m_d(:,8:24:T)));
m_d_09 = sum(sum(m_d(:,9:24:T)));
m_d_10 = sum(sum(m_d(:,10:24:T)));
m_d_11 = sum(sum(m_d(:,11:24:T)));
m_d_12 = sum(sum(m_d(:,12:24:T)));
m_d_13 = sum(sum(m_d(:,13:24:T)));
m_d_14 = sum(sum(m_d(:,14:24:T)));
m_d_15 = sum(sum(m_d(:,15:24:T)));
m_d_16 = sum(sum(m_d(:,16:24:T)));
m_d_17 = sum(sum(m_d(:,17:24:T)));
m_d_18 = sum(sum(m_d(:,18:24:T)));
m_d_19 = sum(sum(m_d(:,19:24:T)));
m_d_20 = sum(sum(m_d(:,20:24:T)));
m_d_21 = sum(sum(m_d(:,21:24:T)));
m_d_22 = sum(sum(m_d(:,22:24:T)));
m_d_23 = sum(sum(m_d(:,23:24:T)));
m_d_24 = sum(sum(m_d(:,24:24:T)));

if (strcmp(load, 'duck')) 
    m_c_high    = m_c_19 + m_c_20 + m_c_21 + m_c_22;     
    m_c_mid     = m_c_01 + m_c_02 + m_c_03 + m_c_04 + ...
                  m_c_05 + m_c_06 + m_c_07 + m_c_08 + ...
                  m_c_09 + m_c_10 + m_c_11 + m_c_16 + ...
                  m_c_17 + m_c_18 + m_c_23 + m_c_24;
    m_c_low     = m_c_12 + m_c_13 + m_c_14 + m_c_15; 
    
    m_d_high    = m_d_19 + m_d_20 + m_d_21 + m_d_22;     
    m_d_mid     = m_d_01 + m_d_02 + m_d_03 + m_d_04 + ...
                  m_d_05 + m_d_06 + m_d_07 + m_d_08 + ...
                  m_d_09 + m_d_10 + m_d_11 + m_d_16 + ...
                  m_d_17 + m_d_18 + m_d_23 + m_d_24;
    m_d_low     = m_d_12 + m_d_13 + m_d_14 + m_d_15;    
elseif (strcmp(load, 'nswd')) 
    m_c_high    = m_c_17 + m_c_18 + m_c_19 + m_c_20;     
    m_c_mid     = m_c_01 + m_c_06 + m_c_07 + m_c_08 + ...
                  m_c_09 + m_c_10 + m_c_11 + m_c_12 + ...
                  m_c_13 + m_c_14 + m_c_15 + m_c_16 + ...                      
                  m_c_21 + m_c_22 + m_c_23 + m_c_24;
    m_c_low     = m_c_02 + m_c_03 + m_c_04 + m_c_05; 
    
    m_d_high    = m_d_17 + m_d_18 + m_d_19 + m_d_20;     
    m_d_mid     = m_d_01 + m_d_06 + m_d_07 + m_d_08 + ...
                  m_d_09 + m_d_10 + m_d_11 + m_d_12 + ...
                  m_d_13 + m_d_14 + m_d_15 + m_d_16 + ...                      
                  m_d_21 + m_d_22 + m_d_23 + m_d_24;
    m_d_low     = m_d_02 + m_d_03 + m_d_04 + m_d_05;          
elseif (strcmp(load, 'mojo'))
    m_c_high    = m_c_15 + m_c_16 + m_c_17 + m_c_28;     
    m_c_mid     = m_c_04 + m_c_05 + m_c_06 + m_c_07 + ...
                  m_c_08 + m_c_09 + m_c_10 + m_c_12 + ...
                  m_c_13 + m_c_13 + m_c_14 + m_c_19 + ...                      
                  m_c_20 + m_c_21 + m_c_22 + m_c_23;
    m_c_low     = m_c_24 + m_c_01 + m_c_02 + m_c_03;     
else
    error('no type of load existing'); 
end

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

if (stage == 0) 
    FileName = sprintf('%s_365_bss_%s_res_%d_seed_%d.xlsx', ...
            datestr(now, 'yymmdd_HHMM'), objective_type, ...
            reserved_batt, test_random_number);
else
    FileName = ...
    sprintf('%s_bss_stage_0%d_evcust_0%d_grid_%d_res_%d_seed_%d.xlsx',  ...
    datestr(now, 'yymmdd_HHMM'), stage, ...
    length(h_customer),has_smart_grid_services, ...
    reserved_batt, test_random_number);    
end
f_mx.filename                       = FileName;      
f_mx.objective_type                 = objective_type;
f_mx.objective                      = sum(objective(:, T));
f_mx.objective_revised              = sum(objective(:, T)) - E_l*hct_x_sum;
f_mx.seed                           = test_random_number;
f_mx.has_grid_scheduling            = has_smart_grid_services;
f_mx.reservation                    = reserved_batt;
f_mx.ev_customers                   = length(h_customer); 
f_mx.battery_count                  = I;
f_mx.sum_m_s                        = sum(sum(m_s));
f_mx.sum_m_e                        = sum(sum(m_e));
f_mx.sum_v                          = hct_x_sum;
f_mx.sum_m_r                        = sum(sum(m_r));
f_mx.sum_m_c                        = sum(sum(m_c));
f_mx.sum_m_c_high                   = sum(m_c_high);
f_mx.sum_m_c_mid                    = sum(m_c_mid);
f_mx.sum_m_c_low                    = sum(m_c_low);
f_mx.sum_m_d                        = sum(sum(m_d));
f_mx.sum_m_d_high                   = sum(m_d_high);
f_mx.sum_m_d_mid                    = sum(m_d_mid);
f_mx.sum_m_d_low                    = sum(m_d_low);    
 

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

