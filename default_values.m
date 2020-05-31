%Filename   : default_values.m
%Description: Default values for common test parameters
%
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-06  1.0   Creation
%======================================================================

test_random_number      = 1; 
test_customer_data      = [7, 17];  %(1 start at 7 am, 
                                    % 1 start at 5 pm )
test_battery_number     = 10; 
test_load_type          = 'nswd';   %NSW electricity demand
test_tot_time_in_hr     = 131400;   %15 years
objective_type          = 'swap_lease_and_grid';
%flags
has_smart_grid_services = 0;