%Filename   : default_values.m
%Description: Default values for common test parameters
%
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-06  1.0   Creation
%william         2016-12-29  1.1   Added reserved_batt
%                                  Corrected random number generation
%william         2017-01-13  1.3   Added stage for NSW Electricity
%william         2017-03-29  1.5   Added hour_slice
%======================================================================
hour_slice   = 2; %1 - per hour, 2 - per 30 minutes, 12 - per 5 minutes
test_customer_data      = [7, 17]*hour_slice;  %(1 start at 7 am, 
                                               % 1 start at 5 pm )
test_battery_number     = 10; 
test_load_type          = 'nswd';   %[nswd, or sa]
test_tot_time_in_hr     = 131400;   %15 years
objective_type          = 'swap_lease_and_grid';

%flags
has_smart_grid_services = 0;
reserved_batt           = 0; 
stage                   = 0; 

%constants
HOURS_IN_A_DAY          = 24; 
DAYS_IN_A_YEAR          = 365; 