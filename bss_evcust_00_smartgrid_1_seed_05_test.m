%Filename   : bss_evcust_00_smartgrid_1_seed_05_test.m
%Description: Test for no customers, has smart grid services, seed 05, 1 yr
% 
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-06  1.0   Creation
%======================================================================

%-----------
default_values; 
%-----------

test_random_number      = 05; 
test_customer_data      = [];
test_tot_time_in_hr     = 8760; %1 year
has_smart_grid_services = 1;

%-----------
bss_thread; 
%-----------