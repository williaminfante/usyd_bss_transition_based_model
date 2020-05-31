%Filename   : bss_pur_grid_1_res_1_sa_test.m
%Description: 
% 
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2017-03-29  1.5   Creation
%======================================================================

%-----------
default_values; 
%-----------
test_customer_data      = [];
has_smart_grid_services = 1;
reserved_batt           = 1; 
test_load_type          = 'sa';
%-----------
bss_thread; 
%-----------