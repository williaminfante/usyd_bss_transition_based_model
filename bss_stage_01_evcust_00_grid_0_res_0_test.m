%Filename   : bss_stage_01_evcust_00_grid_0_res_0_test.m
%Description: Test for no customers, NO smart grid services
% 
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2017-01-13  1.3   Creation
%======================================================================

%-----------
default_values; 
%-----------

test_customer_data      = [];
has_smart_grid_services = 0;
stage                   = 1; 

%-----------
bss_thread; 
%-----------