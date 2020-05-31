%Filename   : bss_stage_02_evcust_00_grid_1_res_0_test.m
%Description: Test for no customers, has smart grid services
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
has_smart_grid_services = 1;
stage                   = 2;

%-----------
bss_thread; 
%-----------