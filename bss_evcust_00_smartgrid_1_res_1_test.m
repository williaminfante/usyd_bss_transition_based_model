%Filename   : bss_evcust_00_smartgrid_1_res_1_test.m
%Description: Test for no customers, has smart grid services
% 
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-30  1.2   Creation
%======================================================================

%-----------
default_values; 
%-----------

test_customer_data      = [];
has_smart_grid_services = 1;
reserved_batt           = 1; 

%-----------
bss_thread; 
%-----------