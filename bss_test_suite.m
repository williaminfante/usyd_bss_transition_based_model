%Filename   : bss_test_suite.m
%Description: 
% 
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-29  1.1   Creation
%======================================================================

test_random_number      = 99; 

%no main control 
bss_obj_type_swap_only_test
bss_obj_type_lease_only_test
bss_obj_type_swap_and_lease_test

%smart grid 
bss_obj_type_lease_and_grid_test
bss_obj_type_swap_and_grid_test
bss_obj_type_swap_lease_and_grid_test

%consider customer loyalty 
bss_obj_type_res_lease_and_grid_test
bss_obj_type_res_swap_and_grid_test
bss_obj_type_res_swap_lease_and_grid_test