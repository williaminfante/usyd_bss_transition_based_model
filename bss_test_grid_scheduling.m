%Filename   : bss_test_grid_scheduling.m
%Description: grid scheduling with different seeds
% 
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-29  1.1   Creation
%william         2016-12-30  1.2   Added Reservation Tests
%======================================================================

test_random_number      = 05; 
bss_evcust_00_smartgrid_0_test
bss_evcust_00_smartgrid_1_test
bss_evcust_00_smartgrid_1_res_1_test

test_random_number      = 29; 
bss_evcust_00_smartgrid_0_test
bss_evcust_00_smartgrid_1_test
bss_evcust_00_smartgrid_1_res_1_test

test_random_number      = 88; 
bss_evcust_00_smartgrid_0_test
bss_evcust_00_smartgrid_1_test
bss_evcust_00_smartgrid_1_res_1_test