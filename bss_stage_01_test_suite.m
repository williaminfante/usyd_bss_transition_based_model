%Filename   : bss_stage_01_test_suite.m
%Description: 
% 
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2017-01-13  1.3   Creation
%======================================================================

%Scenario 0 - No Electric Vehicles 
test_random_number      = 88; 
%%No Scheduling
bss_stage_01_evcust_00_grid_0_res_0_test
%%Strict Scheduling
bss_stage_01_evcust_00_grid_1_res_0_test
%%Scheduling with Battery Reservation
bss_stage_01_evcust_00_grid_1_res_1_test

%Scenario 1 - With Electric Vehicles 
test_random_number      = 99; 
%%No Scheduling
bss_stage_01_evcust_02_grid_0_res_0_test

%%Strict Scheduling
bss_stage_01_evcust_02_grid_1_res_0_test

%%Scheduling with Battery Reservation
bss_stage_01_evcust_02_grid_1_res_1_test