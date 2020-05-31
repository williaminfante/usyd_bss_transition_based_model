%Filename   : bss_stage_02_evcust_02_grid_1_res_1_test.m
%Description: 
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

objective_type          = 'swap_lease_and_grid';
has_smart_grid_services = 1;
reserved_batt           = 1; 
stage                   = 2;

%-----------
bss_thread; 
%-----------