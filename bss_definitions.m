%Filename   : bss_definitions.m
%Description: Definition of names based on transitions
%
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-06  1.0   Creation
%william         2017-03-29  1.5   Removed m3 state
%======================================================================

%INITIALIZE
%Transitions
m0_to_m0 = ones(I,T);    %t = 1, state is waiting
m1_to_m0 = zeros(I, T); 
m2_to_m0 = zeros(I, T);
m4_to_m0 = zeros(I, T);
m0_to_m1 = zeros(I, T);
m1_to_m1 = zeros(I, T);
m2_to_m1 = zeros(I, T);
m4_to_m1 = zeros(I, T);
m0_to_m2 = zeros(I, T);
m1_to_m2 = zeros(I, T);
m2_to_m2 = zeros(I, T);
m0_to_m4 = zeros(I, T);
m1_to_m4 = zeros(I, T);
m2_to_m4 = zeros(I, T);
m4_to_m4 = zeros(I, T);
m0_to_m5 = zeros(I, T);
m1_to_m5 = zeros(I, T);
m2_to_m5 = zeros(I, T);
m5_to_m5 = zeros(I, T);


%Definitions 
%Modes
m0 = m0_to_m0 + m1_to_m0 + m2_to_m0 + m4_to_m0;             %Waiting
m1 = m0_to_m1 + m1_to_m1 + m2_to_m1 + m4_to_m1;             %Charge Grid
m2 = m0_to_m2 + m1_to_m2 + m2_to_m2;                        %Discharge Grid
m4 = m0_to_m4 + m1_to_m4 + m2_to_m4            + m4_to_m4;  %Discharge HEV
m5 = m0_to_m5 + m1_to_m5 + m2_to_m5 + m5_to_m5;             %Disposal

%States
b0      = m0_to_m0 + m1_to_m0 + m2_to_m0 + ...              %Idle State
          + m4_to_m0 + ...
          m0_to_m5 + m1_to_m5 + m2_to_m5 + m5_to_m5; 
bc      = m0_to_m1 + m1_to_m1 + m2_to_m1 + ...              %Charging State
          m4_to_m1;
bd      = m0_to_m2 + m1_to_m2 + m2_to_m2 + ...              %Discharging
          m0_to_m4 + m1_to_m4 + m2_to_m4 + m4_to_m4;

%Functions
m_w     = m0_to_m0 + m1_to_m0 + m2_to_m0 + ... %Waiting
          m4_to_m0; 
m_c     = m0_to_m1 + m1_to_m1 + m2_to_m1 + ... %Charge to grid
          m4_to_m1;           
m_d     = m0_to_m2 + m1_to_m2 + m2_to_m2;      %Discharge to grid
m_s     = m0_to_m4 + m1_to_m4 + m2_to_m4;      %ev Swapping 
m_e     = m4_to_m4;                            %Ev driving
m_r     = m0_to_m5 + m1_to_m5 + m2_to_m5;      %Recycle
m_f     = m5_to_m5;                            %Final state 

%Shortcut Definitions 
m_b     = m0 + m1 + m2;                        %located in the Bss
m_return    = m4_to_m0 + m4_to_m1;             %Return to bss from High use 