%Filename   : bss_definitions.m
%Description: Definition of names based on transitions
%
%Modification History: 
%======================================================================
%Author          Date        Ver   Remarks  
%======================================================================
%william         2016-12-06  1.0   Creation
%======================================================================

%INITIALIZE
%Transitions
m0_to_m0 = ones(I,T);    %t = 1, state is waiting
m1_to_m0 = zeros(I, T); 
m2_to_m0 = zeros(I, T);
m3_to_m0 = zeros(I, T); 
m4_to_m0 = zeros(I, T);
m0_to_m1 = zeros(I, T);
m1_to_m1 = zeros(I, T);
m2_to_m1 = zeros(I, T);
m3_to_m1 = zeros(I, T);
m4_to_m1 = zeros(I, T);
m0_to_m2 = zeros(I, T);
m1_to_m2 = zeros(I, T);
m2_to_m2 = zeros(I, T);
m0_to_m3 = zeros(I, T);
m1_to_m3 = zeros(I, T);
m2_to_m3 = zeros(I, T);
m3_to_m3 = zeros(I, T);
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
m0 = m0_to_m0 + m1_to_m0 + m2_to_m0 + m3_to_m0 + m4_to_m0; %Waiting
m1 = m0_to_m1 + m1_to_m1 + m2_to_m1 + m3_to_m1 + m4_to_m1; %Charge Grid
m2 = m0_to_m2 + m1_to_m2 + m2_to_m2;                       %Discharge Grid
m3 = m0_to_m3 + m1_to_m3 + m2_to_m3 + m3_to_m3;            %Discharge LEV
m4 = m0_to_m4 + m1_to_m4 + m2_to_m4            + m4_to_m4; %Discharge HEV
m5 = m0_to_m5 + m1_to_m5 + m2_to_m5 + m5_to_m5;            %Disposal

%States
b0      = m0_to_m0 + m1_to_m0 + m2_to_m0 + ...              %Idle State
          m3_to_m0 + m4_to_m0 + ...
          m0_to_m5 + m1_to_m5 + m2_to_m5 + m5_to_m5; 
bc      = m0_to_m1 + m1_to_m1 + m2_to_m1 + m3_to_m1 + ...   %Charging State
          m4_to_m1;
bd      = m0_to_m2 + m1_to_m2 + m2_to_m2 + ...              %Discharging
          m0_to_m3 + m1_to_m3 + m2_to_m3 + m3_to_m3 + ...
          m0_to_m4 + m1_to_m4 + m2_to_m4 + m4_to_m4;

%Functions
m_w     = m0_to_m0 + m1_to_m0 + m2_to_m0 + ... %Waiting
          m3_to_m0 + m4_to_m0; 
m_c     = m0_to_m1 + m1_to_m1 + m2_to_m1 + ... %Charge to grid
          m3_to_m1 + m4_to_m1;           
m_d     = m0_to_m2 + m1_to_m2 + m2_to_m2;      %Discharge to grid
m_s     = m0_to_m3 + m1_to_m3 + m2_to_m3 + ... %Swapping 
          m0_to_m4 + m1_to_m4 + m2_to_m4;      
m_sl    = m0_to_m3 + m1_to_m3 + m2_to_m3;      %Swapping (Low ev use)
m_sh    = m0_to_m4 + m1_to_m4 + m2_to_m4;      %Swapping (High ev use)
m_e     = m3_to_m3 + m4_to_m4;                 %Ev driving
m_el    = m3_to_m3;                            %Ev driving (Low use)
m_eh    = m4_to_m4;                            %Ev driving (High use)
m_r     = m0_to_m5 + m1_to_m5 + m2_to_m5;      %Recycle
m_f     = m5_to_m5;                            %Final state 

%Shortcut Definitions 
m_b     = m0 + m1 + m2;                        %located in the Bss
m_rh    = m4_to_m0 + m4_to_m1;                 %Return to bss from High use 