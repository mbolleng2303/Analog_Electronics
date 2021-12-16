%% Analog Electronics final project

clc; 
close all; 
clear;

addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));

load('UMC65_RVT.mat');

%% Initializations
designkitName		= 'umc65';
circuitTitle		= 'Analog Design - Project';
elementList.nmos	= {'Mn3','Mn4','Mn6'};
elementList.pmos	= {'Mp1','Mp2','Mp5','Mp7','Mp8'};

% spec.VDD            = 1.1; 
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;

simulator			= 'spectre';
simulFile			= 0;
simulSkelFile		= 0;
spec				= [];
analog			= cirInit('analog',circuitTitle,'top',elementList,spec,choice,...
						designkitName,NRVT,PRVT,simulator,simulFile,simulSkelFile);
analog			= cirCheckInChoice(analog, choice);


%% Project: circuit
disp('                                                       ');
disp('  VDD          VDD                    VDD              ');
disp('   |             |                      |              ');
disp('  Mp8-+---------Mp7---------------------Mp5            ');
disp('   |--+          |                      |              ');
disp('   |          +--+--+         node 3->  +-----+---OUT  ');
disp('   |          |     |                   |     |        ');
disp('   |    IN1--Mp1   Mp2--IN2             |     |        ');
disp('   |          |     |                   |     |        ');
disp('   | node 1-> |     | <-node 2          |     Cl       ');
disp('   |          |--+  +------+-Cm---Rm----+     |        ');
disp(' Ibias        |  |  |      |    ↑       |     |        ');
disp('   |         Mn3-+-Mn4     |  node 4    |     |        ');
disp('   |          |     |      +-----------Mn6    |        ');
disp('   |          |     |                   |     |        ');
disp('  GND        GND   GND                 GND   GND       ');

%% AI: Implement your OpAmp according to your handcalculations
%% EX: Specs 

spec.fGBW   = 100e6;        % [Hz] GBW frequency 

spec.GBW    = spec.fGBW * 2 * pi; % radians/s 

spec.gain   =10^( 50/20);          % [] voltage gain 

spec.gaindB = 50; %  

spec.CL     = 2e-12;       % [F], Load cap 

spec.RL     = 0;         % [Ohm] 

Vout = 0.52; 

Voutcm = 0.3; 

Mp7.vdsat = -0.18; 

Mp5.vdsat = -0.18; 

%% Miller Capacitance and resistance 

spec.Cm     = spec.CL/4;    % [F], miller cap 

spec.Rm     = 346;     % [F], Load cap 

  

%% Mn6: Second stage NMOS: Design choices and implementation 

Mn6.gm = 2*pi * spec.fGBW * spec.CL; 

   % [V], DC output voltage 

  

Mn6.lg  = 200e-9;       % [m], channel length 

Mn6.vsb = 0; 

Mn6.vds = 0.52; 

Mn6.vth = tableValueWref('vth', NRVT, Mn6.lg, 0, Mn6.vds, Mn6.vsb); 

Mn6.vgs = Voutcm-0.1; 

Mn6.vov = Mn6.vgs- Mn6.vth ; 

Mn6.w = mosWidth('gm', Mn6.gm, Mn6); 

Mn6     = mosNfingers(Mn6); 

Mn6     = mosOpValues(Mn6); 

%% Mp5: Second stage NMOS: Design choices and implementation 

Mp5.lg  = 1000e-9; 

Mp5.vsb = 0; 

Mp5.vov = -0.2; 

Mp5.vds = Mn6.vds - spec.VDD; 

Mp5.ids = Mn6.ids; 

Mp5.vth = tableValueWref('vth', PRVT, Mp5.lg, 0, Mp5.vds, Mp5.vsb); 

Mp5.vgs = Mp5.vov + Mp5.vth; 

Mp5.w = mosWidth('ids', Mp5.ids, Mp5); 

Mp5     = mosNfingers(Mp5); 

Mp5     = mosOpValues(Mp5); 

  

%% Mn1/2: Second stage NMOS: Design choices and implementation 

Mp2.lg = 100e-9; %changer de 65 à 100 

Mp2.gm = (spec.Cm)*spec.GBW; 

Mp2.vov = 0.165;  

Mp2.vsb = 0; % triple-well-technology 

Mp2.vds = Mn6.vgs - (Mp7.vdsat) - spec.VDD; 

Mp2.vth = tableValueWref('vth', PRVT,Mp2.lg, 0, Mp2.vds,Mp2.vsb); 

Mp2.vgs = Mp2.vov + Mp2.vth; 

Mp2.w = mosWidth('gm',Mp2.gm, Mp2); 

Mp2 = mosNfingers(Mp2); 

Mp2 = mosOpValues(Mp2); 

Mp1 = cirElementCopy(Mp2, Mp1); 

%% Mn4/3: Second stage NMOS: Design choices and implementation 

Mn4.lg  = 200e-9;      % [m], channel length CHANGER LG DE MN4 POUR AUGMENETER FDANS LE GRAPHE 

Mn4.vsb = 0.0; 

Mn4.ids = Mp2.ids; 

Mn4.vds = Mn6.vgs; 

Mn4.vgs = Mn6.vgs; 

Mn4.vth = tableValueWref('vth', NRVT, Mn4.lg, 0, Mn4.vds, Mn4.vsb); 

Mn4.vov = Mn4.vgs - Mn4.vth ; 

Mn4.w   = mosWidth('ids', Mn4.ids, Mn4); 

Mn4     = mosNfingers(Mn4); 

Mn4     = mosOpValues(Mn4); 

Mn3 = cirElementCopy(Mn4, Mn3); 

%% Mn7/8: Second stage NMOS: Design choices and implementation 

Mp7.lg  = 1e-6; 

Mp7.ids = Mp2.ids*2; 

Mp7.vsb = 0;    % triple-well-technology 

Mp7.vds = -0.18;  

Mp7.vgs = Mp5.vgs; 

Mp7.vth = tableValueWref('vth',PRVT,Mp7.lg, 0,Mp7.vds,Mp7.vsb); 

Mp7.vov = -Mp7.vth +Mp7.vgs; 

  

Mp7.w = mosWidth('ids', Mp7.ids, Mp7); 

Mp7   = mosNfingers(Mp7); 

Mp7 = mosOpValues(Mp7); 

  

%% Mp8: Second stage NMOS: Design choices and implementation 

Mp8.lg  = 1e-6; 

Mp8.ids = 1e-3; 

Mp8.vov = -0.2; 

Mp8.vsb = 0;    % triple-well-technology 

Mp8.vgs = Mp7.vgs; 

Mp8.vds = Mp7.vgs;  

Mp8.vth = tableValueWref('vth',PRVT,Mp8.lg, 0,Mp8.vds,Mp8.vsb); 

Mp8.vov = -Mp8.vth +Mp8.vgs ; 

  

Mp8.w = mosWidth('ids', Mp8.ids, Mp8); 

Mp8   = mosNfingers(Mp8); 

Mp8 = mosOpValues(Mp8); 
%% AI: Fill out the empty variables required to plot the transfer-function.
%  meaning of each variable see comment and
%  location of nodes see line 31 

AvDC1 = ;  % DC gain 1st stage
AvDC2 = ;  % DC gain 2nd stage
C1    = ;  % Capacitance on node 1
G1    = ;  % Admittance  on node 1
C2    = ;  % Capacitance on node 2
G2    = ;  % Admittance  on node 2
C3    = ;  % Capacitance on node 3
G3    = ;  % Admittance  on node 3
C4    = ;  % Capacitance on node 4
G4    = ;  % Admittance on node 4 (hint: what happens with CL at very high frequencies?)


%% AI: Set-up Rm, Cc and CL and calculate the zero required for the transfer-fct

spec.Cm = ;
spec.Cl = ;
spec.Rm = ; 
z1 = ;

%% AI: Fill out the empty variables required for the performance summary
Vin_cm_min  = ;   
Vin_cm_max  = ; 
Vout_cm_min = ;                         
Vout_cm_max = ;             
Pdiss       = ;

%% Sanity check (do not modify)

disp('======================================');
disp('=      Transistors in saturation     =');
disp('======================================');
if mosCheckSaturation(Mp1)
	fprintf('\nMp1:Success\n')
end
if mosCheckSaturation(Mp2)
	fprintf('Mp2:Success\n')
end
if mosCheckSaturation(Mn3)
	fprintf('Mn3:Success\n')
end
if mosCheckSaturation(Mn4)
	fprintf('Mn4:Success\n')
end
if mosCheckSaturation(Mp5)
	fprintf('Mp5:Success\n')
end
if mosCheckSaturation(Mn6)
	fprintf('Mn6:Success\n')
end
if mosCheckSaturation(Mp7)
	fprintf('Mp7:Success\n')
end
if mosCheckSaturation(Mp8)
	fprintf('Mp8:Success\n\n')
end


%% Summary of sizes and biasing points (do not modify)

disp('======================================');
disp('=    Sizes and operating points      =');
disp('======================================');
analog = cirElementsCheckOut(analog); % Update circuit file with 
% transistor sizes
mosPrintSizesAndOpInfo(1,analog); % Print the sizes of the 
% transistors in the circuit file
fprintf('IBIAS\t= %6.2fmA\nRm\t= %6.2f Ohm\nCm\t= %6.2fpF\n\n',Mp8.ids/1e-3,spec.Rm,spec.Cm/1e-12);

%% Performance summary (do not modify)

disp('======================================');
disp('=        Performance                 =');
disp('======================================');

fprintf('\nmetrik        \t result\n');
fprintf('Vin,cm,min [mV] \t%.0f\n',Vin_cm_min/1e-3);
fprintf('Vin,cm,max [mV] \t%.0f\n',Vin_cm_max/1e-3);
fprintf('Vout,cm,min [mV] \t%.0f\n',Vout_cm_min/1e-3);
fprintf('Vout,cm,max [mV] \t%.0f\n',Vout_cm_max/1e-3);
fprintf('Pdiss [mW]       \t%.1f\n',Pdiss/1e-3);

%% Ploting transfer function (do not modify)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if control toolbox in Matlab is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s   = tf('s');
% transfer function
TF1   = AvDC1*AvDC2*((1+s*C1/(2*G1))*(1-s*(1/z1)))/ ...
                ((1+s*C1/G1)*(1+s*C2/G2)*(1+s*C3/G3)*(1+s*C4/G4));

freq = logspace(1,12,1e3);
figure(1)
bode(TF1,2*pi*freq); grid on;
h = gcr;
setoptions(h,'FreqUnits','Hz');
title('Frequency response Opamp');
hold all


