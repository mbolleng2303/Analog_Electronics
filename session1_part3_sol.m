%% Analog Electronics Session 2: Common source amplifiers
% In this exersice session we will talk about common source amplifiers.
% These are amplifiers with the AC ground at the source of the transistor.
% The first exercise is a CS amplifier with a resistive load. In a second
% exercise this load is replaced by a PMOS transistors.
% The last exercise contains an amplifier with a cascode. Fun fun fun!
% TIP: The key ideas of this lab can also be found in the course notes (4.2
% Common-source stage)
% TIP: If your code does not work in matlab, you can ask Matlab why it is
% not working by typing 'why' in the command window.

%% Adding paths + Loading MOS tables
addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));

clear;
close all;
clc;

load ('UMC65_RVT.mat');

%% Initialize everything
designkitName   = 'umc65';
circuitTitle    = 'Analog Design - Session 2';

%Declaration of the circuit components
elementList.nmos = {'Mn1','Mn2'};
elementList.pmos = {'Mp2'};

spec.VDD        = 1.1;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator       ='spectre';
simulFile       = 0;
simulSkelFile   = 0;
analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice,...
    designkitName, NRVT, PRVT, simulator, simulFile, simulSkelFile);

analog          = cirCheckInChoice(analog, choice);

fprintf('\n--- First Exercise: Common source amplifier with resistor ---\n');
% %% EX1: Circuit
% disp('       VDD           ');
% disp('        |            ');
% disp('        RL           ');
% disp('        |----+-OUT   ');
% disp('  IN---Mn1   |       ');
% disp('        |    CL      ');
% disp('        |    |       ');
% disp('       GND  GND      ');
% 
% 
% %% EX1: Specs
% spec.fGBW       = 100e6;    % [Hz] GBW frequency
% spec.Cl         = 1e-12;   % [F] load capacitance
% 
% %% EX1: Goal (ETA: 45min)
% % Get at least a gain of 18dB (width of transistor < 1mm)
% % Notice that at the end mosCheckSaturation is used to check whether the
% % transistor of interest is biased in saturation!
% 
% %% EX1: Design Choices
% Mn1.lg  = 200e-9;           % [m] Design choice
% Mn1.vov = 0.0;             % [V] Design choice: Choice of inversion
% spec.VDD= 1.1;              % [V] Power supply voltage (1.1 -> 2.1V)
% VDS = spec.VDD/2;
% %% EX1: Implementation, Don't touch!
% Mn1.vsb = 0;
% Mn1.vds = VDS;
% Mn1.vth = tableValueWref('vth', NRVT, Mn1.lg, 0.3, Mn1.vds, Mn1.vsb);
% Mn1.vgs = Mn1.vth + Mn1.vov;
% 
% Mn1.gm = 2*pi*spec.fGBW*spec.Cl*1.1;
% Mn1.w = mosWidth('gm', Mn1.gm, Mn1);
% 
% Mn1 = mosNfingers(Mn1);
% Mn1 = mosOpValues(Mn1);
% 
% RL = (spec.VDD - Mn1.vds)/Mn1.ids;
% 
% fprintf('\nSolution of exercise 1:\n');
% Gain = Mn1.gm/(Mn1.gds+1/RL);
% Gain_dB = 20*log10(Gain);
% pole    = -(Mn1.gds+1/RL)/(spec.Cl+Mn1.cdd+Mn1.cdb);
% freq    = logspace(1,10,1e3);
% GBW = Gain*(abs(pole)/2/pi);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% if control toolbox in Matlab is available
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s       = tf('s');
% TF1     = Gain/(1-s/pole);
% figure;
% bode(TF1,freq); grid on;
% title('Voltage gain (common source with resistor load)');
% 
% %%%%%%%
% % else
% %%%%%%%
% 
% % TF1 = Gain./(1-1j*2*pi*freq/pole);
% % figure;
% % subplot(211);
% % semilogx(freq, 20*log10(abs(TF1)));
% % ylabel('Magnitude [dB]');
% % xlabel('Frequency [Hz]');
% % title('Bode plot');
% % subplot(212);
% % semilogx(freq, (angle(TF1))*180/pi);
% % ylabel('Phase [degree]');
% % xlabel('Frequency [Hz]');
% 
% fprintf('RL      = %.2f Ohm \n', RL);
% fprintf('rds     = %.2f Ohm \n', 1/Mn1.gds);
% fprintf('Gain    = %.2f \n', Gain);
% fprintf('Gain_dB = %.2f dB \n', Gain_dB);
% fprintf('Vds     = %.2f V \n', Mn1.vds);
% fprintf('Vdsat   = %.2f V \n', Mn1.vdsat);
% fprintf('width   = %.2f mm \n', Mn1.w*1000);
% fprintf('GBW     = % 2.2f MHZ \n ', GBW/1e6);
% 
% fprintf('\nGain(in dB) requirement > 18 dB:\n');
% if (Gain_dB > 18)
%    fprintf('your design: Success!\n') % Check Gain requirement
% else
%    fprintf('your design: Fail...\n')
% end
% 
% fprintf('\nTransisors in saturation:\n');
% if mosCheckSaturation(Mn1)
%    fprintf('Mn1: Success!\n') % Check if transistor is in saturation
% end
% 
% fprintf('\nWidth requirement < 1 mm:\n');
% if (Mn1.w < 0.001)
%    fprintf('Mn1: Success!\n') % Check width requirement
% else
%    fprintf('Mn1: Fail...\n')
% end

fprintf('\n--- Second Exercise: Common source amplifier with PMOS ---\n');
%% EX2: Circuit
disp('       VDD           ');
disp('        |            ');
disp('       Mp2           ');
disp('        |----+-OUT   ');
disp('  IN---Mn2   |       ');
disp('        |    CL      ');
disp('        |    |       ');
disp('       GND  GND      ');


%% EX2: Specs
spec.fGBW       = 100e6;    % [Hz] GBW frequency
spec.Cl         = 1e-12;   % [F] load capacitance
spec.VDD        = 1.1;      % [V] Power supply voltage

%% EX2: Goal (ETA: 45min)
% Get at least a gain of 26dB (width of transistor < 1mm)

%% EX2: Design Choices
Mn2.lg  = 200e-9;           % Design choice: flatness
Mn2.vov = 0.0;            % Design choice: inversion level
Mn2.vds = spec.VDD/2;       % Design choice: maximum swing
Mp2.lg  = 500e-9;            % Design choice: flatness
Mp2.vov = 0.0;             % Design choice: inversion level
Mp2.vds = -spec.VDD/2;

%% EX2: Implementation, Don't touch!
Mn2.vsb = 0;
Mn2.vth = tableValueWref('vth', NRVT, Mn2.lg, 0.3,Mn2.vds, Mn2.vsb);
Mn2.vgs = Mn2.vth+Mn2.vov;
Mn2.gm = 2*pi*spec.fGBW*spec.Cl*1.05;
Mn2.w = mosWidth('gm', Mn2.gm, Mn2);

Mn2 = mosNfingers(Mn2);
Mn2 = mosOpValues(Mn2);

Mp2.vsb = 0;
Mp2.vth = tableValueWref('vth', PRVT, Mp2.lg, -0.3, Mp2.vds, Mp2.vsb);
Mp2.vgs = Mp2.vov + Mp2.vth;
Mp2.ids = Mn2.ids;
Mp2.w = mosWidth('ids', Mp2.ids, Mp2);

Mp2 = mosNfingers(Mp2);
Mp2 = mosOpValues(Mp2);

Gain    = Mn2.gm/(Mn2.gds+Mp2.gds);
Gain_dB = 20*log10(Gain);
pole    = -(Mn2.gds+Mp2.gds)/(spec.Cl+Mn2.cdd+Mn2.cdb+Mp2.cdd+Mp2.cdb);
freq    = logspace(1,10,1e3);
GBW = Gain*(abs(pole)/2/pi);

%% EX2: Figures of Merit + plot

fprintf('\nSolution of exercise 2:\n');
fprintf('rds_n = %6.2f Ohm\n', 1/Mn2.gds);
fprintf('rds_p = %6.2f Ohm\n', 1/Mp2.gds);
fprintf('Gain = %.2f\n',Gain);
fprintf('Gain_dB = %.2f dB\n',Gain_dB);
fprintf('width of Mp2   = %.2f mm \n', Mp2.w*1000);
fprintf('width of Mn2   = %.2f mm \n', Mn2.w*1000);
fprintf('GBW            = % 2.2f MHz \n', GBW/1e6);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if control toolbox in Matlab is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s       = tf('s');
TF2     = Gain/(1-s/pole);
figure;
bode(TF2,freq); grid on;
title('Voltage gain (Common source with PMOS load)');

%%%%%%%%
%% else
%%%%%%%%
% TF2 = Gain./(1-1j*2*pi*freq/pole);
% figure;
% subplot(211);
% semilogx(freq, 20*log10(abs(TF2)));
% ylabel('Magnitude [dB]');
% xlabel('Frequency [Hz]');
% title('Bode plot');
% subplot(212);
% semilogx(freq, (angle(TF2))*180/pi);
% ylabel('Phase [degree]');
% xlabel('Frequency [Hz]');


fprintf('\nGain(in dB) requirement > 26 dB:\n');
if (Gain_dB > 26)
    fprintf('your design: Success!\n') % Check Gain requirement
else
    fprintf('your design: Fail...\n')
end

fprintf('\nTransisors in saturation:\n');
if mosCheckSaturation(Mp2)
    fprintf('Mp2: Success!\n') % Check if transistor is in saturation
else
    fprintf('Mp2: Fail...\n')
end
if mosCheckSaturation(Mn2)
    fprintf('Mn2: Success!\n')
else
    fprintf('Mn2: Fail...\n')
end

fprintf('\nWidth requirement < 1 mm:\n');
if (Mp2.w < 0.001)
    fprintf('Mp2: Success!\n') % Check width requirement
else
    fprintf('Mp2: Fail...\n')
end
if (Mn2.w < 0.001)
    fprintf('Mn2: Success!\n')
else
    fprintf('Mn2: Fail...\n')
end

