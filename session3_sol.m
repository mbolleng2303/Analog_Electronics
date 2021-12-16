%% Analog Electronics Session 3: design of differential pairs and of an OTA
% In this MATLAB session we will verify the results of exercise 2 from 
% session 3. 

%% Adding paths + Loading MOS tables
addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));

clear;
%close all;
clc;

load ('UMC65_RVT.mat');

%% Initialization
designkitName   = 'umc65';
circuitTitle    = 'Analog Design - Session 3';

%Declaration of the circuit components
elementList.nmos = {'M1','M2','M4'};
elementList.pmos = {'M5','M6'};

spec.VDD   = 1.1;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator             = 'spectre';
simulFile             = 0;
simulSkelFile         = 0;

analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice,...
    designkitName, NRVT, PRVT, simulator, simulFile, simulSkelFile);

analog          = cirCheckInChoice(analog, choice);

fprintf('\n----------- OTA -----------------------\n');
disp('                                      ');
disp('       VDD         VDD                ');
disp('        |           |                 ');
disp('       M5---+-------M6                ');
disp('        |   |       |                 ');
disp('        +---+       +----------+-OUT2 ');
disp('        |           |          |      ');
disp(' IN1---M1           M2---IN2   |      ');
disp('        +-----+-----+          |      ');
disp('              |                CL     ');
disp('              M4               |      ');
disp('              |                |      ');
disp('             GND              GND     ');
fprintf('\n---------------------------------------\n');

%% AI: Set the given Specs
spec.VDD   = 1.1;       % [V] Power supply voltage
spec.fGBW  = 716e6;     % [Hz] GBW frequency
spec.CL    = 500e-15;   % [F] load capacitance
spec.gain  = 8;         % [ ] voltage gain

%% AI: Fill out the empty spaces
Voutcm   = 1*0.5; % option 1: changing Vout,cm changes VOV, which changes pole-zero-doublet-location
M4.ids   = 400e-6;
M1.vov   = -0.2;             
M5.vgs   = Voutcm-spec.VDD; 

%%
% AI: Code the OTA design. The comments below should help you, how to 
%     type the code

%% Section for setting up M4
%  Set the length
M4.lg  = 1e-6;

%  set the bias conditions
M4.vov = 0.2;
%  Reminder: M4.ids see "Design choices"

M4.vsb = 0;    % triple-well-technology
M4.vds = 0.11; % or M4.vov as pessimistic estimate 
M4.vth = tableValueWref('vth',NRVT,M4.lg, 0,M4.vds,M4.vsb);
M4.vgs = M4.vth + M4.vov;

% set the width
M4.w = mosWidth('ids', M4.ids, M4);
M4   = mosNfingers(M4);

% calculate remaining OP-parameters (fT, Cgs, etc.)
M4 = mosOpValues(M4);

%% Section for setting up M5 and M6
%  set length
M5.lg  = 1*0.2e-6; % option 2: changing L has an impact on Vth, which changes VOV and hence pole-zero-doublet-location

%  set the bias conditions
%  Reminder: M5.vgs see "Design choices"
M5.ids = M4.ids/2;

M5.vsb = 0; % triple-well-technology
M5.vds = M5.vgs;
M5.vth = tableValueWref('vth',PRVT,M5.lg,M5.vgs,M5.vds,M5.vsb);
M5.vov = M5.vgs-M5.vth;

%  set the width
M5.w = mosWidth('ids',M5.ids,M5);
M5   = mosNfingers(M5); %LTSpice

%  calculate remaining OP-parameters (fT, Cgs, etc.)
M5 = mosOpValues(M5);

%  copy M5 to M6
M6 = cirElementCopy(M5, M6);

%% Section for setting up  M1 and M2
%  set the length
M1.lg = 65e-9;

%  set the bias conditions
%  M1.vov see above
M1.ids = M4.ids/2;

M1.vsb = 0; % triple-well-technology
M1.vds = spec.VDD+M5.vgs-M4.vds;
M1.vth = tableValueWref('vth', NRVT,M1.lg, 0, M1.vds,M1.vsb);
M1.vgs = M1.vov + M1.vth;

%  set the width
M1.w = mosWidth('ids',M1.ids, M1);
M1 = mosNfingers(M1);

%  calculate remaining OP-parameters (fT, Cgs, etc.)
M1 = mosOpValues(M1);

%  copy M1 to M2
M2 = cirElementCopy(M1, M2);

%% AI: Fill out the empty equations with the correct formula

AvDC = M2.gm/(M2.gds+M6.gds)    % DC gain
C1  = M5.cgs+M6.cgs+M5.cdb;     % Capacitance on node 1
G1   = M5.gm;                   % Resistance  on node 1
C2  = spec.CL;                  % Capacitance on node 2
G2   = M6.gds + M2.gds;         % Resistance  on node 2

Voutcm      = spec.VDD + M5.vgs
Vincm       = M4.vds + M1.vgs
SpecGain    = 20*log10(spec.gain)
pdom        = G2/(C2)
fGBW        = M1.gm*2*pi*spec.CL % [Hz]    GBW frequency
SR          = M4.ids/spec.CL/1e9 % [mV/ns] slew rate
Vcm_in_min  = M4.vdsat + M1.vgs
Vcm_in_max  = spec.VDD + M5.vgs - M1.vdsat + M1.vgs
Vout_cm_min = M4.vds + M1.vdsat
Vout_cm_max = spec.VDD + M5.vdsat

% Sanitycheck, if everything is in saturation
fprintf('\nTransistors in saturation:\n');
if mosCheckSaturation(M1)
    fprintf('M1:Success\n')
end
if mosCheckSaturation(M2)
    fprintf('M2:Success\n')
end
if mosCheckSaturation(M4)
    fprintf('M4:Success\n')
end
if mosCheckSaturation(M5)
    fprintf('M5:Success\n')
end
if mosCheckSaturation(M6)
    fprintf('M6:Success\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if control toolbox in Matlab is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s   = tf('s');
TF1 = AvDC*(1+s*C1/(2*G1))/((1+s*C1/G1)*(1+s*C2/G2));
freq = logspace(1,10,1e3);
figure(1)
bode(TF1,2*pi*freq); grid on;
h = gcr;
setoptions(h,'FreqUnits','Hz');
hold all
title('Frequency response OTA');
figure(2)
TFfb=feedback(TF1,1);
step(TFfb);
title('Step response OTA');
hold all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if not, uncomment the code below and commet the control toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  freq = logspace(1,10,1e3);
%  Av   = AvDC .*...
%         (1+1j*2*pi*freq.*(Cd1/(2*G1))) .* ...
%         (1./(1+1j*2*pi*freq.*(Cd1/(G1)))) .* ...
%         (1./(1+1j*2*pi*freq.*(Cd2/(G2))));
%
% figure
% subplot(211);
% semilogx(freq, 20*log10(abs(Av)));
% ylabel('Magnitude [dB]');
% xlabel('Frequency [Hz]');
% title('Frequency response OTA');
% grid on;
% hold all;
% subplot(212);
% semilogx(freq,(angle(Av))*180/pi);
% ylabel('Phase [degree]');
% xlabel('Frequency [Hz]');
% grid on;
% hold all;

