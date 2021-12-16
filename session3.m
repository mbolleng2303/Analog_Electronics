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

%-------------------------------------------------------------------------
%% AI: Set the given Specs

spec.VDD   = ;       % [V] Power supply voltage
spec.fGBW  = ;     % [Hz] GBW frequency
spec.CL    = ;   % [F] load capacitance
spec.gain  = ;         % [ ] voltage gain

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%% AI: Fill out the empty spaces

Voutcm   = ;
M4.ids   = ;
M1.vov   = ;             
M5.vgs   = ; 

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% AI: Code the OTA design. The comments below should help you, how to 
%     type the code

%% Section for setting up M4
%  Set the length



%  set the bias conditions



% set the width



% calculate remaining OP-parameters (fT, Cgs, etc.)



%% Section for setting up M5 and M6
%  set length



%  set the bias conditions



%  set the width



%  calculate remaining OP-parameters (fT, Cgs, etc.)



%  copy M5 to M6



%% Section for setting up  M1 and M2
%  set the length



%  set the bias conditions


%  set the width



%  calculate remaining OP-parameters (fT, Cgs, etc.)



%  copy M1 to M2



%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%% AI: Fill out the empty equations with the correct formula

AvDC =     % DC gain
C1  = ;     % Capacitance on node 1
G1   = ;                   % Resistance  on node 1
C2  = ;                  % Capacitance on node 2
G2   = ;         % Resistance  on node 2

Voutcm      = 
Vincm       = 
SpecGain    = 
pdom        = 
fGBW        =  % [Hz]    GBW frequency
SR          =  % [mV/ns] slew rate
Vcm_in_min  = 
Vcm_in_max  = 
Vout_cm_min = 
Vout_cm_max = 

%-------------------------------------------------------------------------

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

