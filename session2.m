%% Analog Electronics Session 2: cascade amplifier with Miller compensation

%% Adding paths + Loading MOS tables
addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));

% clear;
% close all;
% clc;

load ('UMC65_RVT.mat');

%% Initialize everything
designkitName   = 'umc65';
circuitTitle    = 'Analog Design - Session 3';

%Declaration of the circuit components
elementList.nmos = {'Mn1','Mn2'};
elementList.pmos = {'Mp2'};

spec.VDD            = 1.1;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator       ='spectre';
simulFile       = 0;
simulSkelFile   = 0;

analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice,...
    designkitName, NRVT, PRVT, simulator, simulFile, simulSkelFile);

analog          = cirCheckInChoice(analog, choice);

fprintf('\n--- First Exercise: Miller cap ---\n');
%% EX: Circuit
disp('                                      ');
disp('       VDD            VDD                 ');
disp('        |              |                  ');
disp('        R1            Mp2                 ');
disp('        |              |                  ');
disp('        |-----+-Cm--Rm-+----+-----+-OUT2  ');
disp('        |     |        |    |     |       ');
disp(' IN1---Mn1    +-------Mn2   Cl    RL      ');
disp('        |              |    |     |       ');
disp('        |              |    |     |       ');
disp('       GND            GND  GND   GND      ');

%% EX: Specs
spec.fGBW   = 150e6;       % [Hz] GBW frequency
spec.gain   = 80;          % [] voltage gain
spec.CL     = 1e-12;       % [F], Load cap
spec.RL     = 5e3;         % [ohm] load resistor

%% Miller Capacitance and resistance
spec.Cm     = 0.0;    % [F], miller cap
spec.Rm     = 0.0;         % [F], Load cap

%% EX: Second stage NMOS: Design choices
VOUT        = ;   % [V], DC output voltage
Mn2.vov     = ;          % [V], overdrive voltage
Mn2.lg      = ;       % [m], channel length

%% EX: Second stage NMOS: Implementation


Mn2     = mosNfingers(Mn2);
Mn2     = mosOpValues(Mn2);

%% EX: Second stage PMOS: Design choices
Mp2.vov     = ;         % [V], overdrive voltage
Mp2.lg      = ;        % [m], channel length

%% EX: Second stage PMOS: Implementation;

Mp2     = mosNfingers(Mp2);
Mp2     = mosOpValues(Mp2);

%% EX: First stage NMOS: Design choices
Mn1.vov     = ;         % [V], overdrive voltage
Mn1.lg      = ;        % [m], channel length

%% EX: First stage NMOS: Implementation

Mn1     = mosNfingers(Mn1);
Mn1     = mosOpValues(Mn1);

%% First Stage Resistance: given stage 1 parameters


%% EX: Figures of Merit + plot
Av1 = ;
Av2 = ;

p1 = /...
    ;

p2 = / ...
    ;

z1    = ;

% Miller C calculated to achieve specification in GBW
Cccalc = ;
%R calculated for removing non-dominant pole
Rmcalc = ;

gainn = Av1*Av2;
if abs(p1)<abs(p2)
    fGBWn = gainn*abs(p1)/2/pi;
else 
    fGBWn = gainn*abs(p2)/2/pi;
end

fprintf('\n=== Results EX2 ===\n');
fprintf('\nVINT = %gV\n',Mn2.vgs);
fprintf('\nFirst stage: Gain = %g\n',Av1);
fprintf('Second stage: Gain = %g\n1/gdsMn2 = %gOhm, 1/gdsMp2 = %gOhm\n',Av2,1/Mn2.gds,1/Mp2.gds);
fprintf('Total current consumption: %6.2fmA\n',(Mn1.ids+Mn2.ids)/1e-3);

fprintf('\n\t Spec \t\t Actual\n')
fprintf('fGBW: \t %d MHz \t %d MHz\n',...
    spec.fGBW/1e6,round(fGBWn/1e6));
fprintf('gain: \t %g \t\t %g\n',...
    spec.gain,gainn);

freq    = logspace(1,12,1e3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if control toolbox in Matlab is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s   = tf('s');
TF1 = gainn*(1-s*(1/z1))/((1-s/p1)*(1-s/p2));
figure(1);
bode(TF1,freq); grid on;
h = gcr;
setoptions(h,'FreqUnits','Hz');
title('Frequency response cascaded amplifier');
hold all;
figure(2)
TFfb=feedback(TF1,1);
step(TFfb);
title('Step response cascaded amplifier');
hold all

fprintf('\nTransistors in saturation:\n');
if mosCheckSaturation(Mn1)
    fprintf('Mn1:Success\n')
end
if mosCheckSaturation(Mn2)
    fprintf('Mn2:Success\n')
end
if mosCheckSaturation(Mp2)
    fprintf('Mp2:Success\n')
end

%% Print sizes
analog = cirElementsCheckOut(analog); % Update circuit file with 
% transistor sizes
mosPrintSizesAndOpInfo(1,analog); % Print the sizes of the 
% transistors in the circuit file


