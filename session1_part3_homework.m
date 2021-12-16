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
elementList.nmos = {'Mn3','Mn4'};
elementList.pmos = {'Mp3','Mp4'};

spec.VDD        = 1.1;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator       ='spectre';
simulFile       = 0;
simulSkelFile   = 0;
analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice,...
    designkitName, NRVT, PRVT, simulator, simulFile, simulSkelFile);

analog          = cirCheckInChoice(analog, choice);

fprintf('\n--- Homework 1: Common source amplifier with cascodes ---\n');
%% HW1: Circuit
disp('       VDD           ');
disp('        |            ');
disp('   IN--Mp3           ');
disp('        |            ');
disp('       Mp4           ');
disp('        |            ');
disp('        +----+-OUT   ');
disp('        |    |       ');
disp('       Mn4   |       ');
disp('        |    CL      ');
disp('       Mn3   |       ');
disp('        |    |       ');
disp('       GND  GND      ');

%% HW1: Specs
spec.fGBW       = 80e6;    % [Hz] GBW frequency
spec.Cl         = 15e-12;   % [F] load capacitance
spec.VDD        = 1.1;      % [V] Power supply voltage
spec.Vswing     = 0.2;      % [V] minimum peak to peak voltage swing of output
bodyeffecton    = 1;

%% HW1: Goal (ETA: 45min)
% Get the highest gain possible with this amplifier, you can either
% assume triple well technology (Mn2.vsb = 0). The daredevils can solve
% this exercise with body effect taken into account. (widths < 100um)
% TIP: 46dB is possible

%% HW1: Design Choices
% When the body effect is not playing, design is easy. We start with Mn3.
% We calculate the sizes in order to get the wanted gm. All the other
% transistors are sized in order to drive the current that Mn3 is asking
% for.
% We must play with VG of Mn4 in order to have a VSB for Mn4 which is the
% same as VDS of Mn3 that we wanted

Mp3.lg  = 5*tableLmin(PRVT);    % [m] Design choice
Mp3.vov = -0.1;                 % [V] Design choice: Choice of inversion
Mp3.VD  = spec.VDD-(spec.VDD-spec.Vswing)/4; % [V]

Mp4.lg  = 5*tableLmin(PRVT);    % [m] Design choice
Mp4.vov = -0.05;                % [V] Design choice: Choice of inversion
% low vov, such that for same current, Mp4 has bigger gm such that source
% degenerated impedance is bigger.

Mn3.lg  = 4*tableLmin(NRVT);    % [m] Design choice
Mn3.vov = 0.1;                 % [V] Design choice: Choice of inversion

if bodyeffecton
    Mn4.VG  = 0.6; % With bodyeffect we have to choose VG because VTH is 
    % depending on VGB. VS will be chosen accordingly
else
    Mn4.VS  = (spec.VDD-spec.Vswing)/4; % Without bodyeffect we have can 
    % either choose VG, VS will be depending on VOV and VTH. Either we
    % choose VS, together with VTH and VOV we know VG.
end
Mn4.lg  = 4*tableLmin(NRVT);    % [m] Design choice
Mn4.vov = 0.15;                % [V] Design choice: Choice of inversion
Mn4.VD  = spec.VDD/2;

%% HW1: Implementation, Don't touch!
% Mp3
Mp3.VS  = spec.VDD;
Mp3.vsb = 0;                % Bulk connected to lowest potential
Mp3.vds = Mp3.VD-Mp3.VS;    % Maximum swing
Mp3.vth = tableValueWref('vth',PRVT,Mp3.lg,0,Mp3.vds,Mp3.vsb);
Mp3.vgs = Mp3.vov + Mp3.vth;% By definition
Mp3.gm  = 2*pi*spec.fGBW*spec.Cl;   % gm is fixed by GBW
Mp3.w   = mosWidth('gm',Mp3.gm,Mp3);% width for given gm
Mp3     = mosNfingers(Mp3); % Design choice: cut total width in pieces
Mp3     = mosOpValues(Mp3); % Calculating operating point values
% such as DC current, small signal parameters such as gds, parasitic
% capacitances, noise figure

Mp4.VD  = spec.VDD/2;
Mp4.VS  = Mp3.VD;
Mp4.vsb = 0;            % Bulk connected to lowest potential
Mp4.vds = Mp4.VD - Mp4.VS;
Mp4.vth = tableValueWref('vth',PRVT,Mp4.lg,0,Mp4.vds,Mp4.vsb);
Mp4.vgs = Mp4.vov + Mp4.vth;% By definition
Mp4.VG  = Mp4.vgs+Mp3.vds;
Mp4.ids = Mp3.ids;
Mp4.w   = mosWidth('ids',Mp4.ids,Mp4);% width for given ids
Mp4     = mosNfingers(Mp4); % Design choice: cut total width in pieces
Mp4     = mosOpValues(Mp4); % Calculating operating point values
% such as DC current, small signal parameters such as gds, parasitic
% capacitances, noise figure

% Mn4
Mn4.VB  = 0;
if bodyeffecton
    Mn4.vsb = mosVsbBody(Mn4, Mn4.VG, Mn4.VD, Mn4.vov, 0.2);
    Mn4.VS  = Mn4.vsb+Mn4.VB;
else
    Mn4.vsb = 0;                % Bulk connected to lowest potential
end
Mn4.vds = Mn4.VD-Mn4.VS;    % Maximum swing
Mn4.vth = tableValueWref('vth',NRVT,Mn4.lg,0,Mn4.vds,Mn4.vsb);
Mn4.vgs = Mn4.vov + Mn4.vth;% By definition
Mn4.gm  = 2*pi*spec.fGBW*spec.Cl;   % gm is fixed by GBW
Mn4.w   = mosWidth('ids',Mn4.gm,Mn4);% width for given gm
Mn4     = mosNfingers(Mn4); % Design choice: cut total width in pieces
Mn4     = mosOpValues(Mn4); % Calculating operating point values
% such as DC current, small signal parameters such as gds, parasitic
% capacitances, noise figure

% Mn3
Mn3.VD  = Mn4.VS;
Mn3.VS  = 0;
Mn3.vsb = 0;                % Bulk connected to lowest potential
Mn3.vds = Mn3.VD-Mn3.VS;    % Maximum swing
Mn3.vth = tableValueWref('vth',NRVT,Mn3.lg,0,Mn3.vds,Mn3.vsb);
Mn3.vgs = Mn3.vov + Mn3.vth;% By definition
Mn3.VG  = Mn3.vgs + Mn3.VS;
Mn3.ids = Mp3.ids;
Mn3.w   = mosWidth('ids',Mn3.ids,Mn3);% width for given ids
Mn3     = mosNfingers(Mn3); % Design choice: cut total width in pieces
Mn3     = mosOpValues(Mn3); % Calculating operating point values
% such as DC current, small signal parameters such as gds, parasitic
% capacitances, noise figure

%% HW1: Figures of Merit + plot
%
% # Q: Why can stacking of transistors be interesting? (4.8 Common-gate stage,
% 5.2 Cascode stage)
% # Q: What is the difference if the input would be on Mn2 instead of Mn1
% (4.7 Source degeneration)
% # Q: How many low frequent poles do we have with a one stage amplifier
% with cascodes? How many when using a two amplifier stage wihtout?
% # Q: Is the stability better for a two stage amplifier or a one stage
% amplifier with stacking? (no miller compensation)
% # Q: Given a VOV of 0.2V for all transistors, what is approximately the
% lowest supply voltage we can have?
% # Q: What is the formula for the output impedance? A: rout =
% Nro4*Nro3*Ngm4 // Pro4*Pro3*Pgm4
% # Q: What is the formula for the input impedance? A: 1/(s*CGS)
% # Q: What is the formula for dominant pole? A: 1/(2pi*rout*CL)
% # Q: What is the formula for voltage gain? A: Ngm3*rout
% # Q: Give the formula for GBW?  Ngm3/(2pi*CL)
% # Q: If you want to increase the output impedance seen from the PMOS
% branch, do you increase or decrease vov of Mp4? Why? Do you expect to
% change Mn4.gds in first order?
%

routp = 1/Mp4.gds*(1+Mp4.gm/Mp3.gds);
routn = 1/Mn3.gds*(1+Mn3.gm/Mn3.gds);
rout = routp*routn/(routp+routn);
pole   = -1/(rout*(spec.Cl + Mp4.cdb + Mp4.cgd + ...
    + Mp3.cdb + Mp3.cgd));
Gain = Mp3.gm * rout;
Gain_dB = 20*log10(Gain);
fGBW = Gain * pole / (2*pi);
Vswing = spec.VDD - abs(Mp3.vds) - abs(Mp3.vdsat) ...
    - Mp3.vds - Mp4.vdsat;

fprintf('\nSolution of exercise 3:\n');
fprintf('rout_n = %6.2f KOhm\n', routn/1000);
fprintf('rout_p = %6.2f KOhm\n', routp/1000);
fprintf('Gain = %.2f\n',Gain);
fprintf('Gain_dB = %.2f dB\n',Gain_dB);
fprintf('width of Mp3   = %.2f um \n', Mp3.w*1e6);
fprintf('width of Mp4   = %.2f um \n', Mp4.w*1e6);
fprintf('width of Mn4   = %.2f um \n', Mn4.w*1e6);
fprintf('width of Mn3   = %.2f um \n', Mn3.w*1e6);

% fprintf('Mn3.vds = %gV, Mn4.vsb = %gV\n',Mn3.vds,Mn4.vsb);

freq    = logspace(1,10,1e3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if control toolbox in Matlab is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##s       = tf('s');
##TF3     = Gain/(1-s/pole);
##figure;
##bode(TF3, freq); grid on;
##title('Voltage gain (Cascode stage)');
%%%%%%%%
%% else
%%%%%%%%
TF3 = Gain./(1-1j*2*pi*freq/pole);
figure;
subplot(211);
semilogx(freq, 20*log10(abs(TF3)));
ylabel('Magnitude [dB]');
xlabel('Frequency [Hz]');
title('Bode plot');
subplot(212);
semilogx(freq, (angle(TF3))*180/pi);
ylabel('Phase [degree]');
xlabel('Frequency [Hz]');

