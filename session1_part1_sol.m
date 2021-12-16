%% Analog Electronics Session 1 (Intro): Transistor curves + Matlab environment
% Welcome to the lab sessions of the course of Analog Electronics. This 
% matlab file is the starting point for the labs. 
% Good luck!!

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
circuitTitle    = 'Analog Design - Session 1 - IV curve';

%Declaration of the circuit components
elementList.nmos = {'Mn1'};
elementList.pmos = {};

spec.VDD        = 1.1;
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator       ='spectre';
simulFile       = 0;
simulSkelFile   = 0;
analog = cirInit('analog', circuitTitle, 'top', elementList, spec , choice,...
    designkitName, NRVT, PRVT, simulator, simulFile, simulSkelFile);

analog          = cirCheckInChoice(analog, choice);

%% Circuit
disp('        Vd            ');
disp('        |             ');
disp('  Vg---Mnx            ');
disp('        |             ');
disp('        Vs            ');

fprintf('\n--- First Exercise: Designing transistor from scratch ---\n');

%% Implementation
VDS = 0.01:0.05:spec.VDD;
VGS = 0.01:0.05:spec.VDD;

IDS = NaN(length(VGS),length(VDS));
sat = NaN(length(VGS),length(VDS));

Mn1.lg = 65e-9;
Mn1.w = 10*Mn1.lg;
Mn1.nFingers = 1;


%% for matlab
 if ~isfile('IDS.mat')
%% for octave
% if exist('IDS.mat') ~= 2
for jj=1:length(VDS)
    for ii=1:length(VGS)
        Mn1.vgs = VGS(ii);
        Mn1.vds = VDS(jj);
        Mn1.vsb = 0;
        Mn1.vth = tableValueWref('vth', NRVT, Mn1.lg, Mn1.vgs, Mn1.vds, Mn1.vsb);
        Mn1 = mosOpValues(Mn1);
        IDS(ii,jj) = Mn1.ids;
        sat(ii,jj) = mosCheckSaturation(Mn1);

    end
end
[X,Y] = meshgrid(VDS,VGS);
 else
     load('IDS.mat');
 end

%% Figures of Merit + Plot
save('IDS.mat','X','Y','IDS','sat');
figure; 
ax1 = subplot(311); surf(X,Y,IDS); title('IDS in function of VGS and VDS (mag)');
ylabel('VGS [V]'); xlabel('VDS [V]'); zlabel('IDS [A]');

ax2 = subplot(312); surf(X,Y,20*log10(IDS)); title('IDS in function of VGS and VDS (log)');
ylabel('VGS [V]'); xlabel('VDS [V]'); zlabel('IDS [A]');

ax3 = subplot(313); surf(X,Y,sat); title('MOS in saturation as function of VGS and VDS');
ylabel('VGS [V]'); xlabel('VDS [V]'); zlabel('Saturation [0,1]');

hlink = linkprop([ax1,ax2],{'CameraPosition'});
rotate3d on;
