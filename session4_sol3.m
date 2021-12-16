%% Analog Electronics Session 3: Noise 
% This session starts with a brief discussion about noise.

clc; 
close all; 
clear;

addpath(genpath('circuitDesign'));
addpath(genpath('functions'));
addpath(genpath('models'));
load('UMC65_RVT.mat');

%% Initializations
designkitName		= 'umc65';
circuitTitle		= 'Analog Design - Noise Example';
elementList.nmos	= {'Mn1'};
elementList.pmos	= {};
choice.maxFingerWidth = 10e-6;
choice.minFingerWidth = 200e-9;
simulator			= 'spectre';
simulFile			= 0;
simulSkelFile		= 0;
spec				= [];
noisedemo			= cirInit('noisedemo',circuitTitle,'top',elementList,spec,choice,...
						designkitName,NRVT,PRVT,simulator,simulFile,simulSkelFile);
noisedemo			= cirCheckInChoice(noisedemo, choice);


%% Frequencies for noise PSD integration
fnoise_min				= 1;	% [Hz]
fnoise_max				= 1e3;	% [Hz]
fnoise_spot1			= 1;	% [Hz]

fprintf('\n--- First exercise: Noise ---\n');
%% Circuit
disp('                                 ');
disp('          VDD                    ');
disp('           |                     ');
disp('       CurSour<--- I = cte       ');
disp('           |                     ');
disp('           +---- OUT             ');
disp('           |                     ');
disp('   VIN----Mn1                    ');
disp('           |                     ');
disp('          GND                    ');

%% Goal (ETA: 45min)
% The circuit is an NMOS amplifier with a perfect current source (it
% provides exactly the current that the amplifier needs). Reduce the noise
% contribution of the NMOS amplifier as much as possible

%% Design choices
%
% # Q: For the amplifier, do we want big or small gm for better thermal noise?
% A: As we want to decrease input referred voltage noise (dvn^2 =
% 4kT(2/3)/gm df), we want as high gm as possible
% # Q: If you want to decrease thermal noise as much as possible, to which
% inversion level do you go? What is the consequence for power consumption?
% A: Increasing gm (= mu*W/L*Cox*VOV) means we have to go to strong inversion
% # Q: design for constant gm, if you go towards weak inversion, does the
% power spectral density of the input referred thermal noise change
% significantly, does it increase, why? 
% A: no, the thermal noise voltage is 4kT/gm. For a constant gm, the noise
% will not change significantly.
% # Q: design for constant gm, if you go towards weak inversion, does the
% power spectral density of the input referred flicker noise change
% significantly, does it increase, why?
% A: Yes, going to weak inversion will make devices huge. We know that the
% flicker noise decreases with the area of a transistor.
% # Q: design for constant IDS, if you go towards weak inversion, does the
% power spectral density of the input referred thermal noise change
% significantly, does it increase, why? 
% % A: We expect the input referred noise to go down, because, for a
% constant IDS, the gm of a transistor increases when we go towards weak
% inversion. Therefore the input referred noise, which is 4kT/gm, should go
% down. 
% # Q: design for constant IDS, if you go towards weak inversion, does the
% power spectral density of the input referred flicker noise change
% significantly, does it increase, why? 
% A: Yes, going to weak inversion will make devices huge. We know that the
% flicker noise decreases with the area of a transistor.
%



KeepConstant = 'ids'; % {gm , ids}
Mn1.lg			= 5*tableLmin(NRVT);
Mn1.vds			= 0.6;


gm_transistor=ones(20,1);
width_transistor=ones(20,1);
noisethermal_transistor=ones(20,1);
noiseflicker_transistor=ones(20,1);
noise_transistor=ones(20,1);

for i=1:20
    
Mn1.vov		= i*0.01;
if strcmp(KeepConstant, 'gm')
    Mn1.gm      = 1e-3;
else
    Mn1.ids		= 200e-6;    
end

%% Implementation (don't touch)
Mn1.vsb			= 0;
Mn1.vth			= tableValueWref('vth',NRVT,Mn1.lg,0,Mn1.vds,Mn1.vsb);
Mn1.vgs			= Mn1.vov+Mn1.vth;
if strcmp(KeepConstant,'gm')
    Mn1.w		= mosWidth('gm',Mn1.gm,Mn1);
else
    Mn1.w		= mosWidth('ids',Mn1.ids,Mn1);
end
Mn1				= mosNfingers(Mn1);
Mn1				= mosOpValues(Mn1);
mosCheckSaturation(Mn1);

% Small-signal DC gain
Av__Mn1			= Mn1.gm/Mn1.gds;
Av__Mn1_dB		= 20*log10(Av__Mn1);

%% Noise calculations (don't touch!)
% Mn1.di2_fn is the noise spectral density of the flicker noise component
% without the 1/f dependency. So if we want to know what the noise spectral
% density is at frequency f1, we need to divide this Mn1.di2_fn by f1. To
% know the integrated noise power between frequency f1 and f2 we need to
% multiply Mn1.di2_fn with (log(f2)-log(f1)), where log is in base e.

% Input-referred integrated noise power
d_vn2_thermal_in__Mn1	= Mn1.di2_id/(Mn1.gm^2); 
% Noise power spectral density thermal Mn1
Pn_thermal_in__Mn1		= d_vn2_thermal_in__Mn1*...
    (fnoise_max-fnoise_min); % Integrated noise power thermal Mn1
d_vn2_flicker_in__Mn1	= Mn1.di2_fn/(Mn1.gm^2); 
% Noise power spectral density flicker Mn1
Pn_flicker_in__Mn1		= d_vn2_flicker_in__Mn1*...
    (log(fnoise_max)-log(fnoise_min)); % Integrated noise 
% power flicker Mn1
Pn_total_in__Mn1		= Pn_thermal_in__Mn1+Pn_flicker_in__Mn1; 
% Total input referred integrated noise Mn1

% Output-referred integrated noise power 
d_vn2_thermal_out__Mn1	= d_vn2_thermal_in__Mn1*(Av__Mn1^2);
Pn_thermal_out__Mn1		= d_vn2_thermal_out__Mn1*(fnoise_max-fnoise_min);
d_vn2_flicker_out__Mn1	= d_vn2_flicker_in__Mn1*(Av__Mn1^2);
Pn_flicker_out__Mn1		= d_vn2_flicker_out__Mn1*(log(fnoise_max)-log(fnoise_min));
Pn_total_out__Mn1		= Pn_thermal_out__Mn1+Pn_flicker_out__Mn1;

% Input referred noise spectral densities at frequency f1
d_vn2_thermal_in__Mn1_f1= d_vn2_thermal_in__Mn1;
d_vn2_flicker_in__Mn1_f1= d_vn2_flicker_in__Mn1/fnoise_spot1;
d_vn2_total_in__Mn1		= d_vn2_thermal_in__Mn1_f1+d_vn2_flicker_in__Mn1_f1;

% Output-referred noise spectral density at frequency f1
d_vn2_thermal_out__Mn1_f1   = d_vn2_thermal_in__Mn1_f1*(Av__Mn1^2);
d_vn2_flicker_out__Mn1_f1   = d_vn2_flicker_in__Mn1_f1*(Av__Mn1^2);
d_vn2_total_out__Mn1        = d_vn2_thermal_out__Mn1_f1+d_vn2_flicker_out__Mn1_f1;

%% print (don't touch)
% HW: Redo the exercise, but instead of an ideal current source, use a PMOS
% load instead. Try to size the PMOS for smallest input referred noise
% contribution.

noisedemo = cirElementsCheckOut(noisedemo);
mosPrintSizesAndOpInfo(1,noisedemo)

fprintf('\n\n');
fprintf('------------------------------------------------------------------------------\n');
fprintf(' Noise calculations:\n');
fprintf('------------------------------------------------------------------------------\n');
fprintf('\n');
fprintf('Noise density at f1= %gHz:\n',fnoise_spot1);
fprintf('\n');
fprintf('d_vn2_thermal_in__Mn1_f1    = %9.4f aV^2/Hz (%9.4f nVrms/sqrt(Hz))\n',d_vn2_thermal_in__Mn1_f1*1e18,sqrt(d_vn2_thermal_in__Mn1_f1)*1e9);
fprintf('d_vn2_flicker_in__Mn1_f1    = %9.4f pV^2/Hz (%9.4f uVrms/sqrt(Hz))\n',d_vn2_flicker_in__Mn1_f1*1e12,sqrt(d_vn2_flicker_in__Mn1_f1)*1e6);
fprintf('d_vn2_total_in__Mn1_f1      = %9.4f pV^2/Hz (%9.4f uVrms/sqrt(Hz))\n',d_vn2_total_in__Mn1*1e12,sqrt(d_vn2_total_in__Mn1)*1e6);
fprintf('\n');
fprintf('d_vn2_thermal_out__Mn1_f1   = %9.4f aV^2/Hz (%9.4f nVrms/sqrt(Hz))\n',d_vn2_thermal_out__Mn1_f1*1e18,sqrt(d_vn2_thermal_out__Mn1_f1)*1e9);
fprintf('d_vn2_flicker_out__Mn1_f1   = %9.4f pV^2/Hz (%9.4f uVrms/sqrt(Hz))\n',d_vn2_flicker_out__Mn1_f1*1e12,sqrt(d_vn2_flicker_out__Mn1_f1)*1e6);
fprintf('d_vn2_total_out__Mn1_f1     = %9.4f pV^2/Hz (%9.4f uVrms/sqrt(Hz))\n',d_vn2_total_out__Mn1*1e12,sqrt(d_vn2_total_out__Mn1)*1e6);
fprintf('\n');
fprintf('Integrated noise in the [%4.1e, %4.1e] Hz range:\n',fnoise_min,fnoise_max);
fprintf('\n');
fprintf('Pn_thermal_in__Mn1     = %9.4f pV^2 (%9.4f nVrms)\n',Pn_thermal_in__Mn1*1e12,sqrt(Pn_thermal_in__Mn1)*1e9);
fprintf('Pn_flicker_in__Mn1     = %9.4f pV^2 (%9.4f nVrms)\n',Pn_flicker_in__Mn1*1e12,sqrt(Pn_flicker_in__Mn1)*1e9);
fprintf('Pn_total_in__Mn1       = %9.4f pV^2 (%9.4f nVrms)\n',Pn_total_in__Mn1*1e12,sqrt(Pn_total_in__Mn1)*1e9);
fprintf('\n');
fprintf('Pn_thermal_out__Mn1    = %9.4f pV^2 (%9.4f nVrms)\n',Pn_thermal_out__Mn1*1e12,sqrt(Pn_thermal_out__Mn1)*1e9);
fprintf('Pn_flicker_out__Mn1    = %9.4f pV^2 (%9.4f nVrms)\n',Pn_flicker_out__Mn1*1e12,sqrt(Pn_flicker_out__Mn1)*1e9);
fprintf('Pn_total_out__Mn1      = %9.4f pV^2 (%9.4f nVrms)\n',Pn_total_out__Mn1*1e12,sqrt(Pn_total_out__Mn1)*1e9);
fprintf('\n');
fprintf('\n');

gm_transistor(i)=Mn1.gm;
width_transistor(i)=Mn1.w;
noisethermal_transistor(i)=d_vn2_thermal_in__Mn1_f1;
noiseflicker_transistor(i)=d_vn2_flicker_in__Mn1_f1;

end


subplot(2,2,1)
plot(gm_transistor,width_transistor,'LineWidth',2);
title('Transistor width vs gm')
ylabel('Transistor width (m)') 
xlabel('gm (S)')
subplot(2,2,2)
plot(gm_transistor,noisethermal_transistor,'LineWidth',2);
title('Input referred thermal noise vs gm')
ylabel('Input thermal noise V^2/Hz') 
xlabel('gm (S)')
subplot(2,2,3)
plot(gm_transistor,noiseflicker_transistor,'LineWidth',2);
title('Input referred flicker noise vs gm')
ylabel('Input flicker noise V^2/Hz') 
xlabel('gm (S)')







%% Project (ETA: over 9000)
% The goal is to construct a two stage miller compensated opamp. You will
% receive a document with further information.
% We recommend you to start the implementation in a clean script. You can copy
% the code from other exercises. 
% Don't forget to correctly define all used transistors in the beginning of
% the script. Also don't forget that transistors with a nonzero vsb undergo
% the body effect (in order to calculate the threshold voltage vth (use
% tableValueWref), you first need to use function mosVsbBody to derive the
% vsb voltage). You can choose to either start from the
% second stage or from the first. Both approaches give the same results.
% For basic understanding of this project, you can use the AEslides.pdf to
% guide you. 
% As always... good luck!
% Q: In both stages, do you first size the amplifier or the load transistor?

%% What we have learnt this session
%
% # How to reduce 1/f noise
% # How to reduce thermal noise
%