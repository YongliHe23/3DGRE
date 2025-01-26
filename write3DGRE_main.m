%3D GRE sequence with long TR
% with a 3D spin-warp readout
% modified based on write3DGRE.m from Jon's Pulseq mannual v1.0
% by Yongli He on 2025/01/25

% System/design parameters.
sys = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 200, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 0, ...
              'adcRasterTime', 2e-6, ...
              'gradRasterTime', 4e-6, ...
              'blockDurationRaster', 4e-6, ...
              'rfRasterTime',4e-6,...
              'B0', 3.0);

% Acquisition parameters
fov = [240e-3 240e-3 240e-3];   % FOV (m)
%fov=[50e-3 50e-3 50e-3];
Nx = 60; Ny = Nx; Nz = 60;    % Matrix size
TR = 55e-3;                     % sec
%TR=200e-3;
dwell = 8e-6;                  % ADC sample time (s)
alpha = 20;                      % flip angle (degrees)
alphaPulseDuration = 0.2e-3;
nCyclesSpoil = 2;               % number of spoiler cycles
Tpre = 1.0e-3;                  % prephasing trapezoid duration
rfSpoilingInc = 117;            % RF spoiling increment

% Create a new sequence object
seq = mr.Sequence(sys);           

%for beta pulse
sys_beta = mr.opts('maxGrad', 40, 'gradUnit','mT/m', ...
              'maxSlew', 200, 'slewUnit', 'T/m/s', ...
              'rfDeadTime',sys.rfDeadTime, ... %100e-6, ...
              'rfRingdownTime',sys.rfRingdownTime, ... % 60e-6, ...
              'adcDeadTime',sys.adcDeadTime, ... %0, ...
              'adcRasterTime', 2e-6, ...
              'gradRasterTime', 4e-6, ...
              'blockDurationRaster', 4e-6, ...
              'rfRasterTime',4e-6,...
              'B0', 3.0);

%% create saturation pulse
load ../../ivext_ellipse_2D_wholeFOV_alpha=20_seq.mat;


rf_sat=mr.makeArbitraryRf(Rf_sat*4258,3.14,'system',sys_beta,'delay',100e-6);
rf_sat.signal= rf_sat.signal/max(abs(rf_sat.signal))*max(abs(Rf_sat*4258)); % ensure correct amplitude (Hz)

%rotate gradient along x-axis by 90 deg <-- this turns out to be unnecessary
% now rotate along z-axis by 45 deg
gx_sat=mr.makeArbitraryGrad('x',(Gx_sat-Gy_sat)*(1/sqrt(2))*425.8e3,'system',sys_beta,'delay',100e-6);
%gx_sat=mr.makeArbitraryGrad('x',Gx_sat*425.8e3,'system',sys_beta,'delay',100e-6);
gy_sat=mr.makeArbitraryGrad('y',(Gx_sat+Gy_sat)*(1/sqrt(2))*425.8e3,'system',sys_beta,'delay',100e-6);
%gy_sat=mr.makeArbitraryGrad('y',Gy_sat*425.8e3,'system',sys_beta,'delay',100e-6);
gz_sat=mr.makeArbitraryGrad('z',Gz_sat*425.8e3,'system',sys_beta,'delay',100e-6);

%% Create slice-selective pulse
%[rf] = mr.makeBlockPulse(alpha/180*pi, sys, 'Duration', alphaPulseDuration);
[rf,gzRF, gzRF_r] = mr.makeSincPulse(alpha/180*pi, sys, 'Duration', 4e-3, 'SliceThickness', 0.35*fov(3), 'apodization',0.5,'timeBwProduct',8);
%sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
%     'maxSlew', 0.7*sys.maxSlew/sys.gamma/10, ...           % G/cm/ms
%     'maxRF', 0.15);
% mb=1;
% sliceSep = fov(3)/mb;   % center-to-center separation between SMS slices (m)
% % freq = frequency offset (Hz) corresponding to a sliceSep
% ftype = 'min';
% slThick=0.35*fov(3);
% rfTB  = 8;          % RF pulse time-bandwidth product
% rfDur = 4e-3;       % RF pulse duration (s)
% 
% addpath /home/yonglihe/Documents/MATLAB/HarmonizedMRI/SMS-EPI/sequence/Pulseq
% [rf, gzRF, freq] = getsmspulse(alpha, slThick, rfTB, rfDur, ...
%     mb, sliceSep, sysGE, sys, ...
%     'doSim', true, ...    % Plot simulated SMS slice profile
%     'type', 'st', ...     % SLR choice. 'ex' = 90 excitation; 'st' = small-tip
%     'ftype', ftype);      % filter design. 'ls' = least squares

%% Create fat sat pulse 
arg.fatsat=1;
arg.fatFreqSign=-1;

fatsat.flip    = 90;      % degrees
fatsat.slThick = 1e5;     % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
fatsat.tbw     = 3.5;     % time-bandwidth product
fatsat.dur     = 8.0;     % pulse duration (ms)

% RF waveform in Gauss
wav = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, fatsat.dur, 1e-6, toppe.systemspecs(), ...
    'type', 'ex', ...    % fatsat pulse is a 90 so is of type 'ex', not 'st' (small-tip)
    'ftype', 'min', ...
    'writeModFile', false);

% Convert from Gauss to Hz, and interpolate to sys.rfRasterTime
rfp = rf2pulseq(wav, 4e-6, sys.rfRasterTime);

% Create pulseq object
% Try to account for the fact that makeArbitraryRf scales the pulse as follows:
% signal = signal./abs(sum(signal.*opt.dwell))*flip/(2*pi);
flip = fatsat.flip/180*pi;
flipAssumed = abs(sum(rfp));
rfsat = mr.makeArbitraryRf(rfp, ...
    flip*abs(sum(rfp*sys.rfRasterTime))*(2*pi), ...
    'system', sys);
rfsat.signal = rfsat.signal/max(abs(rfsat.signal))*max(abs(rfp)); % ensure correct amplitude (Hz)
rfsat.freqOffset = arg.fatFreqSign*520;  % Hz

%%
% Define other gradients and ADC events
% Cut the redaout gradient into two parts for optimal spoiler timing
deltak = 1./fov;
Tread = Nx*dwell;

gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', Ny*deltak(2)/2, ...   % PE1 gradient, max positive amplitude
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)/2, ...   % PE2 gradient, max positive amplitude
    'Duration', Tpre);

gxtmp = mr.makeTrapezoid('x', sys, ...  % readout trapezoid, temporary object
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread);
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gxtmp.area/2, ...
    'Duration', Tpre);

adc = mr.makeAdc(Nx, sys, ...
    'Duration', Tread,...
    'Delay', gxtmp.riseTime);

adc.deadTime=0;
% extend flat time so we can split at end of ADC dead time
gxtmp2 = mr.makeTrapezoid('x', sys, ...  % temporary object
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread + 20e-6);   
[gx, ~] = mr.splitGradientAt(gxtmp2, gxtmp2.riseTime + gxtmp2.flatTime,sys);

gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
gxSpoil = mr.makeExtendedTrapezoidArea('x', gxtmp.amplitude, 0, gzSpoil.area, sys);

gxSpoil_sat=mr.makeTrapezoid('x', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);

gzSpoil_sat=mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);

% y/z PE steps
pe1Steps = ((0:Ny-1)-Ny/2)/Ny*2;
pe2Steps = ((0:Nz-1)-Nz/2)/Nz*2;

% Calculate TR delay
TRmin = mr.calcDuration(rf)+mr.calcDuration(gzRF_r) + mr.calcDuration(gxPre) ...
   + mr.calcDuration(gx) + mr.calcDuration(gxSpoil)+max(mr.calcDuration(rf_sat),mr.calcDuration(gx_sat))...
   +mr.calcDuration(gxSpoil_sat)+mr.calcDuration(gzSpoil)+(mr.calcDuration(rfsat)+mr.calcDuration(gxSpoil_sat))*arg.fatsat;
delayTR = TR - TRmin;

% Loop over phase encodes and define sequence blocks
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC is turned on and used for receive gain calibration on GE scanners
% iZ > 0: Image acquisition

nDummyZLoops = 1;

rf_phase = 0;
rf_inc = 0;

%for iZ = -nDummyZLoops:2
for iZ=-nDummyZLoops:Nz
    isDummyTR = iZ < 0;

    msg = sprintf('z encode %d of %d   ', iZ, Nz);
    for ibt = 1:length(sprintf('z encode %d of %d   ', iZ - 1, Nz))
        fprintf('\b');
    end
    fprintf(msg);
    
    for iY = 1:Ny
        % Turn on y and z prephasing lobes, except during dummy scans and
        % receive gain calibration (auto prescan)
        yStep = (iZ > 0) * pe1Steps(iY);
        zStep = (iZ > 0) * pe2Steps(max(1,iZ));

        if arg.fatsat
            % RF spoiling for last excitation
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);

            % fat sat
            rfsat.phaseOffset = rf_phase/180*pi;
            seq.addBlock(rfsat, mr.makeLabel('SET', 'TRID', 2-isDummyTR));
            seq.addBlock(gxSpoil_sat, gzSpoil_sat);
            
            %RF spoiling for fatsat
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);

            %beta pulse
            rf_sat.phaseOffset=rf_phase/180*pi;
            seq.addBlock(rf_sat,gx_sat,gy_sat,gz_sat);
            seq.addBlock(gxSpoil_sat,gzSpoil_sat);

            % excitation pulse and RF spoiling
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            seq.addBlock(rf, gzRF);
            seq.addBlock(gzRF_r);
            %seq.addBlock(rf);

        else
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
    
            %beta pulse
            rf_sat.phaseOffset = rf_phase/180*pi;
            seq.addBlock(rf_sat,gx_sat,gy_sat,gz_sat,mr.makeLabel('SET', 'TRID', 2-isDummyTR));
            %beta pulse gradient crusher
            seq.addBlock(gxSpoil_sat,gzSpoil_sat)
    
            % RF spoiling for beta pulse
            rf.phaseOffset = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
            
            % Excitation
            % Mark start of segment (block group) by adding label.
            % Subsequent blocks in block group are NOT labelled.
            seq.addBlock(rf,gzRF);
            seq.addBlock(gzRF_r);
            %seq.addBlock(rf);           
            
        end
            % Encoding
            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre, yStep), ...
                mr.scaleGrad(gzPre, zStep));
            if isDummyTR
                seq.addBlock(gx);
            else
                seq.addBlock(gx, adc);
            end
    
            % rephasing/spoiling and TR delay
            seq.addBlock(gxSpoil, ...
                mr.scaleGrad(gyPre, -yStep), ...
                mr.scaleGrad(gzPre, -zStep));
            seq.addBlock(gzSpoil)
            
            seq.addBlock(mr.makeDelay(delayTR));    
    end
end
fprintf('Sequence ready\n');

% Check sequence timing
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Output for execution
ifn='ov-spr_sinc-alpha.seq';
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'gre');
seq.write(ifn);

%% Optional plots

% Plot sequence
Noffset = Ny*(nDummyZLoops+1);
seq.plot('timerange',[2*Noffset 2*Noffset+4]*TR, 'timedisp', 'ms');

%% covert to .tar file
sysGE = toppe.systemspecs('maxGrad', 5, ...   % G/cm
    'maxRF', 0.15, ...
    'maxSlew', 20, ...                        % G/cm/ms
    'adcDeadTime', 20, ...           % us. Half of 40us since applied both before + after ADC window.
    'psd_rf_wait', 148, ...          % RF/gradient delay (us)
    'psd_grd_wait', 156);            % ADC/gradient delay (us)

ceq = seq2ceq(ifn);
ofn = ['ov-spr_sinc-alpha' '.tar'];
ceq2ge(ceq, sysGE, ofn, 'preserveArea', false);
system(sprintf('tar xf %s', ofn));
return

% Plot k-space (2d)
[ktraj_adc,t_adc,ktraj,t_ktraj,t_excitation,t_refocusing] = seq.calculateKspacePP();
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');
