clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% setup
% All units are Hz and Volts

% test signal setup
fsymb       = 7.5e6;      % test symbols frequency
analog_ovrs = 25;         % oversampling to represent analog signals
cutoff      = 3e6;        % signal filter cutoff frequency (1st order, -3dB)
% ADC setup
fs          = 472.5e6;    % ADC raw sample rate
sd_ovrs     = fs/fsymb;   % ADC oversampling rate
diff_bias1  = +0.01;
diff_bias2  = -0.01;
cmp_bias    = -0.2;
Vref        =  2.5; 
noise_sigma =  0.02; 
% simulation time
sim_time    = 0.1e-3;     % simulation time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate input signal

% time sequence
t = (0:1/(fs*analog_ovrs):sim_time-1/(fs*analog_ovrs));

% input signal - random PAM3
in4 = randi([-1 +1], 1, length(t)/(analog_ovrs*(fs/fsymb)));
in4 = reshape(repmat(in4, analog_ovrs*(fs/fsymb), 1), 1, length(t));

% signal filter
rec_tf   = c2d(tf(1, [1/(2*pi*cutoff) 1]), 1/(fs*analog_ovrs));   
[flt.tx_filter_num, flt.tx_filter_den] = tfdata(rec_tf, 'v');

% adc input signal
adc_input = filter(flt.tx_filter_num, flt.tx_filter_den, in4);

% resample input signal for plot only
adc_input_r = resample(adc_input, 1, analog_ovrs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADC analog part
dac_out     = zeros(1, length(adc_input)/analog_ovrs);
noise       = noise_sigma * randn(3, length(t));
dac         = 0;
integrator1 = 0;
integrator2 = 0;


% filter coefficients
a1 = 2;
b1 = 0.5;
a2 = 1;
b2 = 50;

for i=1:length(adc_input)
    err1        = adc_input(i)- b1*dac + diff_bias1;    % error calculation 
    integrator1 = integrator1 + a1*err1 + noise(1, i);  % integrator 
    
    err2        = integrator1 - b2*dac + diff_bias2;    % error calculation 
    integrator2 = integrator2 + a2*err2 + noise(2, i);  % integrator
    
    if(mod(i, analog_ovrs)==0)               % comparator and sampler
        dac  = Vref*sign(integrator2 + cmp_bias + noise(3, i));
        dac_out(i/analog_ovrs) = dac;    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADC digital part

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% giant filter just for testing
% f_pass=3
% f_block=13
% a_block=-150dB

flt2_num = [  -0.0001573141843083,-8.633271245897e-05,-9.917286695792e-05, -0.0001030285907574, ...
  -9.276608939232e-05, -6.227238464387e-05, -4.259718936231e-06, 8.946081417693e-05, ...
  0.0002281307394882,0.0004216319357445,0.0006805245562375, 0.001015528639767, ...
   0.001437574466143, 0.001957222069436, 0.002584447537782, 0.003328222462836, ...
   0.004196053486178, 0.005193674775058, 0.006324679572635,  0.00759012972115, ...
   0.008988407884082,  0.01051484141936,  0.01216156637634,  0.01391741208812, ...
    0.01576814883838,  0.01769631926666,  0.01968158726936,  0.02170084686872, ...
    0.02372889933714,  0.02573881089455,  0.02770239264592,  0.02959065647781, ...
    0.03137483602187,  0.03302704415484,  0.03452019764641,  0.03582978394968, ...
    0.03693360772938,  0.03781270661771,  0.03845193541496,   0.0388399794176, ...
    0.03897012819567,   0.0388399794176,  0.03845193541496,  0.03781270661771, ...
    0.03693360772938,  0.03582978394968,  0.03452019764641,  0.03302704415484, ...
    0.03137483602187,  0.02959065647781,  0.02770239264592,  0.02573881089455, ...
    0.02372889933714,  0.02170084686872,  0.01968158726936,  0.01769631926666, ...
    0.01576814883838,  0.01391741208812,  0.01216156637634,  0.01051484141936, ...
   0.008988407884082,  0.00759012972115, 0.006324679572635, 0.005193674775058, ...
   0.004196053486178, 0.003328222462836, 0.002584447537782, 0.001957222069436, ...
   0.001437574466143, 0.001015528639767, 0.0006805245562375, 0.0004216319357445, ...
  0.0002281307394882, 8.946081417693e-05, -4.259718936231e-06, -6.227238464387e-05, ...
  -9.276608939232e-05 ,-0.0001030285907574, -9.917286695792e-05, -8.633271245897e-05, ...
  -0.0001573141843083];

adc_output = filter(flt2_num, 1, dac_out)/Vref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot signal before and after ADC
figure; hold all; grid on;
stairs(adc_input_r(end-10000:end), 'r');
stairs(adc_output(end-10000:end),  'b');
legend('ADC input', 'ADC output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot signal spectrum 
[adc_input_r_pwr, f] = periodogram(adc_input_r, [], [], fs);
[adc_output_pwr,  f] = periodogram(adc_output , [], [], fs); 

figure; hold all; grid on;
plot(f, 10*log10(adc_input_r_pwr),  'r');
plot(f, 10*log10(adc_output_pwr),   'b');
legend('ADC input', 'ADC output');


