%% Delta input parameters
delta_input = [1, zeros(1,9999)];
fs=48e3;
%[delta_BM, ~] = BM_passive(delta_input,fs);
% for V_BMtd , x = no of filter, y = no of time steps

%% Pure sine input parameters
amp = 1;             % Amplitude of the excitation
t = 0:1/fs:0.01;     % 481 time steps
f = 1000;
sine_input = amp * sin(2*pi*f*t);
%[V_BMtd_sine, ~]= BM_passive(sine_input,fs);

%% Combination input parameters
amp = 1;             % Amplitude of the excitation
t = 0:1/fs:0.01;
comb_input = amp * sin(2*pi*1000*t) + cos(2*pi*2000*t) + cos(2*pi*9000*t);
%[V_BMtd_comb, ~]= BM_passive(comb_input,fs);

%% Gamma Filter parameters
no_BM_filters = size(delta_BM,2); % number of membrane sections = no of filters
time_steps = size(delta_BM, 1); % resp for every time step of signal
T = 1/fs;
b = 1.14;
n = delta_input;
low_limit= 1; % start of signal time step
high_limit = size(delta_input, 2); % end of signal time step

%% Freq response of BM passive on delta input

% Need to record the pk gain and centre freq of the response (characteristic desired)

x=(low_limit:high_limit)*fs/(high_limit+1); % freq resp x-axis
fr_matrix = zeros(size(x, 2), no_BM_filters); % store value of freq resp at (x,y) point

% 1. find freq resp of every filters' impulse resp in "delta_BM"
for filter=1:no_BM_filters
    filter_fft = fft(delta_BM(:,filter)); % fft a column of response per filter
    fft_decibels = 20*log10(abs(filter_fft(low_limit:high_limit))); % convert the result into decibels
    fr_matrix(:,filter)= fft_decibels; % store decibels values into freq resp matrix
end

% 2. build matrix of the freq response
% record 1st row: peak gain desired, 2nd row: desired peak gain's central
pk_desired = zeros(2, 1000);
for filter=1:no_BM_filters
    pk_desired(1,filter) = max(fr_matrix(:,filter)); % find the maximum gain of this filter
    time_step_index = find(fr_matrix(:,filter) == pk_desired(1,filter), 1, 'first');
    % find the index where max gain happen
    pk_desired(2,filter) = x(1, time_step_index); %get the freq where pkGain happens
end

%3 build gamma filters

% 1000 filter x their corresponding ir of size high_limit
frGamma_matrix = zeros(size(pk_desired, 2), high_limit);
irGamma_matrix = zeros(size(pk_desired, 2), high_limit);

%for recording every gamma filter's freq resp for delta
for i=1:size(pk_desired, 2) % loop over no of filters
    if pk_desired(1, i) ~= 0 && pk_desired(2, i) ~= 0
        %ignore empty columns which occur if freq resp sampled in intervals
        pk_Gain = pk_desired(1,i); % get desired gain for this filter
        pkGain_freq = pk_desired(2,i); % get the central freq for this desired gain
        [fr_gamma, ir_gamma] = fft_gamma(pkGain_freq, b, T, low_limit, high_limit, pk_Gain);
        % return impulse resp
        frGamma_matrix(i,:) = fr_gamma;
        irGamma_matrix(i,:) = ir_gamma;

    end
end

%% sine input
%ir = impulse resp
% alrdy have delta ir, need to conv with sine input to get filtered output
% no. of gamma filters x sine input
conv_sine = zeros(size(pk_desired, 2), size(sine_input, 2));
for filter_no = 1:size(irGamma_matrix, 1)
    c = conv2(irGamma_matrix(filter_no, :), sine_input, 'same'); % convolution
    conv_sine(filter_no, 1:size(c, 2)) = c;
end

%to plot sine impulse resp of bm Passive

 %              time                  filter                    g[n]
   mesh(1:size(V_BMtd_sine, 1), 1:size(V_BMtd_sine, 2) , transpose(V_BMtd_sine));
   xlabel('time')
   zlabel('position')
   zlim([-1e-5 1e-5])
   ylabel('filter no')
   ylim([300 1000]) % only plotting eligible filters
   title('BM Passive on Pure Sinusoid Input')
   figure()

   mat = zeros(size(conv_sine,1), size(sine_input, 2)); % only take conv result of size of input time steps
   for i=1:size(conv_sine,1)
       mat(i,:) = conv_sine(i,1:size(sine_input, 2));
   end
%to plot sine impulse resp of filterbank
%           filter          time            g[n]
   mesh(1:size(mat, 1), 1:size(mat, 2), transpose(mat));
   xlabel('filter no')
   zlabel('position')
   ylabel('convoluted time steps')
   zlim([-1e-5 1e-5])
   xlim([300 1000])
   title('Filterbank on Pure Sinusoid Input')
   figure()

%% combination input

% no. of gamma filters x combination input
conv_comb = zeros(size(pk_desired, 2), size(comb_input, 2));
for filter_no = 1:size(irGamma_matrix, 1)
    c = conv2(irGamma_matrix(filter_no, :), comb_input, 'same'); % convolution
    conv_comb(filter_no, 1:size(c, 2)) = c;
end

%to plot combination impulse resp of bm Passive

 %    time    filter  g[n]
   mesh(1:size(V_BMtd_sine,1), 1:size(V_BMtd_sine,2) , transpose(V_BMtd_comb));
   xlabel('time')
   zlabel('position')
   zlim([-1e-5 1e-5])
   ylabel('filter no')
   ylim([300 1000])
   title('BM Passive on Combination Input')
   figure()

   mat = zeros(size(conv_comb,1), size(comb_input, 2)); % only take conv result of size of input time steps
   for i=1:size(conv_comb,1)
       mat(i,:) = conv_comb(i,1:size(comb_input, 2));
   end
%to plot sine impulse resp of filterbank
%    filter    time  g[n]
   mesh(1:size(mat, 1), 1:size(mat, 2), transpose(mat));
   xlabel('filter no')
   zlabel('position')
   ylabel('convoluted time steps')
   zlim([-1e-5 1e-5])
   xlim([300 1000])
   title('Filterbank on Combination Input');
   figure()

%% with outer_middle_filter(input1)
ear_input = abs(outer_middle_filter(sine_input));
conv_ear = zeros(size(pk_desired, 2), size(ear_input, 2));
for filter_no = 1:size(irGamma_matrix, 1)
    c = conv2(irGamma_matrix(filter_no, :), ear_input, 'same'); % convolution
    conv_ear(filter_no, 1:size(c, 2)) = c;
end

%to plot sine impulse resp of bm Passive
[V_BMtd_ear, ~]= BM_passive(ear_input,fs);
 %       time    filter  g[n]
   mesh(1:size(V_BMtd_ear, 1), 1:size(V_BMtd_ear, 2) , transpose(V_BMtd_ear));
   xlabel('time')
   zlabel('position')
   zlim([-1e-5 1e-5])
   ylabel('filter no')
   ylim([300 1000])
   title('BM Passive on outer ear filtered input');
   figure()

   mat = zeros(size(V_BMtd_ear, 2), size(ear_input, 2));
   for i=1:size(V_BMtd_ear, 2)
       mat(i,:) = conv_ear(i,1:size(ear_input, 2));
   end
%to plot sine impulse resp of filterbank
%    filter    time  g[n]
   mesh(1:size(mat, 1), 1:size(mat, 2), transpose(mat));
   xlabel('filter no')
   zlabel('position')
   ylabel('convoluted time steps')
   zlim([-1e-5 1e-5])
   xlim([300 1000])
   title('Filterbank on outer ear filtered input');
   figure()



%% Freq response of filterbank and BM Passive (delta)

%to plot delta freq resp of bm Passive
hold on
for i=200:10:no_BM_filters % use for i=200:10:no_BM_filters instead for a more readable graph
    plot(x, fr_matrix(:,i));
end
hold off
xlim([0 1e4]);
ylim([-140 -110]);
xlabel('Frequency(Hz)');
ylabel('Gain 20log(dB)');
title('BM Passive Delta input Frequency Response');
figure()

%to plot delta freq resp of bm Passive
hold on
for i=200:10:size(frGamma_matrix, 1) % use for i=200:10:no_BM_filters instead for a more readable graph
    plot(x,20*log10(abs(frGamma_matrix(i,:))));
end
hold off
xlim([0 9000])
ylim([-140 -110])
xlabel('Frequency(Hz)');
ylabel('Gain 20log(dB)');
title('Gamma Delta input Frequency Response');
