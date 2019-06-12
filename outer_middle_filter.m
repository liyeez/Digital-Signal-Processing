

function [middle_y] = outer_middle_filter(input1)
%clear all; clc;
% provided an input, filter through learningA transfer functions:

input_fs = 48e3; % sampling frequency
input_t = 0:1/input_fs:0.01;
%filter_input = amp *sin(2*pi*1000*t);
% ------------ outer ear ------------
% zeros: 0Hz, 10kHZ
% poles: 2kHz

% estimations made for the radius of the poles and zeros:

outer_z1 = 0.93;
outer_z2 = 0.99*exp(j*2*pi*(10/48));
outer_p1 = 0.96*exp(j*2*pi*(2/48));
%
% -----------------------------------

% ------------ middle ear ------------
% zeros: 0Hz, 12kHZ
% poles: 900Hz

% estimations made for the radius of the poles and zeros:

middle_z1 = 0.98;
middle_z2 = 0.7*exp(j*2*pi*(12/48));
middle_p1 = 0.98*exp(j*2*pi*(0.9/48));

% -----------------------------------

K_out = 10^(22/20);
K_mid =10^(28/20);

% generate 1000 pts on zplane unit circle between 100Hz and 10kHz
theta = linspace(2*pi*0.1e3/input_fs,2*pi*10e3/input_fs,1000);
z = exp(1j*theta);

% calculate magnitude(db) vs freq(kHz) in frequency response
H_mid = MagResp(z,middle_z1,middle_z2,middle_p1,K_mid);
H_out = MagResp(z,outer_z1,outer_z2,outer_p1,K_out);
H_cascade = H_mid.*H_out;

% collecting all poles and zeros....
middle_z = [middle_z1, middle_z2, conj(middle_z2)];
middle_p = [middle_p1, conj(middle_p1)];
outer_z =[outer_z1, outer_z2, conj(outer_z2)];
outer_p = [outer_p1, conj(outer_p1)];

% build filter from coefficients
outer_y = filter(outer_z, outer_p, input1); % pass in input through outer ear filter
middle_y = filter(middle_z,middle_p, outer_y); %  then pass in outer ear's output through middle ear filter

% uncomment to plot the filtered output
%{
plot (input1);
title('original signal');
ylabel('Amplitude');
xlabel('time(s)');
figure
plot(abs(outer_y));
title('outer ear filtered signal');
ylabel('Amplitude');
xlabel('time(s)');
figure
plot(abs(middle_y));
title('middle ear filtered signal');
ylabel('Amplitude');
xlabel('time(s)');
%}

% plot magnitude response graph
semilogx(theta*input_fs/(2*pi),20*log10(abs(H_mid)));
hold on
semilogx(theta*input_fs/(2*pi),20*log10(abs(H_out)));
semilogx(theta*input_fs/(2*pi),20*log10(abs(H_cascade)));

xlabel('Frequency(Hz)')
ylabel('Gain(dB)')
legend('Middle Ear','Outer Ear','Cascade')
hold off


end
