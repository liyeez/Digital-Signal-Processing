

function [ir1_fft, ir1_new] = fft_gamma(f1,b,T, low_limit, high_limit, bm_pkGain)
    n = low_limit:high_limit;
    % Calculate denormalised impulse responses
    ir1 = (n*T).^3.*exp(-2*pi*b*(24.7+0.108*f1)*n*T).*cos(2*pi*f1*n*T);
    ir1_fft = fft(ir1);

    a1 = 1/max(abs(ir1_fft)); % normalize
    a2 = 10^(bm_pkGain/20);   % adjust gain to desired pk gain
    % repeat same process after multiplying constant....
    ir1_new = a1*a2*ir1;

    % Calculate frequency responses again...
    ir1_fft = fft(ir1_new);

end
