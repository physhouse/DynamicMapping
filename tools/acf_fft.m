load vx.dat
spec = fft(vx);
mag = abs(spec);
psd = mag.*mag / (2*pi*length(vx));
acf = ifft(psd) * 2 * pi;
x = 1:1:length(vx);
plot(x, acf);