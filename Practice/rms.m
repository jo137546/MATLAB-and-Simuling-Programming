
N = input('Enter the number of frequency-PSD pairs: ');

% Initialize arrays
frequency = zeros(1, N);
psd = zeros(1, N);


for i = 1:N
    frequency(i) = input(['Enter frequency (Hz) for point ', num2str(i), ': ']);
    psd(i) = input(['Enter PSD value (G^2/Hz) for point ', num2str(i), ': ']);
end

Grms = sqrt(trapz(frequency, psd));

% Display the result
fprintf('The Grms value is: %.4f G\n', Grms);