function Plot_RF_Signal(t,signal, plot_title)
%PLOT_RF_SIGNAL We plot a signal in the time domain
%   Detailed explanation goes here
arguments (Input)
    t
    signal
    plot_title
end

figure;
plot(t*1e9, signal, 'r', 'LineWidth', 1.2);
grid on;
xlabel('Χρόνος (ns)'); ylabel('Πλάτος');
title(plot_title);
xlim([20, 30]);
end