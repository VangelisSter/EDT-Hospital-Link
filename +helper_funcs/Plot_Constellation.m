function Plot_Constellation(symbols, plot_title)
%PLOT_CONSTELLATION Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    symbols
    plot_title
end

figure;
scatter(real(symbols), imag(symbols), 20, 'filled');
grid on;
axis equal;
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('In-Phase');
ylabel('Quadrature');
title(plot_title);

end