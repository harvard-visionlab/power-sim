function PlotHeatMap(power_results, prefs)
%plot
figure(2)
clf

power = power_results.power;
power = round(power, 2);

%position figure on screen
set(0,'units','pixels')
Pix_SS = get(0,'screensize');
w = prefs.figure_width;
h = prefs.figure_height;
set(gcf, 'Position',  [Pix_SS(3)-w*.75, Pix_SS(4) - h*.5, w*.75, .6*h])
HeatMapDownloaded(power, power_results.sub_vector, power_results.trial_vector, true, 'GridLines', '-', 'FontSize', 15, 'TickFontSize', 14);
xlabel('Total # of Subjects', 'FontSize', 20)
if prefs.varied_sim_trials
    ylabel('# of Trials in Condition 1', 'FontSize', 20)
else
    ylabel('# of Trials Per Condition', 'FontSize', 20)
end
title('Power by N and # of Trials', 'FontSize', 20)

end
