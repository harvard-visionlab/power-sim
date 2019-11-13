function MakePilotGraph_tTests(prefs)
%condition number labels
bar_color = [.75,.75,.75];

conds = unique(prefs.data(:,3));
nConds = length(conds);

for c = 1:nConds
    condition_means(c) = mean(prefs.data(prefs.data(:,3) == c, 2));
end

figure(1)
clf
hold on
subplot(2,1,1)
hold on
title('Pilot Data')
bar(1:nConds, condition_means, 'FaceColor', bar_color)

range = max(condition_means) - min(condition_means);
if range == 0
   range = max(condition_means);
end
drop = .5;
ymin = min(condition_means) - drop*range;
ymax = max(condition_means) + drop*range;
ylim([ymin, ymax])
label_y = min(condition_means) - (drop/2)*range;
xlim([0, nConds + 1])


range = 0:(nConds + 1);
for c = 1:length(range)
    x_label{c} = '';
end
for c = 1:nConds
    x_label{c+1} = prefs.cond_names{c};
end
set(gca,'XTick',range)
set(gca, 'XTickLabel',x_label, 'fontsize',24)


ylabel(prefs.header{2})

%bar labels
fs = 16;
for b = 1:length(condition_means)
    text(b, label_y, num2str(b), 'FontSize', fs);
end

subplot(2,1,2)
title('Power Simulation Parameters', 'FontSize', 24)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'XColor','none')
set(gca,'YColor','none')
xlim([0,1])
ylim([0,1])
spot = .9;
spot_jump = .1;
text(.1, .9, 'Successful Study Requies:', 'FontSize', 16)
count = 0;

%t-tests
for p = 1:size(prefs.comps, 1)
    count = count + 1;
    spot = spot - spot_jump;
    if prefs.within_between == 1
        txt = [num2str(count), ': ', num2str(prefs.comps(p, 1)), ' > ', num2str(prefs.comps(p, 2)), ' (within subjects)'];
    else
        txt = [num2str(count), ': ', num2str(prefs.comps(p, 1)), ' > ', num2str(prefs.comps(p, 2)), ' (between subjects)'];
    end
    text(.2, spot, txt, 'FontSize', 16);
end

%exlcusion critera
if prefs.exclusion_min > -inf || prefs.exclusion_max < inf
 spot = spot - spot_jump;
 text(.1, spot, 'Exclusion Criteria:', 'FontSize', 16);
end

if prefs.exclusion_min > -inf
    spot = spot - spot_jump;
    txt = ['Replace subjects with overall score ', '\leq ', num2str(prefs.exclusion_min)];
    text(.2, spot, txt, 'FontSize', 16);
end

if prefs.exclusion_max < inf
    spot = spot - spot_jump;
    txt = ['Replace subjects with overall score ', '\geq ', num2str(prefs.exclusion_max)];
    text(.2, spot, txt, 'FontSize', 16);
end


%position figure on screen
set(0,'units','pixels')
Pix_SS = get(0,'screensize');
w = prefs.figure_width;
h = prefs.figure_height;
set(gcf, 'Position',  [Pix_SS(3)-1.75*w, Pix_SS(4) - h, w, h])
end