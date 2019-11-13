function prefs = MakePilotGraph_ANOVA(prefs)


%get condition means
nConds = prefs.f1_num_levels * prefs.f2_num_levels;
sub_means = nan(length(unique(prefs.sub_nums)), nConds);


for s = 1:length(unique(prefs.sub_nums))
    cond = 0;
    for f1 = 1:prefs.f1_num_levels
        for f2 = 1:prefs.f2_num_levels
            cond = cond + 1;
            sub_means(s, cond) = mean(prefs.data(prefs.data(:,1) == s & ...
                prefs.data(:,3) == f1 & prefs.data(:,4) == f2, 2));
        end
    end
end

prefs.sub_means = sub_means;

%condition number labels
color_matrix = [.75 .75 .75;
    1 1 1;
    1 0 0;
    0 1 0;
    0 0 1;
    0 0 0;
    0 1 1];

count = 0;
xspot = 1;
graphx = [];
graph_colors = [];
mid_spot = [];
condition_means = nanmean(sub_means);

for f1 = 1:length(prefs.f1_names)
    for f2 = 1:length(prefs.f2_names)
        count = count + 1;
        level_f1(count) = f1;
        level_f2(count) = f2;
        graph_colors = [graph_colors; color_matrix(f2, :)];
    end
    mid_spot = [mid_spot, mean(xspot:xspot+length(prefs.f2_names)-1)];
    graphx = [graphx, xspot:xspot+length(prefs.f2_names)-1];
    xspot = xspot + length(prefs.f2_names) + 1;
end

figure(1)
clf

clf
hold on
subplot(2,1,1)
hold on
title('Pilot Data')
for c = 1:count
    bar(graphx(c), condition_means(c), 'FaceColor', graph_colors(c,:))
end

range = max(condition_means) - min(condition_means);
if range == 0
    range = max(condition_means);
end
drop = .5;
ymin = min(condition_means) - drop*range;
ymax = max(condition_means) + drop*range;
ylim([ymin, ymax])
label_y = min(condition_means) - (drop/2)*range;
xlim([min(graphx) - 1, max(graphx) + 1])

if rem(length(prefs.f2_names), 2) == 1 %odd amount of f2 levels
    range = min(graphx) - 1: max(graphx) + 1;
    for c = 1:length(range)
        x_label{c} = '';
    end
    for c = 1:length(prefs.f1_names)
        x_label{range == mid_spot(c)} = prefs.f1_names{c};
    end
    set(gca,'XTick',range)
    set(gca, 'XTickLabel',x_label, 'fontsize',24)
else
    range = min(graphx) - 1:.5: max(graphx) + 1;
    for c = 1:length(range)
        x_label{c} = '';
    end
    for c = 1:length(prefs.f1_names)
        x_label{range == mid_spot(c)} = prefs.f1_names{c};
    end
    set(gca,'XTick',range)
    set(gca, 'XTickLabel',x_label, 'fontsize',24)
end

legend(prefs.f2_names, 'FontSize', 16)

ylabel(prefs.header{2})

%bar labels
fs = 16;
for b = 1:length(condition_means)
    text(graphx(b), label_y, num2str(b), 'FontSize', fs);
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
%ME 1
if prefs.sig_ME1
    count = count + 1;
    spot = spot - spot_jump;
    
    txt = [num2str(count), ': Main Effect ', prefs.header{3}];
    
    text(.2, spot, txt, 'FontSize', 16);
end

%ME2
if prefs.sig_ME2
    count = count + 1;
    spot = spot - spot_jump;
    
    txt = [num2str(count), ': Main Effect ', prefs.header{4}];
    text(.2, spot, txt, 'FontSize', 16);
end

%interaction
if prefs.sig_int
    count = count + 1;
    spot = spot - spot_jump;
    
    txt = [num2str(count), ': Interaction of ',  prefs.header{3}, ' x ', prefs.header{4}];
    
    text(.2, spot, txt, 'FontSize', 16);
end

%t-tests
for p = 1:size(prefs.comps, 1)
    count = count + 1;
    spot = spot - spot_jump;
    
    %extra text
    if prefs.between_within_mixed == 3
        if level_f1(prefs.comps(p, 1)) == level_f1(prefs.comps(p, 2))
            extra_text = ' (within subjects)';
        else
            extra_text = ' (between subjects)';
        end
    elseif prefs.between_within_mixed == 2
        extra_text = ' (within subjects)';
    elseif prefs.between_within_mixed == 1
        extra_text = ' (between subjects)';
    end
    
    txt = [num2str(count), ': ', num2str(prefs.comps(p, 1)), ' > ', num2str(prefs.comps(p, 2)), extra_text];
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