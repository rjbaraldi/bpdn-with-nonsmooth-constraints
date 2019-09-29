clear; 
load temp_test

svfile = '/Users/bobby/Dropbox (uwamath)/TextOnly/IEEEversion_revised/missing_r1_7/testing/';

colors = {'b', 'r', 'g', 'c', 'm'};
dot = {'-', '-*', '<-', 'o-', '->'}; 
names = {'Int', 'DN', 'Int+DN' }; 
h1 = figure;
    plot(perc, s_int,[colors{1} dot{2}], 'LineWidth', 1.5)
    pbaspect([5,1,1])
    xlabel('% Sources Omitted', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
    ylabel('SNR', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
    set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'Position',  [500 500 1000 200]);
    set(gcf, 'PaperSize', [1000 200]); 
    set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    set(gcf,'color','w');
    legend(names{1}, 'Location', 'northeast')
    img = frame2im(getframe(h1)); 
    imwrite(img, [svfile, 'int.jpg']); 
close(h1)

h2 = figure;
    semilogx(noise, s_d,[colors{2} dot{3}], 'LineWidth', 1.5)
    pbaspect([5,1,1])
    xlabel('True \sigma', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
    ylabel('SNR', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
    set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold','XScale', 'log');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'Position',  [500 500 1000 250]);
    set(gcf, 'PaperSize', [1000 250]); 
    set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    set(gcf,'color','w');
    legend(names{2}, 'Location', 'northeast')
    img = frame2im(getframe(h2)); 
    imwrite(img, [svfile, 'd.jpg']); 
close(h2)


h3 = figure;
    b = axes('Position', [.1, .2, .8, .7]); 
    set(b, 'Units', 'normalized')
    set(b, 'Color', 'none')
    a = axes('Position', [.1, .2, .8, .7]);
    set(a, 'Units', 'normalized')
    stem(a, perc, noise); 
    set(a, 'xlim', max(perc));
    set(b, 'xlim', max(noise)); 
    plot(noise, sintd,[colors{3} dot{4}], 'LineWidth', 1.5)
    pbaspect([5,1,1])
    xlabel('True \sigma', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
    ylabel('SNR', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
%     set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold');
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf, 'Position',  [500 500 1000 200]);
%     set(gcf, 'PaperSize', [1000 200]); 
%     set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    set(gcf,'color','w'); 
    legend(names{3}, 'Location', 'northeast')
    img = frame2im(getframe(h3)); 
    imwrite(img, [svfile, 'intd.jpg']); 
close(h3)