function l1deno_fig_create(dat, col,nam, svfile, typeset)
h1 = figure;set(h1, 'Visible', 'off')
switch(typeset)
    
    case{'data'}
        h1 = figure;set(h1, 'Visible', 'off')
        plot(dat, col, 'LineWidth', 1.5, 'MarkerSize', 2.5)
        pbaspect([5,1,1])
        xlim([0, 512])
        ylim([-1, 1])
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold', 'xtick',[], 'xticklabel',[])
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'Position', [500 500 1000 200]);
        set(gcf, 'PaperSize', [1000 200]); 
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        set(gcf,'color','w');
        img = frame2im(getframe(h1)); 
        imwrite(img, [svfile, nam '_results.jpg']);
        
    case{'w1'}
        h1 = figure;set(h1, 'Visible', 'off')
        plot(dat, col, 'LineWidth', 1.5,'MarkerSize', 2.5)
        pbaspect([5,1,1])
        xlim([0, 512])
        ylim([-1, 1])
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold', 'xtick',[], 'xticklabel',[])
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'Position', [500 500 1000 200]);
        set(gcf, 'PaperSize', [1000 200]); 
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        set(gcf,'color','w');
        img = frame2im(getframe(h1)); 
        imwrite(img, [svfile, nam '_w1.jpg']);
        
    
    case{'w2'}
        plot(dat, col, 'LineWidth', 1.5,'MarkerSize', 2.5)
        pbaspect([3,1,1])
        xlim([0, 120])
        if strcmp(nam, 'obs')||strcmp(nam, 'True')
            ylim([-.35, 2.6])
        else
            ylim([-2.6, .35])
        end
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold', 'xtick',[], 'xticklabel',[])
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'Position', [500 500 925 290], 'color','w');
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        imwrite(frame2im(getframe(h1)), [svfile, nam, '_res.jpg']);
        
    
    case{'werr'}
        
        semilogy(dat, 'LineWidth', 2,'MarkerSize', 2.5);
        xlabel('k', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        ylabel('l_2(\cdot)^2', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        legend('||w_1-x||^2', '||A(x) - b - w_2||^2', 'Location', 'SouthWest')
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold','XTickMode', 'auto', 'XTickLabelMode', 'auto')
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'Position', [500 500 925 290], 'color','w');
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        imwrite(frame2im(getframe(h1)), [svfile,nam, '_werr.jpg']);
        
    case{'witer'}
       
        semilogy(dat, 'LineWidth', 2,'MarkerSize', 2.5);
        xlabel('k', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        ylabel('||z^+ - z||^2', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        legend('||w_1^+-w_1||^2', '||w_2^+ - w_2||^2','||x^+ - x||^2', 'Location', 'NorthWest')
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold','XTickMode', 'auto', 'XTickLabelMode', 'auto')
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'Position', [500 500 925 500], 'color','w');
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        imwrite(frame2im(getframe(h1)), [svfile,nam, '_err.jpg']);

    case{'obj'}

        plot(dat, 'LineWidth', 2,'MarkerSize', 2.5);
        xlabel('k', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        ylabel('\phi(\cdot)', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        legend('||w_1||_1', '||x||_1', 'Location', 'NorthEast')
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold','XTickMode', 'auto', 'XTickLabelMode', 'auto')
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'Position', [500 500 925 290], 'color','w');
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        imwrite(frame2im(getframe(h1)), [svfile,nam, '_obj.jpg']);
        
    case{'const'}

        loglog(dat, 'LineWidth', 2,'MarkerSize', 2.5);
        xlabel('k', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        ylabel('Parameter Value', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        legend('\alpha', '\eta_1','\eta_2', 'Location', 'NorthEast')
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold','XTickMode', 'auto', 'XTickLabelMode', 'auto')
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'Position', [500 500 925 290], 'color','w');
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        imwrite(frame2im(getframe(h1)), [svfile,nam, '_param.jpg']);
        
    case{'feas'}
        %plot feasibility

        loglog(dat, 'LineWidth', 2,'MarkerSize', 2.5);
        xlabel('k', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        ylabel('\psi(\cdot)', 'FontSize',18, 'FontName','helvetica','FontWeight','bold')
        legend('||A(x) - b||_p - \sigma','||w_2||_p - \sigma',  'Location', 'NorthEast')
        set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold','XTickMode', 'auto', 'XTickLabelMode', 'auto')
        set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'Position', [500 500 925 290], 'color','w');
        set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        imwrite(frame2im(getframe(h1)), [svfile,nam, '_feas.jpg']);
end
close(h1)



end