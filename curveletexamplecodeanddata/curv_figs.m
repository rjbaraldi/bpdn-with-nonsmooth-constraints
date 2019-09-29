function curv_figs(X, nt, ns, svtitle, set_switch)



    
    h1 = figure;set(h1, 'Visible', 'off');
    
    if strcmp(set_switch, 'data')
        imagesc(X,[-0.5,0.5]);
    else
        imagesc(reshape(X,nt,ns),[-1,0]);
    end

    colormap gray;
    set(gca, 'FontSize',18, 'FontName','helvetica','FontWeight','bold', 'xtick',[], 'xticklabel',[], 'ytick', [], 'yticklabel', [])
    set(gcf, 'PaperPositionMode', 'manual', 'PaperUnits', 'inches', 'Position', [200 200 1000 1000], 'color','w');
    set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    saveas(h1, svtitle, 'jpg');
    close(h1)


end