function out = seismic_figs_4d(X, params, amp, fig_switch, svfile)

nrecx=params.nrecx; 
nrecy=params.nrecy; 
nsrcx=params.nsrcx;
nsrcy=params.nsrcy;

if fig_switch(1)==1
    out = reshape(permute(reshape(X,nrecx,nrecy,nsrcx,nsrcy),[1 3 2 4]),nrecx*nsrcx,nrecy*nsrcy);
else
    out = X; 
end

if fig_switch(2)==1 %(1:nsrcx*nrecx,1:nsrcy*nrecy)
    h1 = figure;set(h1, 'Visible', 'off');imagesc(real(out));colormap gray;
    
    if contains(svfile, 'noise')
        caxis([-1, 0]*amp);  
    else
        caxis([-1 1]*amp);
    end
    
    xlabel('source-x, receiver-x','FontSize',18, 'FontName','helvetica','FontWeight','bold');
    ylabel('source-y, receiver-y','FontSize',18, 'FontName','helvetica','FontWeight','bold');
    set(gca,'fontsize',18,'FontName','helvetica','FontWeight','bold');
    saveas(h1, svfile,'jpg');
%     close(h1); 
end


end