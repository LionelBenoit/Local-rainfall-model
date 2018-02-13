function [] = make_gif( M_for_gif,output_file,speed,max_color_scale )
   
[~,~,st]=size(M_for_gif);


for ep_to_plot=1:st
    display(strcat('create gif, step:',num2str(ep_to_plot),'/',num2str(st)))
    h = figure(4);
    set(gcf,'Visible', 'off'); 
    g=subplot(2,1,1);
    imagesc(M_for_gif(:,:,ep_to_plot));
    %axis square
    axis equal tight
    caxis([0 max_color_scale])
    my_map=colormap('jet');
    my_map(1,:)=[1 1 1];
    colormap(my_map)
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    colorbar
    
    p = get(g,'position');
    set(g, 'position', p.*[-0.4 0.8 1.5 1.5]);
    
    g2=subplot(2,1,2);
    plot(reshape(mean(mean(M_for_gif)),[1,st]),'k')
    
    hold on
    plot(ep_to_plot,mean(mean(M_for_gif(:,:,ep_to_plot))),'r.','MarkerSize',10)
    
    p2 = get(g2,'position');
    set(g2, 'position', p2.*[1 1 1 0.8]);
    
    drawnow;
    
    %capture the plot as an image
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    
    % Write to the GIF File 
    if ep_to_plot == 1
        imwrite(imind,cm,output_file,'gif', 'Loopcount',inf); 
    else
        gif_fps = speed; 
        imwrite(imind,cm,output_file,'gif','WriteMode','append','DelayTime',1/gif_fps); 
    end 
    close 4
end

end