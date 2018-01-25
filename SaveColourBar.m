function SaveColourBar(customcmap,cmaprange,cbarsavefilename)

colourbar_fig = figure('Color','w', ...
                      'Renderer', 'OpenGL', ...
                      'Units', 'inches', ...
                      'PaperUnits', 'inches', ...
                      'PaperSize', [10 3], ...
                      'PaperPositionMode', 'manual', ...
                      'PaperPosition', [0 0 10 3], ...
                      'Visible','off');
plot(0,[ min(cmaprange) max(cmaprange)])
caxis(cmaprange);
colormap(customcmap);
c = colorbar('southoutside');
set(gca,'Visible','off');
axis([0 1 min(cmaprange) max(cmaprange)]);
set(gca,'Position', [0.1 0.5 0.8 0.2],'xtick',[],'YAxisLocation','right','Ticklength', [0 0]);
axis off
box off

% ax = gca;
% axpos = ax.Position;
% cpos = c.Position;
% 
% cpos(4) = 5*cpos(4);
% c.Position = cpos;
% ax.Position = axpos;

[~, titlename] = fileparts(cbarsavefilename);
% titlename = strrep(titlename,'_','\_');
% titlename = strsplit(titlename,'_colourbar');
title(titlename,'Interpreter','none');

set(colourbar_fig,'InvertHardcopy','off')
print(colourbar_fig,'-dpng','-r100',[cbarsavefilename,'.png']);
close(colourbar_fig);

end