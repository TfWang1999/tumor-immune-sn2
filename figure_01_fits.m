function figure_01_fits()
% Plot the saved model fits; this entry never fits or integrates.
root=fileparts(mfilename('fullpath'));
saved=load(fullfile(root,'data','fits','fits.mat'),'bundle');
out=fullfile(root,'outputs');if ~isfolder(out),mkdir(out);end
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman', ...
    'defaultAxesFontSize',9,'defaultTextFontSize',10,'defaultAxesLineWidth',.75);
draw_fits(saved.bundle,root,out);
end

function draw_fits(b,folder,figureFolder)
fig=paper_figure(180,130);tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j=1:4
    d=b.data(j);f=b.fits{j,2};ax=nexttile(tl);hold(ax,'on');
    % Recover the unweighted residual sum of squares from the saved RMSE.
    % This uses the original fitted observation residuals, not interpolation
    % of the dense display curve, and does not run the model again.
    residualSumSquares=numel(d.v)*f.rmse^2;
    totalSumSquares=sum((d.v-mean(d.v)).^2);
    if totalSumSquares>0
        rSquared=1-residualSumSquares/totalSumSquares;
    else
        rSquared=NaN;
    end
    observations=errorbar(ax,d.t,d.v,d.sd,'ko','MarkerFaceColor','w','MarkerSize',4,'LineWidth',.8);
    curve=readtable(fullfile(folder,'data','fits',[d.id,'_curves.csv']));
    modelFit=plot(ax,curve.day,curve.C_mm3,'Color',[36,99,162]/255,'LineWidth',1.35);
    xlabel(ax,'Day');ylabel(ax,'Tumour volume (mm^3)');box(ax,'on');
    title(ax,sprintf('(%c) %s','a'+j-1,d.id),'FontWeight','normal');
    xlim(ax,[d.t(1)-.3,d.t(end)+.3]);ylim(ax,[0,max([d.v+d.sd;curve.C_mm3])*1.06]);
    % Preserve supplied error bars; the MC38 source reports SEM, not SD.
    % Dataset-specific error conventions are described in the manuscript.
    legend(ax,[observations,modelFit],{'Control data','Model fit'}, ...
        'Location','northwest','Orientation','vertical','Box','off','FontSize',8.5);
    text(ax,.20,.75,{sprintf('R^2 = %.4f',rSquared), ...
        sprintf('RMSE = %.2f mm^3',f.rmse)},'Units','normalized', ...
        'VerticalAlignment','top','FontSize',8.5,'Interpreter','tex');
    fprintf('%s: RMSE = %.8f, R^2 = %.8f\n',d.id,f.rmse,rSquared);
end
save_figure(fig,figureFolder,'Figure_01');
end

function fig=paper_figure(w,h)
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2,2,w/10,h/10]);
end

function save_figure(fig,folder,name)
exportgraphics(fig,fullfile(folder,[name,'.pdf']),'ContentType','vector','BackgroundColor','white');
print(fig,fullfile(folder,[name,'.svg']),'-dsvg');
exportgraphics(fig,fullfile(folder,[name,'.png']),'Resolution',300);
close(fig);
end
