function [dum] = create_fig(oxdata, oxcol, csname, titlename, pngname)
    dum = 0.0;
    ox_y = oxdata(:,1);
    ox_f = oxdata(:,oxcol);
    csdata = load(csname);
    cs_y = csdata(:,1)-1.0;
    cs_f = csdata(:,2);
    figure
    hold on
    plot(ox_y, ox_f)
    plot(cs_y, cs_f)
    h = legend('WMLES', 'DNS');
    title(titlename);
    saveas(gcf, pngname);
end

