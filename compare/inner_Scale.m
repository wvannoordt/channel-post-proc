function [dum] = inner_Scale(oxdata, oxcol, csname, titlename, pngname, fignum)
    dum = 0.0;
    
    tau = 20.75;
    
    ox_y = oxdata(:,1)+1.0;
    ox_u = oxdata(:,oxcol);
    csdata = csvread(csname);
    cs_y = csdata(:,1);
    cs_u = csdata(:,2);
    
    cs_mu  = csvread('purdue/cs-mu.csv');
    cs_mu = cs_mu(:,2);
    cs_rho = csvread('purdue/cs-rho.csv');
    cs_rho = cs_rho(:,2);
    ox_mu  = oxdata(:,2);
    ox_rho = oxdata(:,3);
    
    utau_ox = sqrt(tau./ox_rho);
    utau_cs = sqrt(tau./cs_rho);
    ox_yp = ox_rho.*ox_y.*utau_ox./ox_mu;
    cs_yp = cs_rho.*cs_y.*utau_cs./cs_mu;
    
    figure('position', [250 250 1200 800])
    semilogx(ox_yp, ox_u, 'linewidth', 3)
    hold on
    semilogx(cs_yp, cs_u, 'linewidth', 3)
    if (oxcol==9)
      axis([10^0 10^2.5 0 1500])
    end
    %h = legend('WMLES', 'DNS');
    title(titlename);
    saveas(gcf, pngname);
end

