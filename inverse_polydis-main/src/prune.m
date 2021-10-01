function [vs_out,nvs_out] = prune(vs,rcut)
    nvs = numel(vs)/2;
    vts = [vs(:,nvs),vs,vs(:,1)];
    inds = [];
    for i=1:nvs
  
                    x1 = vts(1,i);
                    y1 = vts(2,i);
                    x2 = vts(1,i+2);
                    y2 = vts(2,i+2);
                    x0 = vts(1,i+1);
                    y0 = vts(2,i+1);
                    
                    xx = (x2-x1);
                    yy = (y2-y1);
                    rr = sqrt(xx^2+yy^2);
                    xx = xx/rr;
                    yy = yy/rr;
                    dpar = (x0-x1)*xx+(y0-y1)*yy;
                    dper = (x0-x1)^2+(y0-y1)^2 - dpar^2;
                    dper = sqrt(dper);
                    if (dpar >rr)
                        dpar = dpar - rr;
                    elseif (dpar < 0)
                        dpar = - dpar;
                    else
                        dpar = 0;
                    end
                    dtot = sqrt(dper^2+dpar^2);
                    %vmin = min([dtot,vmin]);
                    if (dper < rcut)
                        inds = [inds,i];
                    end
    end
    vs_out = vs;
    vs_out(:,inds) = [];
    nvs_out = numel(vs_out)/2;
end

