function [vmin] = min_dist_pts(pts)
    
   vmin = 10^(16); 
    np   = numel(pts)/2;
    pts2 = [pts,pts(:,1)];
    
    for i=1:np
        for j=1:np
            if (j ~= i && j ~=i+1)
                if (j ~= 1 || i~=np)
                    x1 = pts2(1,i);
                    y1 = pts2(2,i);
                    x2 = pts2(1,i+1);
                    y2 = pts2(2,i+1);
                    x0 = pts2(1,j);
                    y0 = pts2(2,j);
                    
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
                    vmin = min([dtot,vmin]);
                end             
            end
        end
    end
end

