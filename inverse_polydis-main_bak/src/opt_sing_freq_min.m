function [vs,ier,e_new]=opt_sing_freq_min(vs_in,nvs,zk,u,targs, ...
             angs,nangs,rcut)

ier = 0;
vt = [vs_in,vs_in(:,1)];
rx = (vt(1,1:nvs)-vt(1,2:(nvs+1))).^2;
ry = (vt(2,1:nvs)-vt(2,2:(nvs+1))).^2;
rr = sqrt(rx+ry);
rmin = min(rr);

vs = vs_in;

ifflag = false;

ifused = true;
iloop = 0;

while (ifflag == false)
  iloop = iloop + 1;
%  if (iloop >=150)
%      ier = 1;
%      return
%  end
vt = [vs,vs(:,1)];
rx = (vt(1,1:nvs)-vt(1,2:(nvs+1))).^2;
ry = (vt(2,1:nvs)-vt(2,2:(nvs+1))).^2;
rr = sqrt(rx+ry);
rmin2 = min(rr);   
if (rmin2 <rcut)
    ier = 1;
    vs_out = vs_in;
    e_new = 10000;
    ifflag = false;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xy_s = sum(vs')'/nvs;
[u_s,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(vs,zk,targs,angs,xy_s);
plot(chnkr_s)
shg
ugrad_mat = reshape(u_grad,[nangs*nangs,2*nvs]);
vdelt     = reshape(u_s-u,[nangs*nangs,1]);
vgrad = reshape(2*real(vdelt'*ugrad_mat),[2,nvs]);    

vtot = [real(ugrad_mat);imag(ugrad_mat)];
vect = [real(vdelt);imag(vdelt)];
[uu,ss,vv] = svd(vtot,'econ');
ss(ss<10^(-3)*ss(1,1)) = 0;
ss(ss>0) = 1./ss(ss>0);
vguess = real(vv*ss'*uu'*vect);
%vtt = vtot'*vtot;
%vnr = norm(vtt,2);
%vtt = vtt + 10^(-4)*eye(size(vtt))*vnr;
%vgg = vtot'*vect;
%vguess = real(vtt\vgg);
vguess = reshape(vguess,[2,nvs]);

%[vgrad] = get_grad_faster(zk,vs,nvs,angs,targs,u,h);
%[dder]  = get_dder(zk,vs,nvs,vgrad,angs,targs,h,u);
%norm(vgrad,'fro')
%[hess] = get_hess(zk,vs,nvs,angs,targs,u,h);

% try newton
%vguess = reshape(hess,[numel(vgrad),numel(vgrad)])\vgrad(:);
%vguess = reshape(vguess,[2,nvs]);

inewt = true; 
inewted = 0;
ifused = true;
while (inewt == true)  
    
v_update = vs - vguess;
pg = polyshape(v_update','Simplify',false);
dlarge = max(vguess(1,:).^2+vguess(2,:).^2);
dlarge = sqrt(dlarge);
dlarge = 0;
if (issimplified(pg) && dlarge<abs(pi/(4*zk)))    
'newton'
xy_s = sum(v_update')'/nvs;
[u_new,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(v_update,zk,targs,angs,xy_s);
e_old = norm(u-u_s,'fro')
e_new_newt = norm(u-u_new,'fro')
if (e_old >e_new_newt) 
 vs = v_update;  
 u_s = u_new;
 inewt = false;
 ifused = true;
elseif (e_old > e_new_newt*0.9)
  vguess = vguess/2;
  inewted= inewted + 1;
else
   inewted = 8;
end    
else
  inewt = false;
  ifused = false;
end
if (inewted >= 8)
    inewt = false;
    ifused= false;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vgrad_o = vgrad;    
if (ifused == false)
    vgrad_o = vgrad;
    vs_o    = vs;
end
%%%%% finite difference step size
h = 10^(-4);

%disp('beginning gradient - slow')
%[vgrad] = get_grad_faster(zk,vs,nvs,angs,targs,u,h);
%disp('gradient ended - slow')
xy_s = sum(vs')'/nvs;
[u_s,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(vs,zk,targs,angs,xy_s);
ugrad_mat = reshape(u_grad,[nangs*nangs,2*nvs]);
vdelt     = reshape(u_s-u,[nangs*nangs,1]);
vgrad = reshape(2*real(vdelt'*ugrad_mat),[2,nvs]);

vnorm = norm(vgrad,'fro');
vnorm_rel = vnorm/norm(u,'fro')

e_new = norm(u-u_s,'fro');
if (vnorm_rel >= 10^(-3))

[dder]  = get_dder(u_s,zk,vs,nvs,vgrad,angs,targs,h,u);    
    
hstep = -1/dder;
hstep = -abs(hstep);

%if (ifused == false)
%    hstep2 = dot(vs(:)-vs_o(:),vgrad(:)-vgrad_o(:))/dot(vgrad(:)-vgrad_o(:),vgrad(:)-vgrad_o(:));
%    hstep2 = -abs(hstep2);
%    hstep = min([hstep,hstep2]);
%    hstep = min(hstep,-abs(pi/(zk*4)));
%end
vs_o = vs;
vgrad_o = vgrad;

e_new = 1;
e_old = 0;
istep = 0
while (e_new > e_old)
    hstep
    istep = istep + 1;
  if (istep < 5)  
        v_update = vs + vgrad*hstep;
        pg = polyshape(v_update','Simplify',false);
        
        if (issimplified(pg))
        xy_s = sum(v_update')'/nvs;
        [u_new,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(v_update,zk,targs,angs,xy_s);
%[u_new,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(v_update,zk,targs,angs,xy_s);
        e_old = norm(u-u_s,'fro');
        e_new = norm(u-u_new,'fro');
        disp('in grad descent loop')
        e_old
        e_new
        end
  else
      %v_update = vs + vgrad*hstep; 
      ifflag = true;
      e_new = e_old;
      u_new = u_s;
      v_update = vs;
  end
  hstep = hstep/2;
end   


u_s = u_new;
vs = v_update;



[u_new,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(vs,zk,targs,angs,xy_s);
disp('final value')
e_old = norm(u-u_s,'fro')

ugrad_mat = reshape(u_grad,[nangs*nangs,2*nvs]);
vdelt     = reshape(u_s-u,[nangs*nangs,1]);
vgrad = reshape(2*real(vdelt'*ugrad_mat),[2,nvs]);
vnorm = norm(vgrad,'fro');
vnorm_rel = vnorm/norm(u,'fro')
end
if (vnorm_rel < 10^(-3))
    ifflag = true;
end

    ifused = false;
    
end    
end

