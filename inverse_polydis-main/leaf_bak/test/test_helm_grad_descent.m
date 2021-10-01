addpath('../src');
nangs = 40;
zk = 1.1;

if (1>0)
[verts,xyin,angs,targs] = init_shape(nangs);
tic
[u,chnkr,bd_sol,F,err1] = helm_dirichlet_solver(verts,zk,targs,angs,xyin);
toc
end

nvs = 6;
r0  = 2;
[vs] = init_guess(nvs,xyin,r0);
xy_s = sum(vs')'/nvs;
[u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vs,zk,targs,angs,xy_s);

norm(u-u_s,'fro')

vstart = vs;

vs_found = zeros([size(vs),100]);

for i=1:100
    
i    
vs_found(:,:,i) = vs;    
%%%%% finite difference step size
h = 10^(-4);

if (i >0)

xy_s = sum(v_update')'/nvs;
[u_s,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(vs,zk,targs,angs,xy_s);
ugrad_mat = reshape(u_grad,[nangs*nangs,2*nvs]);
vdelt     = reshape(u_s-u,[nangs*nangs,1]);
%vgrad = reshape(2*real(vdelt'*ugrad_mat),[2,nvs]);    

vtot = [real(ugrad_mat);imag(ugrad_mat)];
vect = [real(vdelt);imag(vdelt)];
vguess = vtot\vect;
vguess = reshape(vguess,[2,nvs]);

%[vgrad] = get_grad_faster(zk,vs,nvs,angs,targs,u,h);
%[dder]  = get_dder(zk,vs,nvs,vgrad,angs,targs,h,u);
%norm(vgrad,'fro')
%[hess] = get_hess(zk,vs,nvs,angs,targs,u,h);

% try newton
%vguess = reshape(hess,[numel(vgrad),numel(vgrad)])\vgrad(:);
%vguess = reshape(vguess,[2,nvs]);
v_update = vs - vguess;
pg = polyshape(v_update');

if (issimplified(pg))
'newton'
xy_s = sum(v_update')'/nvs;
[u_new,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(v_update,zk,targs,angs,xy_s);
e_old = norm(u-u_s,'fro')
e_new_newt = norm(u-u_new,'fro')
if (e_old >e_new_newt) 
 vs = v_update;  
 u_s = u_new;
end
end
end

%disp('beginning gradient - slow')
%[vgrad] = get_grad_faster(zk,vs,nvs,angs,targs,u,h);
%disp('gradient ended - slow')
[u_s,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(vs,zk,targs,angs,xy_s);
ugrad_mat = reshape(u_grad,[nangs*nangs,2*nvs]);
vdelt     = reshape(u_s-u,[nangs*nangs,1]);
vgrad = reshape(2*real(vdelt'*ugrad_mat),[2,nvs]);

[dder]  = get_dder(zk,vs,nvs,vgrad,angs,targs,h,u);
vnorm = norm(vgrad,'fro')
vnorm_rel = vnorm/norm(u,'fro')

hstep = -1/dder
hstep = -abs(hstep);

e_new = 1;
e_old = 0;
while (e_new > e_old)
v_update = vs + vgrad*hstep;

pg = polyshape(v_update');
if (issimplified(pg))
xy_s = sum(v_update')'/nvs;
[u_new,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(v_update,zk,targs,angs,xy_s);
e_old = norm(u-u_s,'fro')
e_new = norm(u-u_new,'fro')
end

hstep = hstep/2;
end
u_s = u_new;
vs = v_update;

end
