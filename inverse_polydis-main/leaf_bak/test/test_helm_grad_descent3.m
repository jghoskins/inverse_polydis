if (1>0)
addpath('../src');
nangs = 40;
[verts,xyin,angs,targs] = init_shape(nangs);

nvs = 4;
r0  = 2;
[vs] = init_guess(nvs,xyin,r0);

nvs_out = nvs;
zk = 1.1;

chnkrs = [];
err_curr = 10000;
end

while (zk<20)
disp(zk)

tic
[u,chnkr,bd_sol,F,err1] = helm_dirichlet_solver(verts,zk,targs,angs,xyin);
toc

xy_s = sum(vs')'/nvs;
[u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vs,zk,targs,angs,xy_s);
chnkrs = [chnkrs, chnkr_s];


err_curr = norm(u-u_s,'fro');

vs_in = vs;
rcut = 0;
rcut  = abs(pi/(2*zk));
[vs_out,ier,e_new] = opt_sing_freq_min(vs_in,nvs,zk,u,...
    targs,angs,nangs,rcut);
%[vs_out,e_new] = opt_sing_freq(vs_in,nvs,zk,u,targs,angs,...
%         nangs);
%[vs_out,ier,e_new] = opt_sing_freq_min(vs_in,nvs,zk,u,targs,angs,...
%         nangs,rcut);
if (ier == 0)
    vs = vs_out;
%    if (e_new < err_curr)
%        vs = vs_out;
        err_curr = e_new;
%    end    
end     
%vs = vs_out;

iffixed = false;
bad = [];
while (iffixed == false)
[vs_out,nvs_out,i2add] = add_verts(vs,nvs,zk,u,targs,angs,nangs,bad);

vs_in = vs_out;
nvs_in= numel(vs_out)/2;
rcut  = abs(pi/(2*zk));
[vs_out,ier,e_new] = opt_sing_freq_min(vs_in,nvs_in,zk,u,...
    targs,angs,nangs,rcut);
if (ier == 0)
    bad = [];
    vs = vs_out;
    if (nvs_out == nvs)
        iffixed = true;
    end
    nvs = nvs_out;
else
    vs_out = vs;
    nvs_out = numel(vs_out)/2;
    bad = [bad, i2add];
    if (i2add == 0) 
        iffixed = true;
    end
    %iffixed = true;
%    if (e_new < err_curr)
%        vs = vs_out;
%        err_curr = e_new;
%        if (nvs_out == nvs)
%            iffixed = true;
%        end
%        nvs = nvs_out;
%    end    
end
end

zk = zk + 0.25;
end