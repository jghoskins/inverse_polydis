function [u,u_grad,varargout] = helm_dirichlet_solver_wgrad(verts,zk, ...
   targs,angs,xyin,cparams)
%HELM_DIRICHLET_SOLVER solves the helmholtz Dirichlet value problem
%for an analytic solution and a collection of plane waves for a rounded
%polygon and returns the Frechet derivative of the scattered fields
%as well with respect to the vertex coordinates
%
% Syntax: [u,u_grad,varargout] = helm_dirichlet_solver_wgrad(verts,zk,
%       targs,angs,xyin,cparams)
% Input:
%   verts = xy coordinates of vertices describing the polygon (2,nverts)
%   zk - Helmholtz wave number
%   targs - xy coordinates of targets in exterior for evaluating potential
%   (2,nt)
%   angs - angles of incidence for the incoming plane waves (nangs,1)
%   xyin - location of interior point 
% Optional input:
%   cparams - options structure
%       cparams.rounded = true if corner rounding is to be used.
%                         false if no rounding is used (true)
%       cparams.autowidths = automatically compute widths (true)
%       cparams.autowidthsfac = if using autowidths, set widths
%                             to autowidthsfac*minimum of adjoining
%                             edges (0.4)
% 	    cparams.eps - resolve curve to tolerance eps
%                    resolve coordinates, arclength,
%          	     and first and second derivs of coordinates
%		         to this tolerance (1.0e-6)
% output:
%   u - potential at targets (nt,nangs)
%   u_grad - frechet derivative of the scattered field 
%       evaluated at the targets (nt,nangs,2,nv)
%   optional output arguments:
%   chnkr - discretized chunker used in solver
%   grad - Frechet derivative of the shape
%   bd_sol - solution of integral equation corresponding to incident plane
%       waves (n,nangs)
%   bd_data_grad - boundary data used for the frechet derivatives (n,nangs,2,nv)
%   bd_sol_grad - solution of integral equation corresponding to boundary
%       data for the frechet derivatives (n,nangs,2,nv)
%   F - structure containing compressed linear systems
%      F.C - combined field representation (D + i*S)
%      F.S - single layer potential
%      F.Sprime - derivative of single layer potential
%   err - estimated error in solving dirichlet problem
%

if(nargin<=5)
    cparams_use = [];
    cparams_use.rounded = true;
    cparams_use.autowidths = true;
    cparams_use.autowidthsfac = 0.4;
    cparams_use.eps = 1e-6;
else
    cparams_use = cparams;
end
pref.k = 16;
pref.dim = 2;
% Get chunk info

chnkr = chunkerpoly(verts,cparams_use,pref);
ref_opt = struct();
ref_opt.maxchunklen = pi/abs(zk);
chnkr = refine(chnkr,ref_opt);
rn = normals(chnkr);
rn = reshape(rn,[2,chnkr.k*chnkr.nch]);

grad = vertgrad(chnkr);

opts_chunkerkerneval = [];
opts_chunkerkerneval.flam = false;
opts_chunkerkerneval.forcesmooth = true;

% Compress all the relevant operators
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'C',1);
dval = 0.5;
opts_flam = [];
opts_flam.flamtype = 'rskelf';
opts_flam.rank_or_tol = 1e-10;
F = [];
F.C = chunkerflam(chnkr,fkern,dval,opts_flam);


fkern_s = @(s,t) chnk.helm2d.kern(zk,s,t,'s',1);
dval_s = 0.0;
F.S = chunkerflam(chnkr,fkern_s,dval_s,opts_flam);


fkern_sp = @(s,t) chnk.helm2d.kern(zk,s,t,'sprime',1);
dval_sp = -0.5;
F.Sprime = chunkerflam(chnkr,fkern_sp,dval_sp,opts_flam);



% Test exact solution
srcinfo = []; srcinfo.r = xyin; targinfo = []; targinfo.r = chnkr.r;
targinfo.r = reshape(targinfo.r,2,chnkr.k*chnkr.nch);

bd_test = chnk.helm2d.kern(zk,srcinfo,targinfo,'s');

bd_sol_ex = rskelf_sv(F.C,bd_test);


utest = chunkerkerneval(chnkr,fkern,bd_sol_ex,targs,opts_chunkerkerneval);

targinfo = [];
targinfo.r = targs;
uex = chnk.helm2d.kern(zk,srcinfo,targinfo,'s');

err_exact = norm(utest-uex,'fro')/norm(uex,'fro');


% Get the incident fields
rval = chnkr.r;
rval = reshape(rval,2,chnkr.k*chnkr.nch);
xval = rval(1,:);
yval = rval(2,:);

rnx = rn(1,:)'; rny= rn(2,:)';
uinc = -exp(1i*zk*(xval'*cos(angs)' + yval'*sin(angs)'));
dudninc = -uinc.*1i.*zk.*(rnx.*cos(angs)' + rny.*sin(angs)');

% solve integral equation
bd_sol = rskelf_sv(F.C,uinc);
nt = length(targs);
nangs = length(angs);
u = complex(zeros(nt,nangs));
[~,nv] = size(verts);

% Compute scattered field at target locations
for i=1:nangs
    u(:,i) = chunkerkerneval(chnkr,fkern,bd_sol(:,i),targs,opts_chunkerkerneval);
end

% Compute neumann data corresponding to incident fields
[n,~] = size(uinc);
dp = complex(zeros(n,nangs));
for i=1:nangs
    dp(:,i) = helm_dprime(zk,chnkr,F.S,bd_sol(:,i));
end
sp = rskelf_mv(F.Sprime,bd_sol);
dudn = dp + 1i*sp;
dudn = dudn+dudninc;


u_grad = complex(zeros(nt,nangs,2*nv));
bd_sol_grad = complex(zeros(n,nangs,2*nv));
bd_data_grad = complex(zeros(n,nangs,2*nv));

for i=1:2*nv
    g = grad(:,i);
    gx = g(1:2:end); gy = g(2:2:end);
    gdot = (gx.*rnx + gy.*rny);
    gdotrep = repmat(gdot,[1,nangs]);
    bd_data_grad(:,:,i) = -gdotrep.*dudn;
    bd_sol_grad(:,:,i) = rskelf_sv(F.C,bd_data_grad(:,:,i));
    for j=1:nangs
        u_grad(:,j,i) = chunkerkerneval(chnkr,fkern,bd_sol_grad(:,j,i),targs,opts_chunkerkerneval);
    end
end

u_grad = reshape(u_grad,[nt,nangs,2,nv]);
bd_sol_grad = reshape(bd_sol_grad,[n,nangs,2,nv]);
bd_data_grad = reshape(bd_data_grad,[n,nangs,2,nv]);

varargout{1} = chnkr;
varargout{2} = grad;
varargout{3} = bd_sol;
varargout{4} = bd_data_grad;
varargout{5} = bd_sol_grad;
varargout{6} = F;
varargout{7} = err_exact;
   
end