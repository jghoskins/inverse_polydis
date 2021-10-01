function [vs,nvs_out,iadd] = add_verts(vs_in,nvs,zk,u,targs,angs,nangs,bad)

inds = [];
vs = vs_in;
nadd = 0;
 for i=1:nvs-1
     rad = norm(vs(:,i+1+nadd)-vs(:,i+nadd),'fro');
     if (rad >4*pi/zk)
        disp('here')
        pmid =  (vs(:,i+1+nadd)+vs(:,i+nadd))/2;
        vs = [vs(:,1:(i+nadd)),pmid,vs(:,(i+1+nadd):end)];
        inds = [inds,i+nadd+1];
        nadd = nadd + 1;
     end
 end

rad = norm(vs(:,end)-vs(:,1),'fro');
if (rad >4*pi/zk)
   disp('here')
   pmid =  (vs(:,end)+vs(:,1))/2;
   vs = [vs(:,1:end),pmid];
   nadd = nadd + 1;
   inds = [inds,nvs+nadd];
end

nvs_out = nvs + nadd;

vs
nvs_out

xy_s = sum(vs')'/nvs_out;
[u_s,u_grad,chnkr_s,grad,bd_sol,bd_data_grad,bd_sol_grad,F,err_est] = ...
   helm_dirichlet_solver_wgrad(vs,zk,targs,angs,xy_s);

ugrad_mat = reshape(u_grad,[nangs*nangs,2*nvs_out]);
vdelt     = reshape(u_s-u,[nangs*nangs,1]);
vgrad = reshape(2*real(vdelt'*ugrad_mat),[2,nvs_out]);
ugrads = sum(abs(vgrad).^2,1);
ugrads = sqrt(ugrads)

iadd = 0;
bad
if (numel(inds) >0)
 max = 0;
 i   = 0;
 for ii=1:numel(inds)
     if (sum(bad==ii) == 0)
         if (max < ugrads(inds(ii)))
             max = ugrads(inds(ii));
             i = ii;
         end
     end
 end
 if (i ~= 0)
%[m,i] = max(ugrads(inds))
    iv = inds(i);
    
        iii = iv;
%    if (iii <nvs_out)
%        s1 = ugrads(iii-1);
%        s2 = ugrads(iii+1);
%        ss = s1 + s2;   
%        vs(:,iii) = s1/ss*vs(:,iii-1) + s2/ss*vs(:,iii+1);
%    else
%        s1 = ugrads(iii-1);
%        s2 = ugrads(1);
%        ss = s1 + s2;   
%        vs(:,iii) = s1/ss*vs(:,iii-1) + s2/ss*vs(:,1);    
%    end 
    
    inds(i) = []
    vs(:,inds) = []
    nvs_out = nvs + 1;
    iadd = i;
    
 else
     vs = vs_in;
     nvs_out = nvs;
 end
else
vs = vs_in;
nvs_out = nvs;



end    


