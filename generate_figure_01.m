function generate_figure_01()
% Reintegrate the published fitted model; no optimization is performed.
root=fileparts(mfilename('fullpath'));
inputs=jsondecode(fileread(fullfile(root,'data','parameters.json')));
observations=readtable(fullfile(root,'data','observations.csv'),'TextType','string');
out=fullfile(root,'data','fits');if ~isfolder(out),mkdir(out);end
ids={'LLC','MC38','A375','CT26'};
bundle=struct('data',[],'fits',{cell(4,2)});metrics=table();parameters=table();
for j=1:4
    id=ids{j};v=observations(observations.dataset==string(id) & observations.used==1,:);
    d=struct('id',id,'t',v.day,'v',v.mean_mm3,'sd',v.sd_mm3);
    p=inputs.fits.(id).dimensional;u0=inputs.fits.(id).control_initial(:);
    opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',1,'NonNegative',1:3, ...
        'Jacobian',@(t,u)tumor_immune_dc_jacobian(u,p));
    duration=d.t(end)-d.t(1);
    sol=ode15s(@(t,u)tumor_immune_dc_rhs(t,u,p),[0,duration],u0,opts);
    U=deval(sol,d.t-d.t(1)).';dense=linspace(0,duration,401).';curve=deval(sol,dense).';
    residual=U(:,1)-d.v;scale=max([d.sd,.05*d.v,ones(numel(d.v),1)],[],2);
    f=struct('rmse',sqrt(mean(residual.^2)),'r2',1-sum(residual.^2)/sum((d.v-mean(d.v)).^2));
    bundle.data=[bundle.data;d];bundle.fits{j,2}=f;
    writetable(table(dense+d.t(1),curve(:,1),curve(:,2),curve(:,3), ...
        'VariableNames',{'day','C_mm3','E_per_day','Z'}),fullfile(out,[id,'_curves.csv']));
    metrics=[metrics;table(string(id),f.rmse,sum((residual./scale).^2),f.r2,numel(d.t), ...
        'VariableNames',{'dataset','RMSE_mm3','weighted_SSE','R_squared','N'})];
    parameters=[parameters;table(string(id),p.r,p.K,u0(2),u0(1),u0(3),d.t(1), ...
        'VariableNames',{'dataset','r_per_day','K_mm3','E0_per_day','C0_mm3','Z0','initial_day'})];
    fprintf('%s: RMSE %.8f, R^2 %.8f\n',id,f.rmse,f.r2);
end
save(fullfile(out,'fits.mat'),'bundle');
writetable(metrics,fullfile(out,'fit_metrics.csv'));
writetable(parameters,fullfile(out,'fitted_parameters.csv'));
end

function dydt = tumor_immune_dc_rhs(~,y,p)
%TUMOR_IMMUNE_DC_RHS Unified tumor--T-cell--DC model family.
% State order: y=[C,T,D].  Setting p.sT=0 gives the baseline model;
% p.sT>0 activates antigen-driven T-cell recruitment.

C = y(1);
T = y(2);
D = y(3);

tumorGrowth = p.r*C*(1-C/p.K);
immuneKilling = p.eta*C*T;
antigenRecruitment = 0;
if p.sT ~= 0
    antigenRecruitment = p.sT*C^p.n/(p.KC^p.n+C^p.n);
end
tCellExpansion = p.lambdaT*D/(p.KD+D)*T;
releasedAntigen = C*T;
dendriticActivation = p.lambdaD*releasedAntigen/(p.KA+releasedAntigen);

dydt = [tumorGrowth-immuneKilling; ...
    antigenRecruitment+tCellExpansion-p.dT*T; ...
    dendriticActivation-p.dD*D];
end


function J = tumor_immune_dc_jacobian(y,p)
%TUMOR_IMMUNE_DC_JACOBIAN Analytic Jacobian of the unified model RHS.
% State order: y=[C,T,D] for dimensional parameters or y=[x,y,z] after
% nondimensionalisation. The parameter-field convention is the same as in
% tumor_immune_dc_rhs.m.

C = y(1);
T = y(2);
D = y(3);

hillSlope = p.n*p.KC^p.n*C^(p.n-1)/(p.KC^p.n+C^p.n)^2;
activationSlope = p.lambdaD*p.KA/(p.KA+C*T)^2;

J = [p.r*(1-2*C/p.K)-p.eta*T, -p.eta*C, 0; ...
    p.sT*hillSlope, p.lambdaT*D/(p.KD+D)-p.dT, ...
    p.lambdaT*p.KD*T/(p.KD+D)^2; ...
    activationSlope*T, activationSlope*C, -p.dD];
end
