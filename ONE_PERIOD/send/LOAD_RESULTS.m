clear; close all; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex')

%% load parameters

cdir = pwd;

cd('OUTPUT_1');

load bvec.txt;
load zvec.txt;
load gvec.txt;
load bpol.txt;
load Vnd.txt;
load Vd.txt;
load Vndk.txt;
load Wv.txt;
load Wtv.txt;
load dpol.txt;
load epol.txt;
load Ednp.txt;
load Rsched.txt;
cd(cdir)

Nb = numel(bvec); Nz = numel(zvec); Ng = numel(gvec);

%% reshape

%--plicies and values
bpolxl = bpol; bpol = 189*ones(Nb,Ng,Nz);
Vndxl  = Vnd ; Vnd  = 189*ones(Nb,Ng,Nz);
Vdxl   = Vd  ; Vd   = 189*ones(Nb,Ng,Nz);
Vndkxl = Vnd ; Vndk = 189*ones(Nb,Ng,Nz);
Wvxl   = Wv  ; Wv   = 189*ones(Nb,Ng,Nz);
Wtvxl  = Wtv ; Wtv  = 189*ones(Nb,Ng,Nz);

inn = 0;
for ig = 1:Ng
for iz = 1:Nz
    for ib = 1:Nb
        inn = inn+1;

        %--policies
        bpol(ib,ig,iz) = bpolxl(inn);

        %--values
        Vnd(ib,ig,iz)  = Vndxl(inn);
        Vd(ib,ig,iz)   = Vdxl(inn);
        Vndk(ib,ig,iz) = Vndkxl(inn);
        Wv(ib,ig,iz)   = Wvxl(inn);
        Wtv(ib,ig,iz)  = Wtvxl(inn);
    end
end
end

%--schedule
Rschedxl = Rsched; Rsched = 189*ones(Nb,Ng);
nsched = 189*ones(Nb,Ng);
inn = 0;
for ig = 1:Ng
    for ib_np = 1:Nb
        inn = inn+1;
        Rsched(ib_np,ig) = Rschedxl(inn);
    end

    gg = gvec(ig);
    nsched(:,ig) = (gg.*bvec)./Rsched(:,ig);  % here: bvec = b_np
end

rhosched = Rsched-1;
%% computations

%% Plots: values and policies

ccxg = [0.8 0.4 0.4;...
      0.4 0.4 0.8];

iz=9; ig=1;
figure(101)
plot(bvec,Vnd(:,ig,iz),'Color',ccxg(ig,:),'LineWidth',3,'LineStyle','-'); hold on
plot(bvec,Vd(:,ig,iz) ,'Color',ccxg(ig,:),'LineWidth',3,'LineStyle','--')
set(gca,'Xgrid','on','Ygrid','on','LineWidth',0.50,'Fontsize',21)
title('Values','Interpreter','LaTex','Fontsize',27)
xlabel('$b$ ','Interpreter','LaTex','Fontsize',27)
leg = legend('$v^{nd}$','$v^d$');
set(leg,'Interpreter','LaTex','Fontsize',23)
legend boxoff
hold off

figure(102)
ig=1;plot(bvec,bpol(:,ig,iz),'Color',ccxg(ig,:),'LineWidth',3,'LineStyle','-'); hold on
ig=2;plot(bvec,bpol(:,ig,iz),'Color',ccxg(ig,:),'LineWidth',3,'LineStyle','-'); 
%plot(bvec,bvec,'Color',[0.1 0.1 0.1],'LineWidth',2,'LineStyle','--'); 
set(gca,'Xgrid','on','Ygrid','on','LineWidth',0.50,'Fontsize',21)
title('$b^{\prime}(\cdot)$','Interpreter','LaTex','Fontsize',27)
xlabel('$b$ ','Interpreter','LaTex','Fontsize',27)
leg = legend('$g_L$','$g_H$');
set(leg,'Interpreter','LaTex','Fontsize',23)
legend boxoff
hold off

%% lots: schedule

figure(1001)
ig=1;plot(nsched(:,ig),rhosched(:,ig),'Color',ccxg(ig,:),'LineWidth',3); hold on
ig=2;plot(nsched(:,ig),rhosched(:,ig),'Color',ccxg(ig,:),'LineWidth',3); 
set(gca,'Xgrid','on','Ygrid','on','LineWidth',0.50,'Fontsize',21)
title('Interest rate schedule','Interpreter','LaTex','Fontsize',27)
xlabel('issuance $n$ ','Interpreter','LaTex','Fontsize',27)
ylabel('rate $\rho$ ','Interpreter','LaTex','Fontsize',27)
leg = legend('$g_L$','$g_H$');
set(leg,'Interpreter','LaTex','Fontsize',23)
legend boxoff
hold off