clearvars
clc
close all
do_vis = 0;

cpath = pwd;
path_CoilSetup = strcat(cpath,'/CoilSetup_ProtoMPEX.xlsx');
path_fuction = strcat(cpath,'/Functions');
path_MPEX_DT = strcat(cpath,'/Ex7-MPEX_AI_Digital_Twin');
path_MPEX_Obs = strcat(cpath,'/ml_data.mat');
path_surrogate_obs_data = strcat(cpath,'/ml_surrogate.mat');

load(path_surrogate_obs_data)
load(path_MPEX_Obs)
p = addpath(genpath(path_fuction));
p = addpath(genpath(path_MPEX_DT));


par_range = [[530-100 530+100];[2300-500 2300+500];[3833 6800];[3833 6800];[160 430]];
n_samples = 1000;
parms=lhsu(par_range(:,1),par_range(:,2),n_samples);
parms = round(parms);
parms = [[530   2300    4000 4000    430]; parms];
parms = [[530   2300    6800 4000    160]; parms];
cparms_tune = [1.2871e+19 5.6259e+17 5.5910 3.5513];

us_ind = find(ml_data_series(:,4)==parms(1,3) & ml_data_series(:,6)==parms(1,5));
figure
sgtitle(strcat('Experiment T1=',num2str(parms(1,1)),' T2=',num2str(parms(1,2)),' P1=',num2str(parms(1,3)),' P2=',num2str(parms(1,4)),' P3=',num2str(parms(1,5))))
imagesc(theta, z(2:end),1000*squeeze(mean(ml_data_imag(us_ind,1,2:end,:))))
colormap("hot")

us_ind = find(ml_data_series(:,4)==parms(2,3) & ml_data_series(:,6)==parms(2,5));
figure
sgtitle(strcat('Experiment T1=',num2str(parms(2,1)),' T2=',num2str(parms(2,2)),' P1=',num2str(parms(2,3)),' P2=',num2str(parms(2,4)),' P3=',num2str(parms(2,5))))
imagesc(theta, z(2:end),1000*squeeze(mean(ml_data_imag(us_ind,1,2:end,:))))
colormap("hot")


h_fig = figure;
for jp =1:size(par_all,1)

    [tmpv,us_ind] = min(abs(par_all(:,1)-parms(jp,3))+abs(par_all(:,2)-parms(jp,5)));
    parms_c = parms(jp,:);
    parms_c([3,5]) = par_all(us_ind,:);
    [r_h,z_h,r_t,z_t,Q_tot_h,Q_tot_t,r_bz0,z_bz0,bz0,Qtot_map] = FieldOracle(parms_c,path_CoilSetup,do_vis,h_fig,cparms_tune);
    
    closs = 100*sum(Q_tot_h)/sum(Q_tot_t)+mean((Q_tot_h-mean(Q_tot_h)).^2)/mean(Q_tot_h)^2;
    if jp == 1
        figure(h_fig)
        sgtitle(strcat('Surrogate (loss=',num2str(closs),') T1=',num2str(parms(jp,1)),' T2=',num2str(parms(jp,2)),' P1=',num2str(parms(jp,3)),' P2=',num2str(parms(jp,4)),' P3=',num2str(parms(jp,5))))
        subplot(1,2,1)
        plot(z_h,Q_tot_h(:),'-b',z_h(2:end),1000*mean(img_all(us_ind,:,:),3)'+Q_tot_h(2:end)','-g','LineWidth',2)
        subplot(1,2,2)
        imagesc(theta, z_h(2:end),1000*squeeze(img_all(us_ind,:,:))+repmat(Q_tot_h(2:end)',1,100))
        colormap("hot")
        drawnow
        bloss = closs;
        bind = jp;
    elseif jp == 2
        figure(h_fig)
        sgtitle(strcat('Surrogate (loss=',num2str(closs),') T1=',num2str(parms(jp,1)),' T2=',num2str(parms(jp,2)),' P1=',num2str(parms(jp,3)),' P2=',num2str(parms(jp,4)),' P3=',num2str(parms(jp,5))))
        subplot(1,2,1)
        plot(z_h,Q_tot_h(:),'-b',z_h(2:end),1000*mean(img_all(us_ind,:,:),3)'+Q_tot_h(2:end)','-g','LineWidth',2)
        subplot(1,2,2)
        imagesc(theta, z_h(2:end),1000*squeeze(img_all(us_ind,:,:))+repmat(Q_tot_h(2:end)',1,100))
        colormap("hot")
        drawnow
        if closs<bloss
            bloss = closs;
            bind = jp;
        end
    elseif closs<bloss
        bloss = closs;
        bind = jp;
        figure(h_fig)
        sgtitle(strcat('Surrogate (loss=',num2str(closs),') T1=',num2str(parms(jp,1)),' T2=',num2str(parms(jp,2)),' P1=',num2str(parms(jp,3)),' P2=',num2str(parms(jp,4)),' P3=',num2str(parms(jp,5))))
        subplot(1,2,1)
        plot(z_h,Q_tot_h(:),'-b',z_h(2:end),1000*mean(img_all(us_ind,:,:),3)'+Q_tot_h(2:end)','-g','LineWidth',2)
        subplot(1,2,2)
        imagesc(theta, z_h(2:end),1000*squeeze(img_all(us_ind,:,:))+repmat(Q_tot_h(2:end)',1,100))
        colormap("hot")
        drawnow
    end
    
end

