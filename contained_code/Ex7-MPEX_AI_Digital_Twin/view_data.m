close all
clear all
clc


load('ml_data.mat')

view_exp = 1;% Label to view exp 1-8

slices = find(ml_data_label(:,1)==view_exp);


h_fig = figure;

for j = 1:length(slices)
    cimg = squeeze(ml_data_imag(slices(j),:,:,:));
    imagesc(theta,z,cimg(:,end:-1:1))
    xlabel('radial angle')
    ylabel('z')
    cpar = ml_data_series(slices(j),:);
    str = '';
    for js=1:length(txt_xlsx)
        if js < length(txt_xlsx);
            str = append(str,txt_xlsx{js},':',num2str(round(cpar(js))),' ,');
        else
            str = append(str,txt_xlsx{js},':',num2str(round(cpar(js))),'.');
        end
    end
    str = append('Time: ',num2str(dT(ml_data_label(slices(j),2))),' at r=1.65, control avg parms: ',str);
    title(str)
    drawnow
    pause(.1)

end



