load('D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN\surfactin_RN.mat');
singleData.mzAxis = mzAxis;
singleData.dataCube = dataCube;
figure;
plot(singleData.dataCube);
text(975,2540+75,'\downarrow','fontsize',20)
text(975,2540+75+150,'m/z:1021.8')
text(977-25,986+75,'\downarrow','fontsize',20)
text(977-25-175, 986+75+150, 'm/z:1007.7' );
text(1021-25,2200+75,'\downarrow','fontsize',20);
text(1021-25, 2200+75+150, 'm/z: 1035.8' );
xlabel('# signals');
ylabel('values');
export_fig surfactin.pdf -transparent

load('C:\Users\¦Ð±á\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat');
load('C:\Users\¦Ð±á\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LN\exp_100831_348_136_12_40hr_0_1XLB_LN_7_res.mat');