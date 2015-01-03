addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/glmnet_matlab' ) );
addpath( genpath( '~/Ipopt_3118' ) );
addpath( genpath( '~/SLEP_package_4.1' ) );
addpath( genpath( '~/minFunc_2012' ) );
load( '~/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat' );
exp_100831_348_136_12_40hr_0_1XLB_LN_all_hypers_grid_search( 120 );
