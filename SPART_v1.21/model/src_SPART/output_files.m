function output_files(Out_dir,params_inputs,output,sensor)
if ~exist(Out_dir, 'dir')
 mkdir(Out_dir)
end

params_inputs   =   params_inputs';
% nsim          =   size(parms_inputs,2);
nparams         =   size(params_inputs,2);
headers     =   {'B','lat','lon','SMp',...
                'Cab','Cdm','Cw','Cs','Cca','Cant','N',...
                'LAI','LIDFa','LIDFb','hot',...
                'Pa','aot550','uo3','alt_m','Pa0',...
                'tts','tto','psi'};

T = array2table(params_inputs);
T.Properties.VariableNames(1:nparams) = headers;                    
writetable(T,[Out_dir,'pars_and_input.dat'],'Delimiter',' ')              
% type 'pars_and_input.dat'  


dlmwrite([Out_dir,'wl_sensor.dat'],sensor.center_wvl);
dlmwrite([Out_dir,'wl_spart.dat'],sensor.wl_smac);
dlmwrite([Out_dir,'TOC_reflectance.dat'],output{1});
dlmwrite([Out_dir,'TOA_reflectance.dat'],output{2});
dlmwrite([Out_dir,'TOA_radiance.dat'],output{3});


% fidvs           =   fopen([Out_dir,'pars_and_input.dat'],'w');
% for j =1:nparams
%     fprintf(fidvs,'%s\t',headers{j});
% end
% fprintf(fidvs,' \r');
% fprintf(fidvs,'%d\t\r',params_inputs);
% fclose(fidvs);
% 
% fclose('all');