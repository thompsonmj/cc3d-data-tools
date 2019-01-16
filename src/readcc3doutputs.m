% 
% % sourceFile = ['Y:\multi_simple_liver\' ...
% %     '803108.brown-adm.rcac.purdue.edu\sim_000001\output\screenshots\' ...
% %     'data_log.txt'];
% 
% simDir = ['Y:\multi_simple_liver\803108.brown-adm.rcac.purdue.edu'];
% here = pwd;
% 
% dataFile = ['\data_log.txt'];
% 
% nSims = size( dir(simDir), 1 ) - 5; % 5 files not sim_xxxxxx dirs.
% 
% for iSim = 1:nSims
%     simNo = sprintf('%06d',iSim);
%     sim_id = ['sim_',simNo];
% %     sim_idArray{iSim} = sim_id;
% end
% 
%% How to create a struct with the output data from Jim's sim
% for iSim = 1:nSims
%     simNo = sprintf('%06d',iSim);
%     sim_id = ['sim_',simNo];
%     dataFile = ['\data_log.txt'];
%     dataPath = [simDir,'\',sim_id,'\output\screenshots',dataFile];
%     copyfile(dataPath,here)
%     dataFile = ['data_log.txt'];
%     d = importdata(dataFile,'\t');
%     data(1).(sim_id) = d;
%     fclose all;
%     delete 'data_log.txt'
%     disp(iSim)
% end

%% How to read the struct called data.
% Data columns:
% data.sim_xxxxxx.data(:,colNo);
% mcs:                              col1
% time(hr):                         col2
% PBPKspeciesDictAPAP['CVen']:      col12
% PBPKspeciesDictAPAPG['CVen']:     col13
% PBPKspeciesDictAPAPS['CVen']:     col14


%% How to write the data I want to csv
% for iSim = 1:nSims
%     simNo = sprintf('%06d',iSim);
%     sim_id = ['sim_',simNo];
%     mcs = data(1).(sim_id).data(:,1);
%     time_h = data(1).(sim_id).data(:,2);
%     APAP = data(1).(sim_id).data(:,12);
%     APAPS = data(1).(sim_id).data(:,13);
%     APAPG = data(1).(sim_id).data(:,14);
%     T = table(mcs, time_h, APAP, APAPS, APAPG);
%     filename = ['Z:\cc3d\apps\CC3D_v377\_matlab\multi_simple_liver\',sim_id,'.csv'];
%     writetable(T,filename)
%     disp(iSim)
% end

%% Plot the dynamics
dataDir = ['..\_models\multi_simple_liver\output'];

for iSim = 1:nSims
    simNo = sprintf('%06d',iSim);
    sim_id = ['sim_',simNo];
    
    
    
    mcs = data(1).(sim_id).data(:,1);
    time_h = data(1).(sim_id).data(:,2);
    APAP = data(1).(sim_id).data(:,12);
    APAPS = data(1).(sim_id).data(:,13);
    APAPG = data(1).(sim_id).data(:,14);
    T = table(mcs, time_h, APAP, APAPS, APAPG);
    filename = ['Z:\cc3d\apps\CC3D_v377\_matlab\multi_simple_liver\',sim_id,'.csv'];
    writetable(T,filename)
    disp(iSim)
end
