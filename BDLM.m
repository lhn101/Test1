%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program : BDLM
% Author : James-A. Goulet
% Date : May 15th 2017
% Description : Process control for Bayesian Dynamic Linear model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%RandStream.setGlobalStream(RandStream('mt19937ar','seed',861040000));  %Initialize random stream number based on clock
set(0,'DefaultAxesFontname','Helvetica')                                %Set default font type
set(0,'DefaultAxesFontSize',20)                                         %Set default font size
format short g                                                          %Set display format
% feature 11
disp(' ')
disp('----------------------------------------------------------------------------------------------')
disp('Bayesian Dynamic Linear Model')
disp('  Version: 2.1 (May 2017)')
disp('  Contact: James.A.Goulet@gmail.com')
disp(' ')
disp(' Issues to be resolved:')
disp('   - state_estimation: m*m''')
disp('   - dynamic regression hidden covariate')
disp('   - file loader')
disp(' ')
disp('----------------------------------------------------------------------------------------------')
disp('////////////////////////////////////Start BDLM')
disp(' ')
disp('----------------------------------------------------------------------------------------------')
disp(' ')
disp('0  -> Create a configuration file for a new dataset')
disp(' ')

%% Load menu displaying existing files
file_data=load_saved_files;
disp(' ')
disp('*     Enter a filename to load a new file')
disp('**    Enter ''delete_%'' to delete a file which number correspond to %')
disp(' ')
file_choice = input('Selection : ');                                                       %Menu 1 - User input (Choice of file)
disp(' ')
if ischar(file_choice)
    if strncmp('delete_',file_choice,7) % Delete a preloaded file
        file_choice=str2mat(file_choice);
        file_choice=str2num(file_choice(1,8:end))+1;
        disp('----------------------------------------------------------------------------------------------')
        disp(['    Deleting data file'])
        disp('----------------------------------------------------------------------------------------------')
        disp(' Are you sure you want to delete the following file?')
        disp(' ')
        disp(['     -> ',file_data{file_choice,1},' /' file_data{file_choice,2}])
        disp(' ')
        choice=input((' Yes(1) or No(0):'));
        if choice==1
            delete(strcat(cd,'/saved_files/',file_data{file_choice,1},'_data.mat'));
            file_data(file_choice,:)=[];
            save(strcat(cd,'/saved_files/','file_data'),'file_data');
            disp(' The file has been deleted')
        else
            disp(' The file was not deleted')
        end
        clear choice
        disp('    ->done.')
        disp(' ')
        disp('//////////////////////////////////////End')
        return
    else % Load a new file
        disp('----------------------------------------------------------------------------------------------')
        disp(['    Load result files '])
        disp('----------------------------------------------------------------------------------------------')
        disp('    ...in progress')
        disp(' ')
        run([cd '/config_files/' file_choice]);                           %Load files
        [model,data,option]=model_initialisation(model,data,option); %Initlize the model
        save(strcat(cd,'/saved_files/',option.name,'_data.mat'));
        i_1=2;
        while i_1<=size(file_data,1)
            if strmatch(option.name,cell2mat(file_data(i_1,1)))
                file_data(i_1,:)=[];
            else
                i_1=i_1+1;
            end
        end
        clear i_1
        file_data(end+1,:)={option.name datestr(now)};
        save(strcat(cd,'/saved_files/','file_data.mat'),'file_data');
        disp(' ')
        disp(['    files saved as: ' cd '/saved_files/' option.name '_data.mat']);
        disp('    ->done.')
        disp(' ')
        disp('//////////////////////////////////////End')
        return
    end
elseif file_choice==0
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp([' / Interactive configuration file builder'])
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    Config_file_creator;
    disp(' ')
    disp('    ->done.')
    disp(' ')
    disp('//////////////////////////////////////End')
    return
else % Load choosen file
    load(strcat(cd,'/saved_files/',file_data{file_choice+1,1},'_data.mat'));
end
disp(' ')
disp('-------------------------------------------------------------------------------')
disp(' / Choose from')
disp('-------------------------------------------------------------------------------')
disp(' ')
disp(' 1  ->  Learn unknown parameters (Newton-Raphson EM)')
disp(' 2  ->  Learn initial values x_0 (SKS)')
disp(' 3  ->  Offline estimation using the Switching Kalman Filter (SKF)')
disp(' 4  ->  Offline estimation using the Switching Kalman Smoother (SKS)')
disp(' ')
disp(' 11 ->  Modify current parameter values')
disp(' 12 ->  Modify current initial x_0 values')
disp(' 13 ->  Modify current training period')
disp(' 14 ->  Plots')
disp(' 15 ->  Display model matrices')
disp(' ')
disp(' ')
user_inputs.inp_1 = input('Selection : ');

tic %ititialize timing variable
if  user_inputs.inp_1==1
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp('1/ Learn BDLM''s unknown parameters with EM')
    disp('-------------------------------------------------------------------------------')
    
    option.iteration_limit_calibration=100;
    option.time_limit_calibration=60;%[min]
    
    [optim]=NR_EM(data,model,option);
    
    model.parameter(model.p_ref)=optim.parameter_opt(model.p_ref);
    
    
    save(strcat(cd,'/saved_files/',option.name,'_data.mat'));
    disp(['    files saved as: ' cd '/saved_files/' option.name '_data.mat']);
    disp('    ->done.')
    
elseif user_inputs.inp_1==2
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp('2/ Learn BDLM''s initial values x_0')
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    estim=state_estimation(data,model,option,'smooth',1,'disp_flag',0);
    for i=1:model.nb_class
        model.initX{i}=estim.x_M{i}(:,1);
        model.initV{i}=estim.V_M{i}(:,:,1);
        model.initS{i}=estim.S(1,i);
    end
    save(strcat(cd,'/saved_files/',option.name,'_data.mat'));
    disp(['    files saved as: ' cd '/saved_files/' option.name '_data.mat']);
    disp('    ->done.')
    
elseif user_inputs.inp_1==3
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp('3/ Offline estimation of x_t using the Switching Kalman Filter (SKF)')
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    estim=state_estimation(data,model,option,'smooth',0);
    plot_estimations(estim,data,model,option)
    
elseif user_inputs.inp_1==4
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp('4/ Offline estimation of x_t using the Switching Kalman Smoother (SKS)')
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    estim=state_estimation(data,model,option,'smooth',1);
    plot_estimations(estim,data,model,option)
    
elseif user_inputs.inp_1==11
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp([num2str(user_inputs.inp_1) '/ Modify current parameter values'])
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    disp('#   |Parameter    |Component |Model # |Observation  |Current value |Bounds min/max |Constraint')
    for i=1:length(model.parameter)
        i_ref=model.p_ref(i);
        if i~=i_ref
            contraint=['@' num2str(i_ref)];
        else
            contraint='';
        end
        disp([repmat('0',1,2-length(num2str(i))) num2str(i) '   ' model.param_properties{i,1} repmat(' ',1,14-length(model.param_properties{i,1})) model.param_properties{i,2} repmat(' ',1,11-length(model.param_properties{i,2})) model.param_properties{i,3} repmat(' ',1,9-length(model.param_properties{i,3}))...
            model.param_properties{i,4} repmat(' ',1,14-length(model.param_properties{i,4})) num2str(model.parameter(i_ref)) repmat(' ',1,15-length(num2str(model.parameter(i_ref))))...
            num2str(model.param_properties{i,5}(1)) repmat(' ',1,5-length(num2str(model.param_properties{i,5}(1)))) '/' num2str(model.param_properties{i,5}(2)), repmat(' ',1,9-length(model.param_properties{i,5})) contraint ]);
    end
    disp(' ')
    disp(' 1   ->  Modify a parameter value')
    disp(' 2   ->  Constrain a parameter to another')
    disp(' 3   ->  Export current parameter properties in config file format')
    disp(' ')
    user_inputs.inp_2 =  input('Selection : ');
    if user_inputs.inp_2==1
        user_inputs.inp_2 =  input(' Modify parameter # ');
        user_inputs.inp_3 = input('        New value : ');
        user_inputs.inp_4 = input('       New bounds : ');
        if ~isempty(user_inputs.inp_3)
            model.parameter(user_inputs.inp_2)=user_inputs.inp_3;
        end
        if ~isempty(user_inputs.inp_4)
            model.param_properties{user_inputs.inp_2,5}=user_inputs.inp_4;
        end
        disp(' ')
        clear i space
        save(strcat(cd,'/saved_files/',option.name,'_data.mat'));
        disp(['    files saved as: ' cd '/saved_files/' option.name '_data.mat']);
    elseif user_inputs.inp_2==2
        user_inputs.inp_2 = input('  Constrain parameter # ');
        user_inputs.inp_3 = input('         to parameter # ');
        model.p_ref(user_inputs.inp_2)=user_inputs.inp_3;
        model.param_properties{user_inputs.inp_2,5}=[nan,nan];
        disp(' ')
        clear i space
        save(strcat(cd,'/saved_files/',option.name,'_data.mat'));
        disp(['    files saved as: ' cd '/saved_files/' option.name '_data.mat']);
    elseif user_inputs.inp_2==3
        disp(' ')
        disp('model.param_properties={')
        for i=1:size(model.param_properties,1)
            space=repmat(' ',1,8-length(model.param_properties{i,1}));
            disp(sprintf(['\t''%-s''' space ',\t ''%-s'',\t ''%-s'',\t ''%-s'',\t [ %-5G, %-5G]\t %%#%d'], model.param_properties{i,:},i));
        end
        disp('};')
        
        disp(' ')
        disp('model.parameter=[')
        for i=1:size(model.parameter,1)
            disp(sprintf('%-8.5G \t %%#%d', model.parameter(i),i));
        end
        disp(']; ')
        
        disp(' ')
        disp(['model.p_ref=[' num2str(model.p_ref) '];'])
        disp(' ')
    else
        disp('   !!! -> invalid input...')
    end
    disp(' ')
    disp('    ->done.')
    
elseif user_inputs.inp_1==12
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp([num2str(user_inputs.inp_1) '/ Modify current initial x_0 values'])
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    disp('#  |state variable    |observation             |E[x_0]         |var[x_0] ')
    for i=1:length(model.initX{1})
        disp([repmat('0',1,2-length(num2str(i))) num2str(i) '  ' model.hidden_states_names{1}{i,1} repmat(' ',1,19-length(model.hidden_states_names{1}{i,1})) model.hidden_states_names{1}{i,3} repmat(' ',1,25-length(model.hidden_states_names{1}{i,3})) sprintf('%-10.5G',model.initX{1}(i)) repmat(' ',1,6) sprintf('%-10.5G',model.initV{1}(i,i))])
    end
    disp(' ')
    disp(' 1   ->  Modify a initial values')
    disp(' 2   ->  Export initial values in config file format')
    disp(' ')
    user_inputs.inp_2 =  input('Selection : ');
    if user_inputs.inp_2==1
        user_inputs.inp_3 =     input('   Modify variable # ');
        user_inputs.inp_4 = input('            New E[x_0] : ');
        user_inputs.inp_5 = input('          New var[x_0] : ');
        for i=1:model.nb_class
            model.initX{i}(user_inputs.inp_3)=user_inputs.inp_4;
            D=diag(model.initV{1});
            D(user_inputs.inp_3)=user_inputs.inp_5;
            model.initV{i}=diag(D);
        end
        clear D i
        save(strcat(cd,'/saved_files/',option.name,'_data.mat'));
        disp(['    files saved as: ' cd '/saved_files/' option.name '_data.mat']);
    elseif user_inputs.inp_2==2
        for m=1:model.nb_class
            disp(' ')
            disp(['model.initX{' num2str(m) '}=['])
            for i=1:size(model.initX{m},1)
                disp(sprintf(['\t%-6.3G'], model.initX{m}(i,:)));
            end
            disp('];')
            
            disp(' ')
            disp(['model.initV{' num2str(m) '}=['])
            for i=1:size(model.initV{m},1)
                disp(sprintf(['\t%-8.3G'], model.initV{m}(i,:)));
            end
            disp('];')
            
            disp(' ')
            for i=1:size(model.initS{m},1)
                disp(sprintf('model.initS{%d}=[%-6.3G];', m, model.initS{m}));
            end
        end
        
    else
        disp('   !!! -> invalid input...')
    end
    disp(' ')
    disp('    ->done.')
    
elseif user_inputs.inp_1==13
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp([num2str(user_inputs.inp_1) '/ Modify current training period'])
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    disp(['   Current training period: from:' num2str(option.training_period(1)) ' to:' num2str(option.training_period(2)) ' days'])
    disp(' ')
    option.training_period(1) =     input('   New start point [days]: ');
    option.training_period(2) =     input('     New end point [days]: ');
    disp(' ')
    option.training_start_idx=find(abs(data.timestamps-option.training_period(1))==min(abs(data.timestamps-option.training_period(1))),1,'first');
    option.training_end_idx=find(abs(data.timestamps-option.training_period(2))==min(abs(data.timestamps-option.training_period(2))),1,'first');
    
    save(strcat(cd,'/saved_files/',option.name,'_data.mat'));
    disp(['    files saved as: ' cd '/saved_files/' option.name '_data.mat']);
    disp('    ->done.')
    
elseif user_inputs.inp_1==14
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp([num2str(user_inputs.inp_1) '/ Plot data'])
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    disp(' 1  ->  Plot raw datasets')
    disp(' 2  ->  Plot missing data')
    disp(' 3  ->  Plot time steps size')
    for i=1:size(model.param_properties,1)
        if strcmp(model.param_properties{i,2},'DH')
            msg=[];
            break
        else msg='(No dynamic regression hidden components (DH))';
        end
        
    end
    disp([' 4  ->  Plot DH''s hidden covariate ' msg])
    disp(' [] ->  No plot')
    
    disp(' ')
    user_inputs.inp_2 = input('Selection : ');
    if user_inputs.inp_2==1
        plot_data(data,option)
    elseif user_inputs.inp_2==2
        plot_missing_data(data,option)
    elseif user_inputs.inp_2==3
        plot_time_steps(data,option)
    elseif user_inputs.inp_2==4
        plot_DRHC(model,option)
    end
    
elseif user_inputs.inp_1==15
    disp(' ')
    disp('-------------------------------------------------------------------------------')
    disp([num2str(user_inputs.inp_1) '/ Display model matrices'])
    disp('-------------------------------------------------------------------------------')
    disp(' ')
    user_inputs.inp_2 = input('  Choose a timestasmp index : ');
    
    
    for m=1:model.nb_class
        disp(' ')
        disp('------------------')
        disp([' Model class #' num2str(m)])
        disp('------------------')
        disp(' ')
        name1=[];
        for i=1:size(model.hidden_states_names{1},1)
            name1=[name1,{model.hidden_states_names{m}{i,2},model.hidden_states_names{m}{i,3}(1)}];
        end
        disp(['       ' sprintf('M%s|%s      ',name1{:})])
        disp(['   ' sprintf('%10s',model.hidden_states_names{1}{:,1})])
        M=model.A{m}(model.parameter(model.p_ref),data.dt_ref,user_inputs.inp_2);
        for i=1:size(M,1)
            if i==1
                disp(['    A=[' sprintf('%-10.3G',M(i,:))]);
            elseif i==size(M,1)
                disp(['       ' sprintf('%-10.3G',M(i,:)) '];']);
                
            else
                disp(['       ' sprintf('%-10.3G',M(i,:))]);
            end
        end
        disp(' ')
        M=model.C{m}(model.parameter(model.p_ref),data.dt_ref,user_inputs.inp_2);
        for i=1:size(M,1)
            if i==1
                disp(['    C=[' sprintf('%-10.2G',M(i,:))]);
            elseif i==size(M,1)
                disp(['       ' sprintf('%-10.2G',M(i,:)) '];']);
                
            else
                disp(['       ' sprintf('%-10.2G',M(i,:))]);
            end
        end
        disp(' ')
        M=model.Q{m}{1}(model.parameter(model.p_ref),data.dt_ref,user_inputs.inp_2);
        for i=1:size(M,1)
            if i==1
                disp([' Q_' num2str(m) num2str(m) '=[' sprintf('%-10.2G',M(i,:))]);
            elseif i==size(M,1)
                disp(['       ' sprintf('%-10.2G',M(i,:)) '];']);
                
            else
                disp(['       ' sprintf('%-10.2G',M(i,:))]);
            end
        end
        if model.nb_class>1
            for j=setdiff(1:size(M,1),m)
                disp(' ')
                M=model.Q{m}{j}(model.parameter(model.p_ref),data.dt_ref,user_inputs.inp_2);
                for i=1:size(M,1)
                    if i==1
                        disp([' Q_' num2str(m) num2str(j) '=[' sprintf('%-10.2G',M(i,:))]);
                    elseif i==size(M,1)
                        disp(['       ' sprintf('%-10.2G',M(i,:)) '];']);
                        
                    else
                        disp(['       ' sprintf('%-10.2G',M(i,:))]);
                    end
                end
            end
        end
    end
    disp(' ')
    disp(['     ' sprintf('%11s  ',data.labels{:})])
    M=model.R{m}(model.parameter(model.p_ref),data.dt_ref,user_inputs.inp_2);
    for i=1:size(M,1)
        if i==1
            disp(['    R=[' sprintf('%-11.2G',M(i,:))]);
        elseif i==size(M,1)
            disp(['       ' sprintf('%-11.2G',M(i,:)) '];']);
            
        else
            disp(['       ' sprintf('%-11.2G',M(i,:))]);
        end
    end
    
    disp(' ')
    disp('    ->done.')
end
time_elapsed=toc;
disp(' ')
disp(['/// Time elapsed: ' num2str(time_elapsed/60) ' min'])
disp('//////////////////////////////////////End')