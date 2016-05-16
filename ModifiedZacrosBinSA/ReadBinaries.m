% Script for reading the KMC reweighting SA output files

clear; clc; fclose('all');

format long g

xvec = linspace(1,100,100);
yvec = zeros(1,100);


%% Parmaeters

% get this data from general_output.txt
N_rxns = 6;
N_specs = 2;
n_clusters = 2;

file_read = 'general_output.txt';
fid=fopen(file_read,'r');
eofstat = false;
while ~eofstat
    textLine = fgetl(fid);
    eofstat = feof(fid);
    m1 = strfind(textLine,'Current KMC time:');
    if ~isempty(m1)
        t_end = str2num(textLine(m1+18:end));
    end
    
    m2 = strfind(textLine,'Species number will be reported in file specnum_output.txt every');
    if ~isempty(m2)
        sample_t = str2num(textLine(m2+64:m2+64+24));
    end
    
end

fclose(fid);



%% SA_output.bin

fid=fopen('SA_output.bin','r');
W_mat = fread(fid,'double');
n_samples = length(W_mat)/N_rxns;
t_vec = linspace(0,t_end,n_samples);
W_mat = reshape(W_mat,[N_rxns,n_samples])';

% figure
% hold on
% for i = 1:N_rxns
%     plot(t_vec,W_mat(:,i))
% end
% xlabel('time')
% ylabel('traj deriv')

%% SA_output.bin

fid=fopen('specnum_output.bin','r');
specnum_mat = fread(fid,'int');
specnum_mat = reshape(specnum_mat,[N_specs,n_samples])';

% figure
% hold on
% for i = 1:N_specs
%     plot(t_vec,specnum_mat(:,i))
% end
% xlabel('time')
% ylabel('spec pop')

%% E.bin

fid=fopen('E.bin','r');
E_vec = fread(fid,'double');

% figure
% plot(t_vec,E_vec)
% xlabel('time')
% ylabel('energy (eV)')

%% clusterocc.bin

fid=fopen('clusterocc.bin','r');
cluster_occ_vec = fread(fid,'int');
cluster_occ_vec = reshape(cluster_occ_vec, [n_clusters,n_samples])';

%% Prop_output.bin
fid=fopen('Prop_output.bin','r');
prop_vec = fread(fid,'double');
prop_vec = reshape(prop_vec, [N_rxns,n_samples])';

% prop_vec(:,[1,2]) = prop_vec(:,[1,2]) / 100;
% figure
% hold on
% for i = 1:N_rxns
%     plot(t_vec,prop_vec(:,i))
% end
% xlabel('time')
% ylabel('prop')
% legend('1','2','3','4','5','6')

%% PropCounter_output.bin
fid=fopen('PropCounter_output.bin','r');
propint_vec = fread(fid,'double');
propint_vec = reshape(propint_vec, [N_rxns,n_samples])';

% propint_vec(:,[1,2]) = propint_vec(:,[1,2]) / 100;
% figure
% hold on
% for i = 1:N_rxns
%     plot(t_vec,propint_vec(:,i))
% end
% xlabel('time')
% ylabel('integral prop')
% legend('1','2','3','4','5','6')

%% Hist.bin
fid=fopen('Hist.bin','r');
Hist = fread(fid,'int');
Hist = reshape(Hist,[4,length(Hist)/4])';

% sampling may be different than species numbers
% listed for each site and then onto the next snapshot with no breaks 