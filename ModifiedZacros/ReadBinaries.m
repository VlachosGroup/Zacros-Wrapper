% Script for reading the KMC reweighting SA output files

function ReweightSAbin

clear; clc; fclose('all');

format long g

xvec = linspace(1,100,100);
yvec = zeros(1,100);


%% Parmaeters
cutoff = 0; %for when you do sample on specnum every 1.00e-6 seconds
%kB = 8.617332478e-5;               % eV/K
kB = 1.380648813e-23;               % J/K


N_clusters = 6;                     % This needs to be automated for an arbitrary number of reactions
N_rxns = 2;
site_norm = 1;                       % Doesn't matter because we normalize the sensitivity coefficients anyway

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

sample_points = floor(t_end/sample_t)+1;

fclose(fid);

%% clusterocc.bin
% We do this to get the number of each ECI at each timestep

file_read = 'clusterocc.bin';
fid=fopen(file_read,'r');

W = fread(fid,'int');
W = reshape(W,[N_clusters,length(W)/N_clusters])';
W(1:cutoff,:) = [];           % Cut off transient data
fclose(fid);

%% E.bin
% We do this to get the energy at each timestep

file_read = 'E.bin';
fid=fopen(file_read,'r');
E = fread(fid,'double');    

%% Hist.bin
% We do this to get the lattice state at each timestep

file_read = 'Hist.bin';
fid=fopen(file_read,'r');
hist = fread(fid,'int');
hist = reshape(hist,[4,length(hist)/4])';


%% Prop_output.bin
% We do this to get the sum of propensities for each reaction type at each timestep

file_read = 'Prop_output.bin';
fid=fopen(file_read,'r');

props = fread(fid,'double');
props = reshape(props,[N_rxns,length(props)/N_rxns])';

end
