[filename,path] = uigetfile('C:\Users\seifer\PycharmProjects\Ocean\QA\framework*.mat','Fetch frameowrk.mat file with node mapping');
filen=[path filename];
load(filen);

bW=input('Board dimension: '); %board width
bS=bW*bW;%size of board
b4S=bS*4;

bestspin=1*ones(size(NodeSpin));
indg=1:length(bestspin);

count_bad=0;
for ind_location=1:bS
    ind_option=4*(ind_location-1)+1; %known for all virtual puzzles
    ind_spin=indg(mapnodes_loc==ind_location & mapnodes_ind==ind_option);
    if isempty(ind_spin)
        count_bad=count_bad+1;
    else
        bestspin(min(ind_spin))=-1;
    end
    ancil_ind_spin=indg(mapnodes_loc==ind_location & mapnodes_ind==-1);
    bestspin(min(ancil_ind_spin))=-1;
end
disp(sprintf('Bad count=%g ',count_bad));
filenamebest='C:\Users\seifer\PycharmProjects\Ocean\QA\bestspin.csv';
writematrix(bestspin,filenamebest);
return;

% Add these lines to existing python file instead of the analysis part (after bqm settings):
filename='C:/Users/seifer/PycharmProjects/Ocean/QA/d_wave_answer7239.csv'
spin_array = np.loadtxt(filename, dtype=int, delimiter=',')
spin_dict = dict(zip(bqm.variables, spin_array))
energy = bqm.energy(spin_dict)
print("Hamiltonian value (energy):", energy)
filename2='C:/Users/seifer/PycharmProjects/Ocean/QA/bestspin.csv'
spin_array2 = np.loadtxt(filename2, dtype=int, delimiter=',')
spin_dict2 = dict(zip(bqm.variables, spin_array2))
energy2 = bqm.energy(spin_dict2)
print("Hamiltonian value (energy) known solution :", energy2)
