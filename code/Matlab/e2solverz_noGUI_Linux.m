function e2solverz_noGUI_Linux()
global screenCoord;
global CH;
global e2folder;
global PaidDwaveService;
global bW; %board width
global bS;%size of board
global b4S;
global OceanFolder;
OceanFolder='/gpfs01/home/seifer/eternity30/Ocean';

bW=16;
bS=bW*bW;%size of board
b4S=bS*4;
rng('shuffle'); % Seeds the random number generator based on the current time
%seed = sum(clock);
%rng(seed, 'twister'); % Set the random seed

compromize=false;
extraline=0;
flagloop=false;
e2folder=input('Eternity folder: ','s');
screenCoord=[0 0];
parallelForce=input('Number of CPU cores? ');%19;
clue_string=input('enter 4 digits, 1 for each clue included (example2: 1111, 1010):  ','s');
disp('Eternity II solver written by Shahar Seifer (C) 2024-2025');
disp('Enter operation key, or enter ! and then type the tile coordinates as asked and then insert the operation key.');
disp('c- remove, v- match tiles, f- fix one selection per shown tile, show option counts, ');
disp('t- enter tile as 4-letter text, z- show a list of all open choices, space- fit all general choices, l- load board state file');
disp('q- solver in discrete approach, p- solver in probablitsic approach, x- match tile according to best score by discrete solver');
disp('k- solver base on search replacements of 4 tiles, r- fill out tiles randomly, o- toggle compromize ON/OFF, j-toggle between probablistic and discrete');
disp('h- forms solved frame block that matches pieces already there,  y-cut out all but 1st outer line frame, u- cut all but 1-2 lines, m- machine solver after load')
disp('s- save board state,  e- export txt file compatible with eternityII editor, esc- end program.');
disp('/ - toggle for assistance from dwave while solving in h option (frame block). 1-9,0,= update board with heavey caclulation from dwave. ');
disp('* - turn no/off payed D-Wave quantum computation instead of simulation')
disp('> - consult AI')
disp('Ocean and Graph2seq environment variables should be defined by the user in Windows according to the folders used in python.')
disp(' ');
PaidDwaveService=false;
disp('Quantum annealing calculation will use simulation on the local computer for free');
disp(' ')
original=[];
planv=1;
planv=input('0- native file virtual puzzle,  1-native file commercial puzzle,  2-load text from eternityII program, 3-commercial puzzle without mandatory center piece:  ');
%load('/home/labs/elbaum/Collaboration/Matlab/shahar/origin6.mat','-mat'); %load originial cell array
if planv>0
  load(sprintf('%s//origin6.mat',e2folder),'-mat'); %load originial cell array
elseif planv==0
  load(sprintf('%s//origin_virtual.mat',e2folder),'-mat'); %load originial cell array
end

keep_options=[];
machine_state=0;
slash_mode=false;

isprob=input('Discrete (0) or probablistic (1) ?');


if planv==0 || planv==1 || planv==3
    filename=input('File name in folder: ','s');%[filename,path] = uigetfile(sprintf('%s//*.mat',e2folder),'Fetch mat file after analysis');
    path=[e2folder '/'];
    if sum(filename>0)>0
        filen=[path filename];
        disp(['loading ' filen]);
        best_options=[];
        ghost_options=[];
        load(filen,'-mat');
        if ~isempty(best_options)
            options=best_options;
        end
        if ~isempty(ghost_options)
            options=ghost_options;
        end
        keep_options=options;
        if isempty(isprob)
            isprob=0;
        end
    else
        if planv==0
            flaggen=input('0- empty board,  1- generate full board,  2- generate frame ? ');
            if flaggen==1
                 options=zeros(256,1024);
                  for ind_location=1:256
                      options(ind_location,4*(ind_location-1)+1)=1;
                  end
                  keep_options=options;
            elseif flaggen==0
                %keep_options=zeros(256,1024);
            elseif flaggen==2
                 options=zeros(256,1024);
                  for ind_location=1:16
                      options(ind_location,4*(ind_location-1)+1)=1;
                  end
                  for ind_location=256-15:256
                      options(ind_location,4*(ind_location-1)+1)=1;
                  end
                  for ind_location=16+1:16:256-32+1
                      options(ind_location,4*(ind_location-1)+1)=1;
                  end
                  for ind_location=32:16:256-16
                      options(ind_location,4*(ind_location-1)+1)=1;
                  end
                  keep_options=options;
            end
        else
             %keep_options=zeros(256,1024);
        end
    end
elseif planv==2
    filename=input('File name in folder (txt from EternityII editor): ','s');
    path=[e2folder '/'];
    %isprob=1; %to remove non-matched pieces
    if sum(filename>0)>0
        filen=[path filename];
        lines=readlines(filen);
        linesm=erase(erase(lines,'['),']');
        writelines(linesm,filen);
        op=readmatrix(filen);
        options=external2native(op,original);
        keep_options=options;
    else
       keep_options=zeros(256,1024);
    end
    
end


%load oringal{r:1-16,c=1-16}= [sym_up, sym_right, sym_down, sym_left]
piece=zeros(16*16,4);
for row=1:16
    for col=1:16
        ob=original{row,col};
        ind=(col-1)+(row-1)*16+1;
        piece(ind,1:4)=ob;
    end
end



%SHOWING SESSION
%neworiginal=fill_neworiginal(options,piece);
%show_original(neworiginal);

%center r=9,c=8:  18,6,6,11
%clue1 r=14,c=3:  16,6,21,11
%clue2 r=3,c=14:  22,20,22,16
%clue3 r=14,c=14:  4,21,15,12
%clue4 r=3,c=3:    22,3,11,19


%build tree:  from obj1(piece1,rot1) via leg (1:up, 2:right, 3:down,4:left) to obj2(piece2,rot2)
% if edge symbol write to zero object
% The rotations are all anticlockwise !
%Initial state, all connections included. [piece2]=floor((obj2-1)/4)+1, [rot2]=(obj2-1)%4
tree=-1*ones(256*4,4,57);  % from which object according to position (1 to 256*4) to all possible objects(1 of 256*4) up to 128 (actually there are 49 max),  according to conntection direction (1 of 4).
treesize=zeros(256*4,4);
for direct1=1:4
    for piece1=1:256
        for rot1=1:4   %rotation(1,2,3,4)=(no rot, 1 CW, 2 CW, 3 CW),  direction(1,2,3,4)=up, right,down, left
            
            treeind1=(piece1-1)*4+rot1;
            sym1=piece(piece1,mod(direct1-1-(rot1-1),4)+1);
            if sym1==0
                tree(treeind1,direct1,1)=0;
                treesize(treeind1,direct1)=1;
                continue;
            end
            for piece2=1:256
                if piece2==piece1
                    continue;
                end
                for rot2=1:4
                    sym2=piece(piece2,mod(direct1+2-1-(rot2-1),4)+1); %must match with the reverse direction of direct1, so direct2=direct1+2
                    if sym2==sym1
                        treeind2=(piece2-1)*4+rot2;
                        treesize(treeind1,direct1)=treesize(treeind1,direct1)+1;
                        tree(treeind1,direct1,treesize(treeind1,direct1))=treeind2;
                    end
                end
            end
        end
    end
end

%load('C:\Users\seifer\Documents\Matlab\eternity\tree3d.mat','-mat'); %overwrite tree with more narrow network generated in selfware3.m


%clue1 r=14,c=3:  16,6,21,11
%clue2 r=3,c=14:  22,20,22,16
%clue3 r=14,c=14:  4,21,15,12
%clue4 r=3,c=3:    22,3,11,19

clue1=[16 6 21 11];
clue2=[22 20 22 16];
clue3=[4 21 15 12];
clue4=[22 3 11 19];
knowncenter=[18 6 6 11];


options=double(ones(256,1024)); %location=(row-1)*16+col,  index valid table


tempv=1:256;

if clue_string(3)=='1'

    index=whichisClue(clue3,original);
    options(location(14,14),:)=0;
    options(location(14,14),index)=1;
    indexv=floor((index-1)/4)*4+1:floor((index-1)/4)*4+4;
    options(tempv(tempv~=location(14,14)),indexv)=0;
end
if clue_string(1)=='1'
    index=whichisClue(clue1,original);
    options(location(14,3),:)=0;
    options(location(14,3),index)=1;
    indexv=floor((index-1)/4)*4+1:floor((index-1)/4)*4+4;
    options(tempv(tempv~=location(14,3)),indexv)=0;
   
end

if clue_string(2)=='1'
    index=whichisClue(clue2,original);
    options(location(3,14),:)=0;
    options(location(3,14),index)=1;
    indexv=floor((index-1)/4)*4+1:floor((index-1)/4)*4+4;
    options(tempv(tempv~=location(3,14)),indexv)=0;
end

if clue_string(4)=='1'
    index=whichisClue(clue4,original);
    options(location(3,3),:)=0;
    options(location(3,3),index)=1;
    indexv=floor((index-1)/4)*4+1:floor((index-1)/4)*4+4;
    options(tempv(tempv~=location(3,3)),indexv)=0;
end

if planv~=3
    index=whichisClue(knowncenter,original);
    options(location(9,8),:)=0;
    options(location(9,8),index)=1;
    indexv=floor((index-1)/4)*4+1:floor((index-1)/4)*4+4;
    options(tempv(tempv~=location(9,8)),indexv)=0;
end

for row=2:15
    for col=2:15
        options(location(row,col),1:16*4)=0;
        options(location(row,col),(256-16)*4+1:256*4)=0;
        for row2=2:15  %corrected 13sep
            options(location(row,col),(location(row2,1)-1)*4+1:(location(row2,1)-1)*4+4)=0;
            options(location(row,col),(location(row2,16)-1)*4+1:(location(row2,16)-1)*4+4)=0;
        end
    end
    
end
%define options to margin locations
for r=1:16
    options(location(r,1),:)=0;
    options(location(r,16),:)=0;
end
for c=2:15
    options(location(1,c),:)=0;
    options(location(16,c),:)=0;
end
options(:,1:4:15*4+1)=0;
options(:,4*16*15+3: 4 :4*16*15+15*4+3)=0;
options(:,2: 4*16 :4*16*15+2)=0;
options(:,4*15+4: 4*16 :4*16*15+4*15+4)=0;

%corners may be exchanged after rotation
options(location(1,1),1)=1; %use left up
options(location(1,1),4*15+4)=1; %use right up
options(location(1,1),4*16*15+2)=1; %use left down
options(location(1,1),4*16*15+4*15+3)=1; %use right down
options(location(1,16),2)=1; %use left up
options(location(1,16),4*15+1)=1; %use right up
options(location(1,16),4*16*15+3)=1; %use left down
options(location(1,16),4*16*15+4*15+4)=1; %use right down
options(location(16,1),4)=1; %use left up
options(location(16,1),4*15+3)=1; %use right up
options(location(16,1),4*16*15+1)=1; %use left down
options(location(16,1),4*16*15+4*15+2)=1; %use right down
options(location(16,16),3)=1; %use left up
options(location(16,16),4*15+2)=1; %use right up
options(location(16,16),4*16*15+4)=1; %use left down
options(location(16,16),4*16*15+4*15+1)=1; %use right down
%margin-non-corner pieces may be excahnge after rotation
for c=2:15
    options(location(1,c),4+1:4:14*4+1)=1;
    options(location(1,c),4*16*15+4+3: 4 :4*16*15+14*4+3)=1;
    options(location(1,c),4*16+2: 4*16 :4*16*14+2)=1;
    options(location(1,c),4*16+4*15+4: 4*16 :4*16*14+4*15+4)=1;
    
    options(location(16,c),4*16*15+4+1:4:4*16*15+14*4+1)=1;
    options(location(16,c),4+3:4:14*4+3)=1;
    options(location(16,c),4*16+4: 4*16 :4*16*14+4)=1;
    options(location(16,c),4*16+4*15+2: 4*16 :4*16*14+4*15+2)=1;
end
for r=2:15
    options(location(r,1),4+4:4:14*4+4)=1;
    options(location(r,1),4*16*15+4+2: 4 :4*16*15+14*4+2)=1;
    options(location(r,1),4*16+1: 4*16 :4*16*14+1)=1;
    options(location(r,1),4*16+4*15+3: 4*16 :4*16*14+4*15+3)=1;
    
    options(location(r,16),4*16*15+4+4:4:4*16*15+14*4+4)=1;
    options(location(r,16),4+2:4:14*4+2)=1;
    options(location(r,16),4*16+3: 4*16 :4*16*14+3)=1;
    options(location(r,16),4*16+4*15+1: 4*16 :4*16*14+4*15+1)=1;
end
%for ind_location=1:256
%    if sum(options(ind_location,:)>0)==1
%        indf=1:1024;
%        indu=1:256;
%        indchoice=indf(options(ind_location,:)>0);
%        indstart=floor((indchoice-1)/4)*4+1;
%        options(indu~=ind_location,indstart:indstart+3)=0;
%    end
%end

baseline_options=options;
baseline_tree=tree;

if ~isempty(keep_options)
    options=keep_options;
else
    [~,options]=analyze(0,tree,treesize,options); %analyze up to level 0 (logics)
    baseline_options=options;
end

if planv==0
    required_options=zeros(256,1024);
    for t=1:256
        required_options(t,(t-1)*4+1)=1;
    end
    RRR=lsqminnorm(double(options),ones(256,1));
    chosenRRR=RRR;
    vec4=1:4;
    for ind=1:4:1024
        ttt=RRR(ind:ind+3);
        t=max(vec4(ttt==max(ttt)));
        yyy=[0 0 0 0];
        yyy(t)=1;
        chosenRRR(ind:ind+3)=yyy;
    end
    firstcheck_sol=sum(sum((options>0.5),2))
    firstcheck_actual_sol=sum(sum((options>0.5).*(required_options>0.5),2))
    firstcheck_algebraic_solved=sum(sum(((options>0.5).*(chosenRRR>0)'),2))
    firstcheck_actual_algebraic_solved=sum(sum(((options>0.5).*(chosenRRR>0)').*(required_options>0.5)))
    
end


%options=options.*(chosenRRR>0)';



for ind_location=1:256
    if sum(options(ind_location,:)>0)==1
        indf=1:1024;
        indu=1:256;
        indchoice=indf(options(ind_location,:)>0);
        indstart=floor((indchoice-1)/4)*4+1;
        options(indu~=ind_location,indstart:indstart+3)=0;
    end
end



if planv==0
    actual_sol=sum(sum((options>0.5).*(required_options>0.5),2))

end
sol=sum(sum(options>0.5,2))

indf=1:1024;
imglin=sum(options>0,2);
show_state=1*(imglin==1);
show_index=zeros(size(show_state));
for ind_location=1:256
    if show_state(ind_location)==1
       show_index(ind_location)=min(indf(options(ind_location,:)>0));
    end
end


keep_options=options;
keep_show_index=show_index;
keep_show_state=show_state;
ghost_options=[];

parallelstate=false;
delete(gcp('nocreate'));

symbollabel=['@' 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];

%RRR=lsqminnorm(double(options),ones(256,1));
%agebraic_solved=sum(RRR>0.5)

%h.HTMLSource=htmltxt;

flagloop=false;

for testno=1:10000


    %SHOWING SESSION
   % neworiginal=fill_neworiginal(show_state,show_index,options,piece);
   % show_original(h,fig,neworiginal,testno);

    %CH = getkey(1); 
    %screenCoord(testno,:) = get(0,'PointerLocation');
    if flagloop==false && machine_state==0
       % while (sum(screenCoord)==0)
       %     pause(0.5);
       % end
        %disp(screenCoord);
    
        CH=input('input: ','s');
        if CH=='!'
            row=input('row= ');
            col=input('col=');
            CH=input('applied input: ','s');
        end
    else
        col=1;
        row=1;
        CH='';
    end

    

    if ~flagloop
        if machine_state==5 
            machine_state=1;
        elseif machine_state==4
            machine_state=machine_state+1;
            h_number=3;
        elseif machine_state==3
            machine_state=machine_state+1;
            h_number=2;
        elseif machine_state==2
            machine_state=machine_state+1;
            %do W
            %h_number=1;
        elseif machine_state==1
            %for indlocation=1:256
            %    c=mod((indlocation-1),16)+1;
            %    r=floor((indlocation-1)/16)+1;
            %    if   true% ~((r==1 && c==1) || (r==1 && c==16) || (r==16 && c==1) || (r==16 && c==16))
            %        show_state(indlocation)=0;
            %        show_index(indlocation)=0;
            %        options(indlocation,:)=baseline_options(indlocation,:);
            %    end
            %end
            options=keep_options;
            show_index=keep_show_index;
            show_state=keep_show_state;
            for t=1:256
                if show_state(t)==0
                    options(t,:)=baseline_options(t,:);
                end
            end
            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            machine_state=machine_state+1;
        end
    end

    if strcmp(CH,'~')
        close(fig);
        return;
    end


    if strcmp(CH,'period')
        disp('Filtering options by selecting random orientations and using AI to predict actual connections. ')
        noiterations=input('Enter number of iterations: ');
        for iter=1:noiterations
            tree=update_tree(tree,show_state,show_index);
            [tree,treesize]=consult_AI(piece,tree,treesize,e2folder,planv); 
            [options]=correct_options(options,tree);
            imglin=sum(options>0,2);
            for indlocation=1:256
                if imglin(indlocation)==0
                    options(indlocation,:)=baseline_options(indlocation,:);
                end
            end
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                    show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
        end
    end



    if col>=1 && col<=16 && row>=1 && row<=16
        ind_location=(col-1)+(row-1)*16+1;
        if strcmp(CH,'t') 
            text_clue=input('Enter piece as 4- capital letter text: ','s');
            clue=[0 0 0 0];
            for t=1:4
                clue(t)=strfind(symbollabel,text_clue(t))-1;
            end
            indchosen=whichisClue(clue,original);
            if indchosen>0
                indstart=floor((indchosen-1)/4)*4+1;
                options(:,indstart:indstart+3)=0;
                options(ind_location,:)=0;
                options(ind_location,indchosen)=1;
                show_state(ind_location)=3;  %imposed state
                show_index(ind_location)=min(indf(options(ind_location,:)>0));

            else
                disp('Such piece does not exist')
            end
        end
        if strcmp(CH,'space') 
            show_state(ind_location)=2;  %open all natural options
            if sum(options(ind_location,:)>0)==0
                %options(ind_location,:)=baseline_options(ind_location,:);
            else
                show_index(ind_location)=min(indf(options(ind_location,:)>0));
            end
        end
        if strcmp(CH,'slash') 
            slash_mode=~slash_mode;
            disp(sprintf('Using d-wave=  %d',slash_mode))
            pause(1);
        end
        if strcmp(CH,'n')
            folderPath=input('Choose a folder with candidate mat files to score: ','s');
            filename = input('File name framework.mat in Ocean/QA folder: ','s');%[filename,path] = uigetfile('/gpfs01/home/seifer/eternity30/Ocean/QA/framework*.mat','Fetch frameowrk.mat file with node mapping');
            filen=[OceanFolder '/QA/' filename];
            filename30 = input('File name _server.py in Ocean folder: ','s');% uigetfile('/gpfs01/home/seifer/eternity30/Ocean/*.py','Fetch python server to run in background, based on modified ask_d_wave.py that you should modify according to comment in this script');
            pythonserverfile=[OceanFolder '/' filename30];      %[OceanFolder '/HamiltonianServerLinuxQA7269.py'];

            timestamp = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
            writematrix('stop',[OceanFolder '/command2server' timestamp '.txt']);
            writematrix('0',[OceanFolder '/replyfromserver' timestamp '.txt']);
            outnumbers_count=1;
            load(filen); %load (NodeSpin), mapnodes_ind, mapnodes_loc, options (not needed)
		    commandext = "bash -c ""module load mathlib mkl/2023 && module load python/2023 && source " + OceanFolder +"/venv/bin/activate && nohup " + OceanFolder + "/venv/bin/python "+ pythonserverfile + " "+ timestamp+ " > output.log 2>&1 & "" &";
		    disp(commandext);
		    [status, cmdout] = system(commandext);
		    disp(status);
		    disp('Output:');
		    disp(cmdout);

            fileList = dir(fullfile(folderPath, '*'));  % Use '*.*' to include all files
            Hamilton_values=zeros(1,length(fileList));
            for kk = 1:length(fileList)
                if fileList(kk).isdir
                   continue;
                else
                    fullFilePath = fullfile(folderPath, fileList(kk).name);
                    %fprintf('Processing file: %s\n', fullFilePath);
                    options=[];
                    best_options=[];
                    load(fullFilePath);
                    if ~isempty(best_options)
                        toptions=best_options;
                    end
                    if ~isempty(options)
                        toptions=options;
                    end
                    %leave only 2 lines of frame
                    for indlocation=1:256
                        c=mod((indlocation-1),16)+1;
                        r=floor((indlocation-1)/16)+1;
                        if ~(r==1 || r==16 || c==1 ||c==16 || r==2 || r==15 || c==2 || c==15 || (c==3 && r==3) || (c==3 && r==14) || (c==14 && r==3) || (c==14 && r==14))
                            toptions(indlocation,:)=baseline_options(indlocation,:);
                        end
                    end
                    timglin=sum(toptions>0,2);
                    tshow_state=1*(timglin==1);
                    tshow_index=zeros(size(tshow_state));
                    for indlocation=1:256
                        if tshow_state(indlocation)==1
                           tshow_index(indlocation)=min(indf(toptions(indlocation,:)>0));
                        end
                    end
                    ttree=update_tree(baseline_tree,tshow_state,tshow_index);
                    % W operation
                    [~,toptions]=analyze(0,ttree,treesize,toptions); 

                    %Here the python program should be already
                    %running, just communicate with files
                    %use: indchosen, indlocation, options,
                    bestspin=1*ones(size(mapnodes_loc'));
                    for nd=1:length(mapnodes_loc)
                        tloc=mapnodes_loc(nd);
                        tind=mapnodes_ind(nd);
                        bestspin(nd)=1-2*toptions(tloc,tind);
                    end
                    filenamebest=[OceanFolder '/QA/spin' timestamp '.csv'];
                    writematrix(bestspin,filenamebest);
                    writematrix('run',[OceanFolder '/command2server' timestamp '.txt']);
                    while true
                        outnumbers=readmatrix([OceanFolder '/replyfromserver' timestamp '.txt']);
                        if length(outnumbers)>outnumbers_count
                            outnumbers_count=length(outnumbers);
                            Hamiltonian=outnumbers(outnumbers_count);
                            break;
                        end
                        pause(3);
                    end
                    Hamilton_values(kk)=Hamiltonian;
                    disp(['Hamiltonian='  Hamiltonian '  File=' fileList(kk).name]);
                end

            end %for kk

            [~,ind_list]=sort(Hamilton_values);
            fidH = fopen([OceanFolder '/H_competition_' timestamp '.txt'], 'w');
            for k = 1:length(fileList)
                kk=ind_list(k);
                if fileList(kk).isdir
                   continue;
                else
                    fprintf(fidH, '%g, %s\n',Hamilton_values(kk),fileList(kk).name );
                end
            end
            fclose(fidH);

            pause(1);
        end
        if strcmp(CH,'return') && sum(options(ind_location,:)>0)>0
            show_state(ind_location)=1; %open recommeded option
            options(ind_location,:)=keep_options(ind_location,:);
            show_index(ind_location)=min(indf(options(ind_location,:)>0));
        end
        if strcmp(CH,'o') 
            if compromize
                compromize=false;
                disp('Compromize= NO');
            else
                compromize=true;
                disp('Compromize= YES');
            end
        end
        if strcmp(CH,'m') 
            if machine_state==0
                machine_state=1;
                disp('machine state= Yes');
                pause(0.5);
            else
                machine_state=0;
                disp('machine state= No');
                h_number=-1;
                flagloop=false;
                pause(0.5);
            end
        end
        if strcmp(CH,'y') || machine_state==3 
            for indlocation=1:256
                c=mod((indlocation-1),16)+1;
                r=floor((indlocation-1)/16)+1;
                if ~(r==1 || r==16 || c==1 ||c==16)
                    show_state(indlocation)=0;
                    show_index(indlocation)=0;
                    options(indlocation,:)=baseline_options(indlocation,:);
                end
            end
            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
        end
        if strcmp(CH,'u') 
            for indlocation=1:256
                c=mod((indlocation-1),16)+1;
                r=floor((indlocation-1)/16)+1;
                if ~(r==1 || r==16 || c==1 ||c==16 || r==2 || r==15 || c==2 || c==15 || (c==3 && r==3) || (c==3 && r==14) || (c==14 && r==3) || (c==14 && r==14))
                    show_state(indlocation)=0;
                    show_index(indlocation)=0;
                    options(indlocation,:)=baseline_options(indlocation,:);
                end
            end
            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
        end
        if strcmp(CH,'j') %toggle between probablistic and discrete
            if isprob==1
                options=prob_filter(options,baseline_options,treesize,tree);
                options=options>0;
                isprob=0;
                disp('discrete mode');
            else
                sbuff=sum(options,2);
                options=options./(sum(options,2)+(sbuff==0));
                isprob=1;
                disp('probablistic mode');
            end

        end
        if strcmp(CH,'e')   %export txt file to eternityII program
            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            indu=1:256;
            for indlocation=1:256
                if show_state(indlocation)==0
                    options(indlocation,:)=baseline_options(indlocation,:);
                end
            end
            for indlocation=1:256
                if show_state(indlocation)>0
                    options(indu~=indlocation,floor((show_index(indlocation)-1)/4)*4+1:floor((show_index(indlocation)-1)/4)*4+4)=0;
                end
            end
            for indlocation=1:256
                if show_state(indlocation)==0
                    indxcc=indf(options(indlocation,:)>0);
                    findrandom=floor(rand*length(indxcc))+1;
                    options(indlocation,:)=0;
                    options(:,floor((indxcc(findrandom)-1)/4)*4+1:floor((indxcc(findrandom)-1)/4)*4+4)=0;
                    options(indlocation,indxcc(findrandom))=1;
                    show_state(indlocation)=1;
                    show_index(indlocation)=indxcc(findrandom);
                end
            end
            filesave=sprintf('%s//savestate%d.txt',e2folder,testno);
            op=translate2external(piece,show_state,show_index);
            writematrix(op,filesave,'Delimiter','space');
            disp(sprintf('saving savestate%d.txt',testno));
            pause(1);
        end
        if strcmp(CH,'s') 
            %options=force_to_options(options,baseline_options);
            filesave=sprintf('%s//savestate%d_C%s.mat',e2folder,testno,clue_string);
            save(filesave,'options','-mat');
            disp(sprintf('saving savestate%d',testno));
            pause(1);
        end
        if strcmp(CH,'hyphen') || strcmp(CH,'1') || strcmp(CH,'2') || strcmp(CH,'3') || strcmp(CH,'4')  || strcmp(CH,'5')  || strcmp(CH,'6') || strcmp(CH,'7') || strcmp(CH,'8') || strcmp(CH,'9')
            %options=force_to_options(options,baseline_options);
            disp('Effecting options table of entire board..');
            deepfill=true;
            filenumber_ready=0;
            if strcmp(CH,'1') 
                problemsize=1024;
            elseif strcmp(CH,'2')
                problemsize=2048;
            elseif strcmp(CH,'3')
                problemsize=4096;
            elseif strcmp(CH,'4')
                problemsize=10000;
            elseif strcmp(CH,'5')
                problemsize=20000;
            elseif strcmp(CH,'6')
                problemsize=30000;
            elseif strcmp(CH,'7')
                problemsize=60000;
            elseif strcmp(CH,'8')
                problemsize=300000;
            elseif strcmp(CH,'9')
                problemsize=1000000;
            elseif strcmp(CH,'hyphen')
                filenumber_ready=input('Enter QA file number to retrieve:  ');
                problemsize=200000; %dummy
            end
            shuffled_choices_pick=indf(options(ind_location,:)>0);
            [~, ghost_options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,filenumber_ready);
            if ~slash_mode %if not during slash mode
                options=ghost_options;
                imglin=sum(options>0,2);
                show_state=1*(imglin==1);
                show_index=zeros(size(show_state));
                for indlocation=1:bS
                    if show_state(indlocation)==1
                       show_index(indlocation)=min(indf(options(indlocation,:)>0));
                    end
                end
                tree=update_tree(tree,show_state,show_index);
                disp('options table updated by dwave computation.');
            else %slash mode, so use dwave to prepare for h function
                disp('ghost_options updated by dwave calculation, which can be used in h function.')
                filesave=sprintf('%s//dwave_latest.mat',e2folder);
                save(filesave,'ghost_options','-mat');
                disp('ghost_options saved.');
                slash_mode=false;
                dwimglin=sum(ghost_options>0,2);
                dwshow_state=1*(dwimglin==1);
                dwshow_index=zeros(size(dwshow_state));
                for dwindlocation=1:256
                    if dwshow_state(dwindlocation)==1
                       dwshow_index(dwindlocation)=min(indf(ghost_options(dwindlocation,:)>0));
                    end
                end
                [ghost_options,dwshow_state,dwshow_index]=force_to_options(ghost_options,options,dwshow_state,dwshow_index);
                disp('ghost_options regularized');        
                disp('slash function deactivated:  slash_mode=false.');
            end
            pause(1);
        end
        if strcmp(CH,'equal')
            QAcounter=input('Enter QAcounter of answer from dwave, to project on options: ');
            options=reproject_dwave_answer(options,QAcounter);
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
            tree=update_tree(tree,show_state,show_index);
            disp('options table updated according to dwave computation, including mat table of the QA slot.');
            pause(1);
        end
        if strcmp(CH,'multiply')
            if PaidDwaveService==false
                PaidDwaveService=true;
                disp('Quantum annealing calculation will use the paid subscription');
            else
                PaidDwaveService=false;
                disp('Simulated quantum annealing will use the local computer for free');
            end
            pause(1);
        end
        if strcmp(CH,'0')
            disp('Effecting only the chosen piece, and check for consistency');
            deepfill=true;
            bakup_options=options;
            row_in_options0=options(ind_location,:);
            shuffled_choices_pick=indf(options(ind_location,:)>0);
            problemsize=input('How many qubits to invest in the start run? ');
            problemsize0=problemsize;
            [~, options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,0);
            row_in_options1=options(ind_location,:);
            options=bakup_options;
            if sum(row_in_options1==row_in_options0)<1024
                problemsize=problemsize0*2;
                [~, options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,0);
                row_in_options2=options(ind_location,:);
                options=bakup_options;
                if sum(row_in_options1==row_in_options2)==1024
                    options(ind_location,:)=row_in_options2; %accept verified update in selections of the pointed piece
                    disp('Options table has been updated in the location of the piece according to verified dwave computation.');
                else
                    disp('No update')
                    do_more=input('Check with more qutbits? (0-no, 1-yes): ');
                    if do_more==1
                        problemsize=problemsize0*4;
                        [~, options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,0);
                        row_in_options3=options(ind_location,:);
                        options=bakup_options;
                        if sum(row_in_options2==row_in_options3)==1024
                            options(ind_location,:)=row_in_options3; %accept verified update in selections of the pointed piece
                            disp('Options table has been updated in the location of the piece according to verified dwave computation.');
                        else
                            disp('No update also in further trial')
                        end
                    end
                end
            else %NOT if sum(row_in_options1==row_in_options0)<1024, meaning that QM did not eliminate any option
                disp('Trial shows no benefit from QM.')
            end

            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
            tree=update_tree(tree,show_state,show_index);
            pause(1);
        end
        if strcmp(CH,'p') 
            disp(sprintf('Start count= %d',sum(show_state>0)));
            NumLevels=input('Enter number of levels to analyze? ');
            if ~parallelstate
                startSafeParpool(parallelForce);
                parallelstate=true;
            end
            if isprob==0
                sbuff=sum(options,2);
                options=options./(sum(options,2)+(sbuff==0));
                tree=baseline_tree;
                tree=update_tree(tree,show_state,show_index);
            end
            [score,options]=analyze_prob(NumLevels,tree,treesize,options);
            %%options=prob_filter(options,baseline_options,treesize,tree);
            %%options=options>0;
            isprob=1;
            disp('probablistic mode');
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            disp(sprintf('End count= %d',sum(show_state>0)));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end

            RRR=lsqminnorm(double(options),ones(256,1));
            extra=max(sum(RRR),0);
            algb_score=abs(256-extra)+1;
            rank_score=rank(double(options));
            disp(sprintf('algebraic score=%d,  Rank=%d',algb_score,rank_score));

        end
        if all(strcmp(CH,'x')  || strcmp(CH,'v'))%search
            if all( show_state(ind_location)==0 || show_state(ind_location)==4)
                if sum(options(ind_location,:))==0
                    options(ind_location,:)=baseline_options(ind_location,:);
                end
            end
            show_state(ind_location)=4;
            tree=baseline_tree;
            tree=update_tree(tree,show_state,show_index);

            for ind=1:1024
                if options(ind_location,ind)>0
                    piece1=1+floor((ind-1)/4);
                    rot1=ind-4*(piece1-1);
                    flag=true;
                    for direction=1:4
                        sym=piece(piece1,mod(direction-1-(rot1-1),4)+1);
                        if direction==1
                            next_row=row-1;
                            next_col=col;
                        elseif direction==2
                            next_row=row;
                            next_col=col+1;
                        elseif direction==3
                            next_row=row+1;
                            next_col=col;
                        elseif direction==4
                            next_row=row;
                            next_col=col-1;
                        end
                        if all(next_row==0 || next_row==17 || next_col==0 || next_col==17)
                            next_sym=0;
                        else
                            next_location=(next_col-1)+(next_row-1)*16+1;
                            indf=1:1024;
                            next_index_solved_vector=indf(options(next_location,:)>0);
                            if show_state(next_location)>0 && sum(next_index_solved_vector==show_index(next_location))>0 %user selected one choice
                                index_piece_v=show_index(next_location);
                            else
                                index_piece_v=indf(options(next_location,:)>0); %take all the options of neighbor pieces
                            end
                            piece1_v=1+floor((index_piece_v-1)/4);
                            rot1_v=index_piece_v-4*(piece1_v-1);
                            next_sym=piece(sub2ind(size(piece),piece1_v,mod(direction-1+2-(rot1_v-1),4)+1));
                        end
                        flag=flag && sum(next_sym==sym)>0;
                    end
                    if flag==false
                        options(ind_location,ind)=0;
                    end
                end
            end
            index_solved_vector=indf(options(ind_location,:)>0);
            indu=1:256;
            %remove choices of pieces already shown in other positions
            for tspos=1:length(index_solved_vector)
                indt=index_solved_vector(tspos);
                indt_rot1=floor((indt-1)/4)*4+1;
                for indt_rot=indt_rot1:indt_rot1+3
                    if sum(show_index(indu~=ind_location)==indt_rot & show_state(indu~=ind_location)>0)>0
                        options(ind_location,indt)=0;
                        index_solved_vector(tspos)=-1;
                    end
                end
            end
            index_solved_vector=index_solved_vector(index_solved_vector>0);
            if strcmp(CH,'v')
                if length(index_solved_vector)>0
                    show_index(ind_location)=min(index_solved_vector);
                else
                    show_state(ind_location)=0;
                    show_index(ind_location)=0;
                end
            else
                maxchoices= length(index_solved_vector);
                after_options=options;
                score=zeros(1,maxchoices);
                if ~parallelstate
                    startSafeParpool(parallelForce);
                    parallelstate=true;
                end
                parfor testchoice=1:maxchoices
                    indchosen=index_solved_vector(testchoice);
                    indstart=floor((indchosen-1)/4)*4+1;
                    toptions=options;
                    toptions(:,indstart:indstart+3)=0;
                    toptions(ind_location,:)=0;
                    toptions(ind_location,indchosen)=1;
                    [toptions,score_nostop]=new_options(compromize,toptions,tree);
                    score(testchoice)=sum(sum(toptions,2))*(score_nostop>0); 
    
                end
                best_score=inf;
                best_indx_x=1;
                for testchoice=1:maxchoices
                    indchosen=index_solved_vector(testchoice);
                    if score(testchoice)==0
                        after_options(ind_location,indchosen)=0;
                    else
                        if abs(score-256)<abs(best_score-256)
                            best_score=score;
                            best_indx_x=indchosen;
                        end
                    end
    
                end
                disp(sprintf('best score=%d ',best_score));
                options=after_options; %updates list of valid choices checked by the perturbation analayizis
                show_index(ind_location)=best_indx_x;
            end
            if sum(options(ind_location,:)>0)==0
                show_index(ind_location)=0;
                show_state(ind_location)=0;
                disp('No matching piece found');
            end
        end
        if strcmp(CH,'l')
            filename=input('File name in folder: ','s');
            path=[e2folder '/'];
            loadsave=[path filename];
            disp(['loading ' loadsave]);
            if ~slash_mode
                best_options=[];
                ghost_options=[];
                load(loadsave,'-mat');
                if ~isempty(best_options)
                    options=best_options;
                end
                if ~isempty(ghost_options)
                    options=ghost_options;
                end
                keep_options=options;
                imglin=sum(options>0,2);
                show_state=1*(imglin==1);
                show_index=zeros(size(show_state));
                for indlocation=1:256
                    if show_state(indlocation)==1
                       show_index(indlocation)=min(indf(options(indlocation,:)>0));
                    end
                end
            elseif slash_mode
                remoptions=options;
                remshow_state=show_state;
                remshow_index=show_index;
                best_options=[];
                ghost_options=[];
                load(loadsave,'-mat');
                if ~isempty(best_options)
                    options=best_options;
                end
                if ~isempty(ghost_options)
                    options=ghost_options;
                end
                %imglin=sum(options>0,2);
                %show_state=1*(imglin==1);
                %show_index=zeros(size(show_state));
                %for indlocation=1:256
                %    if show_state(indlocation)==1
                %       show_index(indlocation)=min(indf(options(indlocation,:)>0));
                %    end
                %end
                %[options,~,~]=force_to_options(options,remoptions,show_state,show_index);
                ghost_options=options;
                options=remoptions;
                show_state=remshow_state;
                show_index=remshow_index;
                disp('File loaded to ghost_options because this is a slash mode. Use h function next.')
                slash_mode=false;
                disp('Slash mode deactivated.');
            end

        end
        if strcmp(CH,'d')
            filename=input('File name in folder: ','s');%[filename,path] = uigetfile(sprintf('%s//*.mat',e2folder),'Fetch mat file after analysis');
            path=[e2folder '/'];
            loadsave=[path filename];
            disp(['adding ' loadsave]);
            adding_options=options;
            load(loadsave,'-mat');
            if isprob==0
                options=1*(adding_options & options);
            else
                sbuff=sum(options,2);
                options=options./(sum(options,2)+(sbuff==0));
                adding_options=adding_options./(sum(options,2)+(sbuff==0));
                options=adding_options + options;
            end
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end

        end
        if strcmp(CH,'c') 
            show_state(ind_location)=0; %hide, cancel
            show_index(ind_location)=0;
            options(ind_location,:)=baseline_options(ind_location,:);
        end
        if strcmp(CH,'z') && show_state(ind_location)>0
            indf=1:1024;
            index_solved_vector=indf(options(ind_location,:)>0);
            fn=[];
            for t=1:length(index_solved_vector)
                index_piece=index_solved_vector(t);
                piece1=1+floor((index_piece-1)/4);
                rot1=index_piece-4*(piece1-1);
                temp=[0 0 0 0];
                for direct1=1:4
                    temp(1,direct1)=piece(piece1,mod(direct1-1-(rot1-1),4)+1);
                end
                textt=[(symbollabel(temp(1,1)+1)) (symbollabel(temp(1,2)+1)) (symbollabel(temp(1,3)+1)) (symbollabel(temp(1,4)+1))];
                fn{t}=textt;
            end
          
            [indx,tf] = listdlg('ListString',fn,'SelectionMode','single'); 
            if ~isempty(indx)
                show_index(ind_location)=index_solved_vector(min(indx));
            end
        end
        if strcmp(CH,'r')
            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            indu=1:256;
            for indlocation=1:256
                if show_state(indlocation)==0
                    options(indlocation,:)=baseline_options(indlocation,:);
                end
            end
            for indlocation=1:256
                if show_state(indlocation)>0
                    options(indu~=indlocation,floor((show_index(indlocation)-1)/4)*4+1:floor((show_index(indlocation)-1)/4)*4+4)=0;
                end
            end
            for indlocation=1:256
                if show_state(indlocation)==0
                    indxcc=indf(options(indlocation,:)>0);
                    findrandom=floor(rand*length(indxcc))+1;
                    options(indlocation,:)=0;
                    options(:,floor((indxcc(findrandom)-1)/4)*4+1:floor((indxcc(findrandom)-1)/4)*4+4)=0;
                    options(indlocation,indxcc(findrandom))=1;
                    show_state(indlocation)=1;
                    show_index(indlocation)=indxcc(findrandom);
                end
            end
        end
        if strcmp(CH,'f')
            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            tree=update_tree(baseline_tree,show_state,show_index);
            [options]=correct_options(options,tree);

        end
        if strcmp(CH,'k')
            timlimit=input('Time limit in minutes? ')*60;
            options=byscore_reorder(original,options,show_state,show_index,timlimit);
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            disp(sprintf('End count= %d',sum(show_state>0)));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
        end


        if strcmp(CH,'backquote')

            disp('Starting relaxed full solution ');
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
            tree=update_tree(baseline_tree,show_state,show_index);


            if ~parallelstate
                startSafeParpool(parallelForce);
                parallelstate=true;
            end
            score_array=zeros(1,parallelForce);
            parfor threadno=1:parallelForce
                kseed = RandStream.create('mrg32k3a','Seed',10*threadno+floor(10*rand(1)));
                RandStream.setGlobalStream(kseed);
                [temp_options]=relaxed_solver(threadno,options,tree,treesize);
                temp_imglin=sum(temp_options>0,2);
                temp_show_state=1*(temp_imglin==1);
                temp_show_index=zeros(size(temp_show_state));
                for temp_indlocation=1:256
                    if temp_show_state(temp_indlocation)==1
                       temp_show_index(temp_indlocation)=min(indf(temp_options(temp_indlocation,:)>0));
                    end
                end
                [temp_options,temp_show_state,~]=force_to_options(temp_options,baseline_options,temp_show_state,temp_show_index);                
                grand_temp_options(:,:,threadno)=temp_options;
                score_array(threadno)=sum(temp_show_state==1);
            end
            filesave=sprintf('%s\\grand_options.mat',e2folder);
            save(filesave,'grand_temp_options','-mat');
            disp('Grand_options file saved');

            vect_th=1:parallelForce;
            disp(score_array);
            thbest=min(vect_th(score_array==max(score_array)));
            options=grand_temp_options(:,:,thbest);

            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end

            pause(1);
        end

        if strcmp(CH,'a')
            for indlocation=1:256
                if sum(options(indlocation,:)>0)==0
                    options(indlocation,:)=baseline_options(indlocation,:);
                end
            end
            counteliminate=input('How many pieces to eliminate at once? ');
            if counteliminate>1
                timlimit=input('Time limit in minutes? ')*60;
            else
                timlimit=30*60;
            end
            if ~parallelstate
                startSafeParpool(parallelForce);
                parallelstate=true;
            end

            indu=1:256;
            indlocation_group=indu(show_state(:)>0 & indu'~=35 & indu'~=46 & indu'~=136 & indu'~=211 & indu'~=222);
            indlocation_group_count=length(indlocation_group);
            target_count=factorial(indlocation_group_count)/(factorial(indlocation_group_count-counteliminate)*factorial(counteliminate));
            best_options=options;

            topar_elimind=zeros(parallelForce,counteliminate);
            topar_bestcount=zeros(parallelForce,1);
            parind_v=1:parallelForce;
            parfor testind=parind_v 
                showcount=0;
                temp_show_state=zeros(256,1);
                temp_show_index=zeros(256,1);
                temp_baseline_options=baseline_options;
                indff=1:1024;
                tic;
                for trynumber=1:min(target_count,1e9)
                    if toc>timlimit
                        break;
                    end
                    temp_options=options;
                    temp_indlocation_group=indlocation_group;
                    temp_elim_ind=[];
                    %try to remove pieces, hoping to get only true pieces
                    for elim_number=1:counteliminate
                        ind_pos=floor(rand*length(temp_indlocation_group))+1;
                        ind=temp_indlocation_group(ind_pos);
                        temp_options(ind,:)=temp_baseline_options(ind,:);
                        temp_indlocation_group=[temp_indlocation_group(1:ind_pos-1) temp_indlocation_group(ind_pos+1:end)];
                        temp_elim_ind=[temp_elim_ind ind];
                    end
    
                    imglin=sum(temp_options>0,2);
                    temp_show_state=1*(imglin==1);
                    for indlocation=1:256
                        if temp_show_state(indlocation)==1
                            temp_show_index(indlocation)=min(indff(temp_options(indlocation,:)>0));
                        end
                    end
                    % Q0 once
                    temp_tree=baseline_tree;
                    temp_tree=update_tree(temp_tree,temp_show_state,temp_show_index);
                    [temp_score,temp_options]=analyze(0,temp_tree,treesize,temp_options); %analyze up to level 0 (logics)
                    %temp_tree=update_tree(temp_tree,temp_show_state,temp_show_index);
                    %[temp_score,temp_options]=analyze(0,temp_tree,treesize,temp_options); %analyze up to level 0 (logics)

                    imglin=sum(temp_options>0,2);
                    temp_show_state=1*(imglin==1);
                    showcount=sum(temp_show_state>0);
                    if temp_score>0 && showcount>topar_bestcount(testind)
                        topar_bestcount(testind)=showcount;
                        topar_elimind(testind,:)=temp_elim_ind;
                    end
                end
            end %parfor
            parindbest=min(parind_v(topar_bestcount==max(topar_bestcount)));
            bestcount=topar_bestcount(parindbest);
            score=0;
            if bestcount>0
                best_elimind=topar_elimind(parindbest,:);
                for t=1:counteliminate
                    ind=best_elimind(1,t);
                    options(ind,:)=baseline_options(ind,:);
                end
                imglin=sum(options>0,2);
                show_state=1*(imglin==1);
                show_index=zeros(size(show_state));
                for indlocation=1:256
                    if show_state(indlocation)==1
                        show_index(indlocation)=min(indf(options(indlocation,:)>0));
                    end
                end
                tree=baseline_tree;
                tree=update_tree(tree,show_state,show_index);
                [score,options]=analyze(0,tree,treesize,options); %analyze up to level 0 (logics)
            end
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            showcount=sum(show_state>0);
            for indlocation=1:256
                if show_state(indlocation)==1
                    show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
            disp(sprintf('count= %d,  score=%d',showcount,score));
        end
        if strcmp(CH,'h') 
            h_number=input('Enter number of frame lines to fill = ');
            flag_snail=input('Which search pattern? (0-brick, 1-snail): ');

            writelines('start',sprintf('%s//break.txt',e2folder));
			flag_useHamiltonian=input('Use Ising Hamiltonian? (0-no, 1-yes): ');
			if flag_useHamiltonian
				filename = input('File name framework.mat in Ocean/QA folder: ','s');%[filename,path] = uigetfile('/gpfs01/home/seifer/eternity30/Ocean/QA/framework*.mat','Fetch frameowrk.mat file with node mapping');
				filen=[OceanFolder '/QA/' filename];
				filename30 = input('File name _server.py in Ocean folder: ','s');% uigetfile('/gpfs01/home/seifer/eternity30/Ocean/*.py','Fetch python server to run in background, based on modified ask_d_wave.py that you should modify according to comment in this script');
				pythonserverfile=[OceanFolder '/' filename30];      %[OceanFolder '/HamiltonianServerLinuxQA7269.py'];
			else
				filen='';
				pythonserverfile='';
			end

            disp('Starting to solve for frame block ...');

            loop_options=options;
            loop_show_state=show_state;
            loop_show_index=show_index;

            streams = RandStream.create('mrg32k3a','NumStreams',parallelForce,'Seed',12345);

            if flag_useHamiltonian
					threadno=1;
                    kseed = RandStream.create('mrg32k3a','Seed',threadno);
                    RandStream.setGlobalStream(kseed);
                    %[options,show_state,show_index,parallelstate,parallelForce,failblock]=form_block(threadno,options,baseline_options,show_state,show_index,piece,baseline_tree,parallelstate,parallelForce,e2folder,compromize,isprob,h_number,slash_mode,ghost_options);
                    [paroptions,parshow_state,parshow_index,failblock]=form_block(threadno,kseed,options,baseline_options,show_state,show_index,piece,baseline_tree,parallelstate,parallelForce,e2folder,compromize,isprob,h_number,slash_mode,ghost_options,flag_useHamiltonian,filen,pythonserverfile,OceanFolder,flag_snail);
            else
				if ~parallelstate
					startSafeParpool(parallelForce);
					parallelstate=true;
                end

                parfor threadno=1:parallelForce  
                    kseed = RandStream.create('mrg32k3a','Seed',threadno);
                    RandStream.setGlobalStream(kseed);
                    %[options,show_state,show_index,parallelstate,parallelForce,failblock]=form_block(threadno,options,baseline_options,show_state,show_index,piece,baseline_tree,parallelstate,parallelForce,e2folder,compromize,isprob,h_number,slash_mode,ghost_options);
                    [paroptions,parshow_state,parshow_index,failblock]=form_block(threadno,kseed,options,baseline_options,show_state,show_index,piece,baseline_tree,parallelstate,parallelForce,e2folder,compromize,isprob,h_number,slash_mode,ghost_options,flag_useHamiltonian,filen,pythonserverfile,OceanFolder,flag_snail);
                    %if ~failblock
                    %    writelines('end',sprintf('%s//break.txt',e2folder));
                    %    s=struct("option",paroptions,"show_state",parshow_state,"show_index",parshow_index);
                    %    save(sprintf('%s//passed.mat',e2folder),'-fromstruct',s);
                    %end
                end
            end


            
            %if sum(show_state>0)>=230 || failblock
            %    flagloop=false;
            %end
            %disp(sprintf('Solved= %d',sum(show_state>0)));
            %if sum(show_state>0)>150
            %    filesave=sprintf('%s\\savesolve%d_C%s.mat',e2folder,testno,clue_string);
            %    save(filesave,'options','-mat');
            %    disp(sprintf('saving savestate%d',testno));
            %    pause(1);
            %end

        end
        if strcmp(CH,'w') || machine_state==3
            disp('W');
            if isprob==0
                %temp_tree=baseline_tree;
                %tree=update_tree(temp_tree,show_state,show_index);
                tree=update_tree(tree,show_state,show_index);
                [score,options]=analyze(0,tree,treesize,options); %analyze up to level 0 (logics)
            else
                [score,options]=analyze_prob(0,tree,treesize,options); %analyze up to level 0 (logics)
            end
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
            if isprob==0
                [options,show_state,show_index]=force_to_options(options,options,show_state,show_index);
                tree=update_tree(tree,show_state,show_index);
            end
            %show_state=1*(imglin==1);
            disp(sprintf('End count= %d',sum(show_state>0)));
            RRR=lsqminnorm(double(options),ones(256,1));
            extra=max(sum(RRR),0);
            algb_score=abs(256-extra)+1;
            rank_score=rank(double(options));
            disp(sprintf('algebraic score=%d,  Rank=%d',algb_score,rank_score));
            if planv==0
                required_options=zeros(256,1024);
                for t=1:256
                    required_options(t,(t-1)*4+1)=1;
                end
                firstcheck_actual_sol=sum(sum((options>0.5).*(required_options>0.5),2));
                disp(sprintf('Still could be solved= %d  ',firstcheck_actual_sol));
            end

        end
        if strcmp(CH,'q')

            %%[options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            disp(sprintf('Start count= %d',sum(show_state>0)));
            NumLevels=input('Enter number of levels to analyze? ');
            if ~parallelstate
                startSafeParpool(parallelForce);
                parallelstate=true;
            end
            tree=baseline_tree;
            tree=update_tree(tree,show_state,show_index);
            [score,options]=analyze(NumLevels,tree,treesize,options); %analyze up to level=1 (chose two elminiations as well)
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
            %%[options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            show_state=1*(imglin==1);
            disp(sprintf('End count= %d',sum(show_state>0)));
            RRR=lsqminnorm(double(options),ones(256,1));
            extra=max(sum(RRR),0);
            algb_score=abs(256-extra)+1;
            rank_score=rank(double(options));
            disp(sprintf('algebraic score=%d,  Rank=%d',algb_score,rank_score));

        end
    end
    screenCoord=[0 0];
    pause(0.5);
end

end %main function


function [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index)
            %force the displayed pieces chosen by user to be one option
            %only. Meaning they are now clues.
            indf=1:1024;
            for indlocation=1:256
                if show_state(indlocation)>0 
                    index_solved_vector=indf(options(indlocation,:)>0);
                    if sum(index_solved_vector==show_index(indlocation))>0 
                        indchosen=show_index(indlocation);
                        indstart=floor((indchosen-1)/4)*4+1;
                        options(:,indstart:indstart+3)=0;
                        options(indlocation,:)=0;
                        options(indlocation,indchosen)=1;
                    else
                        show_state(indlocation)=0;
                        show_index(indlocation)=0;
                    end
                end
            end
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:256
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                else
                    show_index(indlocation)=0;
                end
            end

            %for indlocation=1:256
            %    if show_state(indlocation)==0 
            %        options(indlocation,:)=baseline_options(indlocation,:);
            %    end
            %end
            for indlocation=1:256
                if show_state(indlocation)>0 
                    indf=1:1024;
                    index_solved_vector=indf(options(indlocation,:)>0);
                    if sum(index_solved_vector==show_index(indlocation))>0 
                        indchosen=show_index(indlocation);
                        indstart=floor((indchosen-1)/4)*4+1;
                        options(:,indstart:indstart+3)=0;
                        options(indlocation,:)=0;
                        options(indlocation,indchosen)=1;
                    else
                        show_state(indlocation)=0;
                        show_index(indlocation)=0;
                    end
                end
            end

            for rpt=1:5
                for indlocation=1:256
                    if sum(options(indlocation,:)>0,2)==0 %useful after dwave to fill empty lines
                        options(indlocation,:)=baseline_options(indlocation,:);
                        show_state(indlocation)=0;
                    end
                end
                for indlocation=1:256
                    if sum(options(indlocation,:)>0)==1 && show_state(indlocation)>0
                        indf=1:1024;
                        indu=1:256;
                        indchoice=indf(options(indlocation,:)>0);
                        show_index(indlocation)=indchoice;
                        indstart=floor((indchoice-1)/4)*4+1;
                        options(indu~=indlocation,indstart:indstart+3)=0;
                    end
                end
            end

end

function [options,show_state,show_index]=force_to_home(options,baseline_options,show_state,show_index)
    %useful for h function
    for indlocation=1:256
        if show_state(indlocation)==0
            options(indlocation,:)=baseline_options(indlocation,:);
        end
    end
    for rpt=1:5
        for indlocation=1:256
            if sum(options(indlocation,:)>0,2)~=1
                options(indlocation,:)=baseline_options(indlocation,:);
                show_state(indlocation)=0;
            end
        end
        for indlocation=1:256
            if sum(options(indlocation,:)>0)==1 && show_state(indlocation)>0
                indf=1:1024;
                indu=1:256;
                indchoice=indf(options(indlocation,:)>0);
                show_index(indlocation)=indchoice;
                indstart=floor((indchoice-1)/4)*4+1;
                options(indu~=indlocation,indstart:indstart+3)=0;
            end
        end
    end

end


%now show in solution only correct choices
%options=options.*(actualRRR>0.5)';
%neworiginal=fill_neworiginal(options,piece);
%show_original(neworiginal);

    
    %save('C:\Users\seifer\Documents\Matlab\eternity\perturb3_raw',"perturb_vectors","perturb_param")
    
    %cov=zeros(index_vectors,index_vectors);
    %for row=1:index_vectors
    %    for col=1:index_vectors
    %        cov(row,col)=mean(mean(perturb_vectors(row,:).*perturb_vectors(col,:)));
    %    end
    %end
    %save('C:\Users\seifer\Documents\Matlab\eternity\perturb3_cov',"cov");
    
    %SVD -> U(:,1) and U(:,2) are the most important modes (eigen vectors)
    %[U,S,V] = svd(cov);
    %f=diag(S);
    %save('C:\Users\seifer\Documents\Matlab\eternity\perturb3_SVD',"f","U","V");
    
    %Invsigma=zeros(size(S));
    %maxsigma=max(max(S));
    %for ind=1:min(length(sigma(1,:)),length(sigma(:,1)))
    %    if sigma(ind,ind)>maxsigma*1e-10
    %        Invsigma(ind,ind)=(1/S(ind,ind));
    %        lastind=ind;
    %    else
    %        Invsigma(ind,ind)=0;
    %    end
    %end
    %PseudoInvS=V*Invsigma'*U';
    %result=PseudoInvS*Csample;
    

%end %end of main function (functions cannot be nested)




function ind=whichisClue(clue,original)
    for r=1:16
        for c=1:16
            for rot=1:4
                spiece=original{r,c};
                if rot==2
                   spiece=[spiece(4) spiece(1:3)]; %corrected 1 sep 24
                elseif rot==3
                    spiece=[spiece(3:4) spiece(1:2)];
                elseif rot==4
                    spiece=[spiece(2:4) spiece(1)];
                end
                if sum(spiece==clue)==4
                    ind=((r-1)*16+c-1)*4+rot;
                    return;
                end
            end
        end
    end
    ind=-1;
end

function locationr=location(r,c)
    locationr=16*(r-1)+c;
end



function neworiginal=fill_neworiginal(show_state,show_index,options,piece)
    indf=1:1024;
    for row=1:16
        for col=1:16
            ind_location=(col-1)+(row-1)*16+1;
            if show_state(ind_location)>0
                index_solved_vector=indf(options(ind_location,:)>0);
                if sum(index_solved_vector==show_index(ind_location))>0
                    index_piece=show_index(ind_location);
                    piece1=1+floor((index_piece-1)/4);
                    rot1=index_piece-4*(piece1-1);
                    temp=[0 0 0 0];
                    for direct1=1:4
                        temp(1,direct1)=piece(piece1,mod(direct1-1-(rot1-1),4)+1);
                    end
                    neworiginal{row,col}=temp;
                else
                    show_index(ind_location)=0;
                    show_state(ind_location)=0;
                    neworiginal{row,col}=[0 0 0 0];
                end
            else
                neworiginal{row,col}=[0 0 0 0];
            end
        end
    end

end




function [new_optionsr,score_nostop]=new_options(compromize,options,tree)

    new_optionsr=options>0;
    score_nostop=1;
    imglin=sum(new_optionsr,2);
    score_solved_prev=sum(imglin==1);

    indf=1:1024;
    imglin=sum(new_optionsr>0,2);
    show_state=1*(imglin==1);
    show_index=zeros(size(show_state));
    for indlocation=1:256
        if show_state(indlocation)==1
            show_index(indlocation)=min(indf(new_optionsr(indlocation,:)>0));
        end
    end
    tree=update_tree(tree,show_state,show_index); %I think already done, consider remove

    for clk=2:255   
       %%prev_sumoptions=sum(new_optionsr(:));
       prev_options=new_optionsr;
       vect_ind=1:256*4;
       %prev_options=options;
       for location_ind=1:256
            %if sum(new_optionsr(location_ind,1:4:end) | new_optionsr(location_ind,2:4:end) | new_optionsr(location_ind,3:4:end) | new_optionsr(location_ind,4:4:end) )==1 %certinty in one location
            %    ind=max(vect_ind(new_optionsr(location_ind,:)));
            %    remain_loc=[1:location_ind-1 location_ind+1:256];
            %    ind_p=floor((ind-1)/4)*4+1:floor((ind-1)/4)*4+4;
            %    new_optionsr(remain_loc,ind_p)=0; %remove inds for the piece already determined from all the other locations 
            %end
            for direction=1:4
                if direction==1 && floor((location_ind-1)/16)+1>1
                    next_location=location_ind-16;
                elseif direction==2 && mod(location_ind-1,16)+1<16
                    next_location=location_ind+1;
                elseif direction==3 && floor((location_ind-1)/16)+1<16
                    next_location=location_ind+16;
                elseif direction==4 && mod(location_ind-1,16)+1>1
                    next_location=location_ind-1;
                else
                    continue;
                end
                
                %$if next_location==43
                 %   dummy=0;
                %end

                ind_flags=false(1,256*4);
                v_t=vect_ind(prev_options(location_ind,:));
                if ~isempty(v_t) %|| ~compromize
                    for jj=1:length(tree(1,1,:))
                        indnext=tree(v_t,direction,jj)';
                        indnextn0=indnext(indnext>0);
                        tempflag=false(1,256*4);
                        tempflag(indnextn0)=true;
                        ind_flags=ind_flags | tempflag;
                    end
                    temp_v_next=new_optionsr(next_location,:) & ind_flags;
                    new_optionsr(next_location,:)=temp_v_next;
                else
                    score_nostop=0; %also, do not use the empty location to add constriant on the next location
                end
                %if isempty(v_t)
                %    score_nostop=0;
                %end
            end
       end
       imglin=sum(new_optionsr,2);
       score_solved=sum(imglin==1);

       if (sum(new_optionsr(:)~=prev_options(:))==0 || all((score_nostop==0) && ~compromize) ) && clk>3
           if score_nostop>0 || compromize || score_solved>250
               score_nostop=score_solved;
           end
           break;
       end
       if score_solved>score_solved_prev  %update tree if number of definite solved pieces increased
            score_solved_prev=score_solved;
            show_state=1*(imglin==1);
            for indlocation=1:256
                if show_state(indlocation)==1
                    show_index(indlocation)=min(indf(new_optionsr(indlocation,:)>0));
                end
            end
            %%[new_optionsr,show_state,show_index]=force_to_options(new_optionsr,options>0,show_state,show_index);
            tree=update_tree(tree,show_state,show_index);
       end
    end

end

function [new_optionsr]=correct_options(options,tree)
    bW=16; %board width
    bS=bW^2;%size of board
    b4S=bS*4;

    new_optionsr=options>0;
    imglin=sum(new_optionsr,2);

    indf=1:b4S;
    imglin=sum(new_optionsr>0,2);
    show_state=1*(imglin==1);
    show_index=zeros(size(show_state));
    for indlocation=1:bS
        if show_state(indlocation)==1
            show_index(indlocation)=min(indf(new_optionsr(indlocation,:)>0));
        end
    end
    tree=update_tree(tree,show_state,show_index); %I think already done, consider remove

    prev_options=new_optionsr;
    vect_ind=1:bS*4;
    for location_ind=1:bS
        for direction=1:4
            if direction==1 && floor((location_ind-1)/bW)+1>1
                next_location=location_ind-bW;
            elseif direction==2 && mod(location_ind-1,bW)+1<bW
                next_location=location_ind+1;
            elseif direction==3 && floor((location_ind-1)/bW)+1<bW
                next_location=location_ind+bW;
            elseif direction==4 && mod(location_ind-1,bW)+1>1
                next_location=location_ind-1;
            else
                continue;
            end

            ind_flags=false(1,bS*4);
            v_t=vect_ind(prev_options(location_ind,:));
            if ~isempty(v_t) 
                for jj=1:length(tree(1,1,:))
                    indnext=tree(v_t,direction,jj)';
                    indnextn0=indnext(indnext>0);
                    tempflag=false(1,bS*4);
                    tempflag(indnextn0)=true;
                    ind_flags=ind_flags | tempflag;
                end
                temp_v_next=new_optionsr(next_location,:) & ind_flags;
                new_optionsr(next_location,:)=temp_v_next;
            end
        end
    end


end





function new_tree=update_tree(tree,show_state,show_index)
       new_tree=tree;
       slen=length(tree(1,1,:));
       vect_ind=1:256*4;
       for location_ind=1:256
            if show_state(location_ind)>0 && show_index(location_ind)>0
                ind=show_index(location_ind);
                indpiecelist=floor((ind-1)/4)*4+1:floor((ind-1)/4)*4+4;
                indremovelist=indpiecelist(indpiecelist~=ind);
                new_tree(indremovelist,:,:)=-1;
            end
       end

       for location_ind=1:256
            %if sum(new_optionsr(location_ind,1:4:end) | new_optionsr(location_ind,2:4:end) | new_optionsr(location_ind,3:4:end) | new_optionsr(location_ind,4:4:end) )==1 %certinty in one location
            %    ind=max(vect_ind(new_optionsr(location_ind,:)));
            %    remain_loc=[1:location_ind-1 location_ind+1:256];
            %    ind_p=floor((ind-1)/4)*4+1:floor((ind-1)/4)*4+4;
            %    new_optionsr(remain_loc,ind_p)=0; %remove inds for the piece already determined from all the other locations 
            %end
            if all(show_state(location_ind)==0 ||show_index(location_ind)==0)
                continue;
            end
            here_index=show_index(location_ind);
            for direction=1:4
                if direction==1 && floor((location_ind-1)/16)+1>1
                    next_location=location_ind-16;
                elseif direction==2 && mod(location_ind-1,16)+1<16
                    next_location=location_ind+1;
                elseif direction==3 && floor((location_ind-1)/16)+1<16
                    next_location=location_ind+16;
                elseif direction==4 && mod(location_ind-1,16)+1>1
                    next_location=location_ind-1;
                else
                    new_tree(here_index,direction,:)=permute([0; -1*ones(slen-1,1)],[3 2 1]);
                    continue;
                end
                if all(show_state(next_location)==0 || show_index(next_location)==0)
                    continue;
                end
                next_index=show_index(next_location);
                indjj= 1:slen;
                indnext_v=tree(here_index,direction,indjj);
                if sum(indnext_v==next_index)==1
                    new_tree(here_index,direction,:)=permute([next_index; -1*ones(slen-1,1)],[3 2 1]);
                end
            end
       end

end


function score_nostop=WillStop(options,tree)

    new_optionsr=options>0;
    score_nostop=1;

    for clk=2:32   
       prev_options=new_optionsr;
       vect_ind=1:256*4;
       for location_ind=1:256
            for direction=1:4
                if direction==1 && floor((location_ind-1)/16)+1>1
                    next_location=location_ind-16;
                elseif direction==2 && mod(location_ind-1,16)+1<16
                    next_location=location_ind+1;
                elseif direction==3 && floor((location_ind-1)/16)+1<16
                    next_location=location_ind+16;
                elseif direction==4 && mod(location_ind-1,16)+1>1
                    next_location=location_ind-1;
                else
                    continue;
                end
                
                ind_flags=false(1,256*4);
                v_t=vect_ind(prev_options(location_ind,:));
                if ~isempty(v_t)
                    for jj=1:length(tree(1,1,:))
                        indnext=tree(v_t,direction,jj)';
                        indnextn0=indnext(indnext>0);
                        tempflag=false(1,256*4);
                        tempflag(indnextn0)=true;
                        ind_flags=ind_flags | tempflag;
                    end
                    temp_v_next=new_optionsr(next_location,:) & ind_flags;
                    new_optionsr(next_location,:)=temp_v_next;
                else
                    score_nostop=0; %also, do not use the empty location to add constriant on the next location
                end
            end
       end
       if sum(new_optionsr(:)~=prev_options(:))==0 || score_nostop==0
           if score_nostop>0
               imglin=sum(new_optionsr,2);
               score_nostop=sum(imglin==1);
           end
           break;
       end
    end
end



function test()

    for oo=1:256
    for pp=oo+1:256
    if sum(options(oo,:)==options(pp,:))==1024
    disp(sprintf('oo=%d  pp=%d  deg=%d',oo,pp,sum(options(oo,:)>0)));
    end
    end
    end
    
    indf=1:1024;
    indf(options(123,:)>0)
    
    indt=1:256;
    imglin=sum(options,2);
    
    RRR= orth(single(options'));
    RRR(abs(RRR)<0.1)=0;

    RRR=linsolve(single(options),ones(256,1));
    sum(RRR)
end



function show_original(h,fig,original,testno)
    return

end


% Capture clicks on uiaxes
% Capture key presses


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%    ANALYZE    %%%%%%%%%%%
function [score1,options]=analyze(level,tree,treesize,options)
   compromize=false;
   indf=1:1024;
   last_options=options;
   backup_tree=tree;

   for level=0:level     
    

        %Perturbations
        %perturb_vectors
        perturb_vectors=zeros(256*4,1000,256);
        perturb_param=zeros(256*4,1000,5);
        perturb_score=zeros(256*4,1000);
        prev_min_score=zeros(256*4,1);
        
        if level>=1

    
            parfor ind=1:256*4
                small0big1=1;
                [temp_score,temp_index_vectors,temp_perturb_vectors,temp_perturb_param]=learn(ind,tree,treesize,options,small0big1);
                score(ind,:,:)=temp_score;
                min_score(ind)=min(temp_score(temp_score(:)>0));
                index_vectors(ind)=temp_index_vectors;
                perturb_vectors(ind,:,:)=temp_perturb_vectors;
                perturb_param(ind,:,:)=temp_perturb_param;
            end
            
            tempv=1:256*4;
            candidates=length(tempv(min_score==min(min_score)));
            
            ind=min(tempv(min_score==min(min_score)));
            tempt=1:1000;
            t=min(tempt(score(ind,:,1)==min(min_score)));
            
            if perturb_param(ind,t,2)==-1
                small0big1=1;
            else
                small0big1=0;
            end
            
            if small0big1==0
                for direction= perturb_param(ind,t,2)
                    for jj= perturb_param(ind,t,4)
                        ind2=tree(ind,direction,jj);
                        if ind2>0
                            %I:  small supressor
                            tree(ind,direction,jj)=-1; %temporarily remove option
                            direction2=mod(direction-1+2,4)+1;
                            for jj2=1:treesize(ind2,direction2)
                                if tree(ind2,direction2,jj2)==ind
                                    tree(ind2,direction2,jj2)=-1;
                                    break;
                                end
                            end
            
            
                        end
            
                    end
                end
            end %if small
            
            if small0big1==1
            
                %large promoter
                for ind4=[floor((ind-1)/4)*4+1:ind-1 ind+1:floor((ind-1)/4)*4+4]
                    for direction=1:4
                        for jj=1:treesize(ind4,direction)
                            ind2=tree(ind4,direction,jj);
                            if ind2>0
                                tree(ind4,direction,jj)=-1; %temporarily remove option
                                direction2=mod(direction-1+2,4)+1;
                                for jj2=1:treesize(ind2,direction2)
                                    if tree(ind2,direction2,jj2)==ind4
                                        tree(ind2,direction2,jj2)=-1;
                                        break;
                                    end
                                end
                            end
                        end
                    end
                end
                
            
            end %if large
            

        end %if level>=1 
    
        [options,score_nostop]=new_options(compromize,options,tree);
		if score_nostop==0
            disp('break predicted')
        end


        imglin=sum(options>0,2);
        show_state=1*(imglin==1);
        show_index=zeros(size(show_state));
        for indlocation=1:256
            if show_state(indlocation)==1
                show_index(indlocation)=min(indf(options(indlocation,:)>0));
            else
                show_index(indlocation)=0;
            end
        end

        [options,show_state,show_index]=force_to_options(options,options,show_state,show_index);

    end %for level

    score1=sum(sum(options,2)); %to be minimized
    %tree=backup_tree;
    %disp(sprintf('exp(Entropy) =%d',score1))
    %delete(gcp('nocreate'));
end


%ASSOCIATED with analyze
function [temp_score,index_vectors,perturb_vectors,perturb_param]=learn(ind,tree,treesize,options,small0big1)
   compromize=false;

    index_vectors=0;
    backup_tree=tree;
    if small0big1==0
        for direction=1:4
            for jj=1:treesize(ind,direction)
                ind2=tree(ind,direction,jj);
                if ind2>0
                    %I:  small supressor
                    tree(ind,direction,jj)=-1; %temporarily remove option
                    direction2=mod(direction-1+2,4)+1;
                    for jj2=1:treesize(ind2,direction2)
                        if tree(ind2,direction2,jj2)==ind
                            tree(ind2,direction2,jj2)=-1;
                            break;
                        end
                    end
                    
                    [ttemp_options,score_nostop]=new_options(compromize,options,tree);
                    imglin=sum(ttemp_options,2);
                    

                    if true%sum(abs(imglin))>0
                        index_vectors=index_vectors+1;
                        perturb_vectors(index_vectors,:)=imglin';
                        perturb_param(index_vectors,:)=[ind direction ind2 jj -1]; %small supressor
                        temp_score(index_vectors,1)=(min(imglin)==0)*256000+sum(imglin);
                        %figure(1)    
                        %imshow(reshape(abs(imglin),[16 16])*256/max(abs(imglin)));
                    end
                    tree=backup_tree;
    
                    %II:  small promoter
                    if false
                        for jjj=[1:jj-1 jj+1:treesize(ind,direction)]
                            ind3=tree(ind,direction,jjj);
                            if ind3>0
                                tree(ind,direction,jjj)=-1; %temporarily remove option
                                direction2=mod(direction-1+2,4)+1;
                                for jj2=1:treesize(ind3,direction2)
                                    if tree(ind3,direction2,jj2)==ind
                                        tree(ind3,direction2,jj2)=-1;
                                        break;
                                    end
                                end
                            end
                        end
                        [ttemp_options,score_nostop]=new_options(compromize,options,tree);
                        imglin=sum(ttemp_options,2);
        
                        if true% sum(abs(imglin))>0
                            index_vectors=index_vectors+1;
                            perturb_vectors(index_vectors,:)=imglin';
                            perturb_param(index_vectors,:)=[ind direction ind2 jj +1]; %small promoter
                            temp_score(index_vectors,1)=(min(imglin)==0)*256000+sum(imglin);
                            %figure(1)    
                            %imshow(reshape(abs(imglin),[16 16])*256/max(abs(imglin)));
                        end
                        tree=backup_tree;
                    end
                end
            
            end
        end
    end %if small

    if small0big1==1
        %large supressor
        if false
            for direction=1:4
                for jj=1:treesize(ind,direction)
                    ind2=tree(ind,direction,jj);
                    if ind2>0
                        tree(ind,direction,jj)=-1; %temporarily remove option
                        direction2=mod(direction-1+2,4)+1;
                        for jj2=1:treesize(ind2,direction2)
                            if tree(ind2,direction2,jj2)==ind
                                tree(ind2,direction2,jj2)=-1;
                                break;
                            end
                        end
                    end
                end
            end
            [ttemp_options,score_nostop]=new_options(compromize,options,tree);
            imglin=sum(ttemp_options,2);
        
            if true%sum(abs(imglin))>0
                index_vectors=index_vectors+1;
                perturb_vectors(index_vectors,:)=imglin';
                perturb_param(index_vectors,:)=[ind -1 -1 -1 -1]; %large supressor
                temp_score(index_vectors,1)=(min(imglin)==0)*256000+sum(imglin);
                %figure(1)
                %imshow(reshape(abs(imglin),[16 16])*256/max(abs(imglin)));
            end
            tree=backup_tree;
        end

        %large promoter
        for ind4=[floor((ind-1)/4)*4+1:ind-1 ind+1:floor((ind-1)/4)*4+4]
            for direction=1:4
                for jj=1:treesize(ind4,direction)
                    ind2=tree(ind4,direction,jj);
                    if ind2>0
                        tree(ind4,direction,jj)=-1; %temporarily remove option
                        direction2=mod(direction-1+2,4)+1;
                        for jj2=1:treesize(ind2,direction2)
                            if tree(ind2,direction2,jj2)==ind4
                                tree(ind2,direction2,jj2)=-1;
                                break;
                            end
                        end
                    end
                end
            end
        end
        [ttemp_options,score_nostop]=new_options(compromize,options,tree);
        imglin=sum(ttemp_options,2);
    
        if true%sum(abs(imglin))>0
            index_vectors=index_vectors+1;
            perturb_vectors(index_vectors,:)=imglin';
            perturb_param(index_vectors,:)=[ind -1 -1 -1 +1]; %large promoter
            temp_score(index_vectors,1)=(min(imglin)==0)*256000+sum(imglin);
            %figure(1)
            %imshow(reshape(abs(imglin),[16 16])*256/max(abs(imglin)));
        end

    end %if large

    if index_vectors>1000
        perturb_vectors=perturb_vectors(1:1000,:);
        perturb_param=perturb_param(1:1000,:);
        temp_score=temp_score(1:1000,1);
    elseif index_vectors<1000 && index_vectors>0
        perturb_vectors=[perturb_vectors ; zeros(1000-index_vectors,256)];
        perturb_param=[perturb_param ; zeros(1000-index_vectors,5)];
        temp_score=[temp_score; zeros(1000-index_vectors,1)];
        if max(temp_score)==0
            temp_score(1,1)=256000;
        end
    else
        perturb_vectors= zeros(1000,256);
        perturb_param=zeros(1000-index_vectors,5);
        temp_score= zeros(1000-index_vectors,1);
        temp_score(1,1)=256000;
    end
    


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter for probablistic generated options %%%

function options=prob_filter(options,baseline_options,treesize,tree)

    
    next_options=options;
    if true
        for st=1:256
            tempu=1:1024;
            tmp=max(tempu(options(st,:)==max(options(st,:))));
            tempu_el=tempu(tempu~=tmp);
            next_options(st,tempu_el)=0;
        end
        for st=1:4:1024
            tempv=(1:256)'*[1 1 1 1];
            tempu=ones(256,1)*[1 2 3 4];
            vp=tempv(options(:,st:st+3)==max(max(options(:,st:st+3))));
            up=tempu(options(:,st:st+3)==max(max(options(:,st:st+3))));
            f=min(length(up),length(vp));
            fp=1;%floor(rand*f)+1;
            u=up(fp);
            v=vp(fp);
            tmplt=zeros(256,4);
            tmplt(v,u)=1;
            next_options(:,st:st+3)=next_options(:,st:st+3).*tmplt; %eliminate competing locations for the piece
        end
    
    end
    
    if true
        indf=1:1024;
        score_possible=zeros(1,256);
        for row=1:16
            for col=1:16
                next_imglin=sum(next_options>0.5,2);
                ind_location=(col-1)+(row-1)*16+1;
                if next_imglin(ind_location)==1
                    current_index=indf(next_options(ind_location,:)>0);
                    index_possible=current_index;
                    count=0;
                    for dir=1:4
                        flag=false;
                        if dir==1
                            if row>1
                                if next_imglin(ind_location-16)==1
                                    next_index=indf(next_options(ind_location-16,:)>0);
                                    for t=1:treesize(current_index,dir)
                                        if tree(current_index,dir,t)==next_index
                                            flag=true;
                                            break;
                                        end
                                    end
                                else
                                    flag=true;
                                end
                            else
                                flag=true;
                            end
                        end
                        if dir==2
                            if col<16
                                if next_imglin(ind_location+1)==1
                                    next_index=indf(next_options(ind_location+1,:)>0);
                                    for t=1:treesize(current_index,dir)
                                        if tree(current_index,dir,t)==next_index
                                            flag=true;
                                            break;
                                        end
                                    end
                                else
                                    flag=true;
                                end
                            else
                                flag=true;
                            end
                        end
    
                        if dir==3
                            if row<16
                                if next_imglin(ind_location+16)==1
                                    next_index=indf(next_options(ind_location+16,:)>0);
                                    for t=1:treesize(current_index,dir)
                                        if tree(current_index,dir,t)==next_index
                                            flag=true;
                                            break;
                                        end
                                    end
                                else
                                    flag=true;
                                end
                            else
                                flag=true;
                            end
                        end
                        if dir==4
                            if col>1
                                if next_imglin(ind_location-1)==1
                                    next_index=indf(next_options(ind_location-1,:)>0);
                                    for t=1:treesize(current_index,dir)
                                        if tree(current_index,dir,t)==next_index
                                            flag=true;
                                            break;
                                        end
                                    end
                                else
                                    flag=true;
                                end
                            else
                                flag=true;
                            end
                        end
                        if flag
                            count=count+1;
                        end
    
                    end
                    if count==4
                        score_possible(ind_location)=1;
                    else
                        score_possible(ind_location)=0;
                    end
    
                end
            end
        end
    
    
        for ind_location=1:256
            if score_possible(ind_location)==0
                options(ind_location,:)=baseline_options(ind_location,:);
            end
        end
    
    end


end






%Assistance for probablistic approach

function [score1,options]=analyze_prob(levelmax,tree,treesize,options)
    score1=0;
    compromize=true;

    backup_tree=tree;
    
    %startSafeParpool(parallelForce)
    

    for level=1:levelmax     
    
        ind_vector=1:16*16*4;      
        
        
        %Perturbations
        %perturb_vectors
        perturb_vectors=zeros(256*4,1000,256);
        perturb_param=zeros(256*4,1000,5);
        score=zeros(256*4,1000,1);
        prev_min_score=zeros(256*4,1);
        grand_tree=-1*ones(256*4,1024,4,57);
        
        parfor ind=1:256*4  %%%%%%%%%%%%%    parfor
            temp_tree=tree;
            temp_score=1000000*ones(1000,1);
            temp_index_vectors=zeros(1000,256);
            temp_perturb_vectors=zeros(1000,256);
            temp_perturb_param=zeros(1000,5);
            
            small0big1=1;
            [temp_tree,temp_score,temp_index_vectors,temp_perturb_vectors,temp_perturb_param]=learn_prob(ind,tree,treesize,options,small0big1,compromize);
            grand_tree(ind,:,:,:)=temp_tree;
            score(ind,:,:)=temp_score;
            min_score(ind)=min(temp_score(temp_score(:)>0));
            index_vectors(ind)=temp_index_vectors;
            perturb_vectors(ind,:,:)=temp_perturb_vectors;
            perturb_param(ind,:,:)=temp_perturb_param;
        end
        
        tempv=1:256*4;
        candidates=length(tempv(min_score==min(min_score)));
        
        ind=min(tempv(min_score==min(min_score)));
        tempt=1:1000;
        t=min(tempt(score(ind,:,1)==min(min_score)));
        
        
        if perturb_param(ind,t,2)==-1 && min(min_score)<10000
            small0big1=-2;
            tree=permute(grand_tree(ind,:,:,:),[2 3 4 1]);
        elseif perturb_param(ind,t,2)>0 && min(min_score)<10000
            disp('unexpected');
            small0big1=0;
        else
            small0big1=2;
            disp('condition not known');
        end 
       
        if small0big1==0 
            for direction= perturb_param(ind,t,2)
                for jj= perturb_param(ind,t,4)
                    ind2=tree(ind,direction,jj);
                    if ind2>0
                        %I:  small supressor
                        tree(ind,direction,jj)=-1; %temporarily remove option
                        direction2=mod(direction-1+2,4)+1;
                        for jj2=1:treesize(ind2,direction2)
                            if tree(ind2,direction2,jj2)==ind
                                tree(ind2,direction2,jj2)=-1;
                                break;
                            end
                        end
        
                        %imglin=sum(new_options(compromize,options,tree),2);
                    end
        
                end
            end
        end %if small
        
        if small0big1==1
        
            %large promoter
            for ind4=[floor((ind-1)/4)*4+1:ind-1 ind+1:floor((ind-1)/4)*4+4]
                for direction=1:4
                    for jj=1:treesize(ind4,direction)
                        ind2=tree(ind4,direction,jj);
                        if ind2>0
                            tree(ind4,direction,jj)=-1; %temporarily remove option
                            direction2=mod(direction-1+2,4)+1;
                            for jj2=1:treesize(ind2,direction2)
                                if tree(ind2,direction2,jj2)==ind4
                                    tree(ind2,direction2,jj2)=-1;
                                    break;
                                end
                            end
                        end
                    end
                end
            end
            %imglin=sum(new_options(compromize,options,tree),2);
        
        end %if large
        
        
        locind_full=1:256;
        locind_vector=locind_full(options(:,ind)>0);
        locind=perturb_param(ind,t,5);
        if length(locind_vector)>1 && locind>0
            locind_el=locind_vector(locind_vector~=locind);
            options(locind_el,ind)=0;
        end

        [options,score]=new_options_prob(compromize,options,tree);
        
        %imglin=sum(options>0.5,2);
        %score_sum=sum(imglin);

        %predesigned puzzle
        %required_options=zeros(256,1024);
        %for tt=1:256
        %    required_options(tt,(tt-1)*4+1)=1;
        %end
        %actual_score=sum(sum((options>0).*required_options));
        
        Rrank=rank(double(options));
        if Rrank==256
            RRR=linsolve(double(options),ones(256,1));
            added_points=sum(RRR);
        else
            RRR=lsqminnorm(double(options),ones(256,1));
            added_points=sum(RRR);
        end
        score1=added_points;
    
    end %for level

    tree=backup_tree;
    disp(sprintf('Solution sum =%d',score1))
    %delete(gcp('nocreate'));


end %function



function [temp_tree,temp_score,index_vectors,perturb_vectors,perturb_param]=learn_prob(ind,tree,treesize,options,small0big1,compromize)
    index_vectors=0;
    backup_tree=tree;
    required_options=zeros(256,1024);
    for t=1:256
       required_options(t,(t-1)*4+1)=1;
    end
    temp_tree=tree;
    
    perturb_vectors=zeros(1000,256);
    perturb_param=zeros(1000,5);
    temp_score=zeros(1000,1);
    
    if small0big1==0
        for direction=1:4
            for jj=1:treesize(ind,direction)
                ind2=tree(ind,direction,jj);
                if ind2>0
                    %I:  small supressor
                    tree(ind,direction,jj)=-1; %temporarily remove option
                    direction2=mod(direction-1+2,4)+1;
                    for jj2=1:treesize(ind2,direction2)
                        if tree(ind2,direction2,jj2)==ind
                            tree(ind2,direction2,jj2)=-1;
                            break;
                        end
                    end
                    
                    [temp_options,~]=new_options_prob(compromize,options,tree);
                    imglin=sum(temp_options,2);
                    %actual_score=sum(sum((temp_options>0).*required_options));
                   

                    if true %actual_score>=256%rank(single(temp_options))>=Rrank-1
                        index_vectors=index_vectors+1;
                        perturb_vectors(index_vectors,:)=imglin';
                        perturb_param(index_vectors,:)=[ind direction ind2 jj -1]; %small supressor
                        extra=0;
                        %if Rrank==256
                        %   RRR=linsolve(single(temp_options),ones(256,1));
                        %else
                        RRR=lsqminnorm(double(temp_options),ones(256,1));
                        %end
                        extra=max(sum(RRR),0);



                        temp_score(index_vectors,1)=(min(imglin)==0)*256000+abs(256-extra)+1;
                        %figure(1)    
                        %imshow(reshape(abs(imglin),[16 16])*256/max(abs(imglin)));
                    end
                    tree=backup_tree;
    
                end
            
            end
        end
    end %if small

    if small0big1==1

        %large promoter
        for ind4=[floor((ind-1)/4)*4+1:ind-1 ind+1:floor((ind-1)/4)*4+4]
            for direction=1:4
                for jj=1:treesize(ind4,direction)
                    ind2=tree(ind4,direction,jj);
                    if ind2>0
                        tree(ind4,direction,jj)=-1; %temporarily remove option
                        direction2=mod(direction-1+2,4)+1;
                        for jj2=1:treesize(ind2,direction2)
                            if tree(ind2,direction2,jj2)==ind4
                                tree(ind2,direction2,jj2)=-1;
                                break;
                            end
                        end
                    end
                end
            end
        end
        [temp_options,~]=new_options_prob(compromize,options,tree);
        imglin=sum(temp_options,2);
        %actual_score=sum(sum((temp_options>0).*required_options));
    
        if true %actual_score>=256 %rank(single(temp_options))>=Rrank-1
            temp_tree=tree;
            index_vectors=index_vectors+1;
            perturb_vectors(index_vectors,:)=imglin';
            perturb_param(index_vectors,:)=[ind -1 -1 -1 +1]; %large promoter
            extra=0;
            %if Rrank==256
            %     RRR=linsolve(single(temp_options),ones(256,1));
            %else
            RRR=lsqminnorm(double(temp_options),ones(256,1));
            %end
            extra=max(sum(RRR),0);

            temp_score(index_vectors,1)=(min(imglin)==0)*256000+abs(256-extra)+1;
            %figure(1)
            %imshow(reshape(abs(imglin),[16 16])*256/max(abs(imglin)));
        end

    end %if large

    if index_vectors>1000
        perturb_vectors=perturb_vectors(1:1000,:);
        perturb_param=perturb_param(1:1000,:);
        temp_score=temp_score(1:1000,1);
    elseif index_vectors<1000 && index_vectors>0
        if max(temp_score)==0
            temp_score(1,1)=256000;
        end
    else
        temp_score(1,1)=256000;
    end
    


end



function [new_options,score]=new_options_prob(compromize,options,tree)

    score=0;
   new_options=double(options);
    jjvect=1:length(tree(1,1,:));
    maxclk=40;
    for clk=2:maxclk   
       vect_ind=1:256*4;
       next_new_options=double(ones(256,1024));
       for location_ind=1:256
            
            %if sum(new_options(location_ind,1:4:end)>0 | new_options(location_ind,2:4:end)>0 | new_options(location_ind,3:4:end)>0 | new_options(location_ind,4:4:end)>0 )==1 %certinty in one location
            %    ind=max(vect_ind(new_options(location_ind,:)>0));
            %    remain_loc=[1:location_ind-1 location_ind+1:256];
            %    ind_p=floor((ind-1)/4)*4+1:floor((ind-1)/4)*4+4;
            %    new_options(remain_loc,ind_p)=0; %remove inds for the piece already determined from all the other locations 
            %    sbuf=sum(new_options(remain_loc,:));
            %    new_options(remain_loc,:)=new_options(remain_loc,:)./(sbuf+(sbuf==0));
            %end
            for direction=1:4
                if direction==1 && floor((location_ind-1)/16)+1>1
                    next_location=location_ind-16;
                elseif direction==2 && mod(location_ind-1,16)+1<16
                    next_location=location_ind+1;
                elseif direction==3 && floor((location_ind-1)/16)+1<16
                    next_location=location_ind+16;
                elseif direction==4 && mod(location_ind-1,16)+1>1
                    next_location=location_ind-1;
                else
                    continue;
                end
                
                ind_prob=double(zeros(1,256*4));
                actvect=vect_ind(new_options(location_ind,:)>0); %0.000001
                for tt=1:length(actvect)
                    indnext=tree(actvect(tt),direction,jjvect);
                    indnextn0=indnext(indnext>0); %indices connected by index actvect(tt)
                    count0=double(length(indnextn0));
                    if count0>0
                        ind_prob(indnextn0)=ind_prob(indnextn0)+new_options(location_ind,actvect(tt))/count0; %add probability according to number of coneected in tree and according to probablity of the local index
                        %ind_prob(indnextn0)=ind_prob(indnextn0)+1/count0; %add probability according to number of coneected in tree and according to probablity of the local index
                    end
                end
                next_new_options(next_location,:)=next_new_options(next_location,:).*ind_prob;
            end

       end %for locations
       sbuff=sum(next_new_options,2);
       next_new_options=next_new_options./(sum(next_new_options,2)+(sbuff==0));
       if  min(next_new_options(:))<1e-100 || clk==maxclk %sum(sum(abs(new_options-next_new_options)<1e-12))==1024*256 ||
            new_options=next_new_options;
            RRR=lsqminnorm(double(new_options),ones(256,1));
            extra=max(sum(RRR),0);
            imglin=sum(new_options>0,2);
            score=256*256/((min(imglin)==0 && compromize==0)*256000+abs(256-extra)+1);

           break;
       else
           new_options=next_new_options;
       end
    end


end


%####################################################
%#####      Consult AI    ###########################

function [tree,treesize]=consult_AI(piece,tree,treesize,e2folder,planv)
    global pp_choice;
    global pp_avail;
    global pp_or;
    global pp_list;
    global maxn2;
    pt0=1:257;
    dir0=1:4;
    vectavail0=1:64;

    pp_list0=1:256;
    maxn2=0;
    pp_choice=zeros(256,257);
    pp_avail=zeros(256,257,64);
    pp_or=zeros(257,1);%chosen orientation, to be consistent
    pp_or(257)=1;
    %prepare avail matrix based on tree
    for piece1=1:256
        for ttdir=1:4
            for ttor1=1:4
                ind1=(piece1-1)*4+ttor1;
                for temp_p=1:treesize(ind1,ttdir)
                    ind2=tree(ind1,ttdir,temp_p);
                    if ind2==0
                        piece2=257;
                        ttor2=1;
                    elseif ind2>0
                        piece2=floor((ind2-1)/4)+1;
                        ttor2=mod(ind2-1,4)+1;
                    else
                        continue;
                    end
                    pos=(ttdir-1)*16+(ttor1-1)*4+ttor2;
                    pp_avail(piece1,piece2,pos)=1;
                end
            end
        end
    end
   
    if planv==0
        only_orientation_1=true;
    else
        only_orientation_1=false;
    end
    tic;
    keep_pp_or=zeros(257,1);
    pp_or=ones(257,1);
    useway=2;
    if useway==1
        for testcount=1:1000000
            pp_choice=zeros(256,257);
            pp_or=zeros(257,1);%chosen orientation, to be consistent
            pp_or(257)=1;
            pp_list=pp_list0(randperm(length(pp_list0)));
            check=recurs_fillchoice(1,2);  %choose only one consistent option for pairing between every two pieces (randomly if needed)
            if check
                disp("Succeeded");
                break
            else
                if sum(pp_or>0)==257
                    keep_pp_choice=pp_choice;
                    keep_pp_or=pp_or;
                    keep_pp_list=pp_list;
                end
                if toc>4*60*60 && sum(keep_pp_or>0)==257
                    pp_or=keep_pp_or;
                    pp_list=keep_pp_list;
                    pp_choice=keep_pp_choice;
                end
                if toc>4*60*60
                    disp("Ended as a compromize due to length of time")
                    break;
                end
            end
        end
        disp(sprintf('Max n2= %g',maxn2));
        if only_orientation_1
            disp(sprintf('Number of correct orientations= %g',sum(pp_or==1)));
        end
        %Prepare matrix to be saved as csv
        features=zeros(256,4);
        edge_pairs_mat=zeros(256,256);
        for tp1=1:256
            tor1=1;
            temp=[0 0 0 0];
            for tp2=1:256
                tval=pp_choice(tp1,tp2);
                if tval>0 && tval<=64
                    tdir=floor((tval-1)/16)+1;
                    edge_pairs_mat(tp1,tp2)=tdir;
                    if sum(temp)==0
                        tor1=floor((tval-1-(tdir-1)*16)/4)+1;
                        tor2=mod(tval-1,4)+1;
                        for tdir=1:4
                            temp(tdir)=piece(tp1,mod(tdir-1-(tor1-1),4)+1);
                        end
                    end
                end
            end
            features(tp1,:)=temp;
        end
    elseif useway==2

        %Prepare matrix to be saved as csv
        features=zeros(256,4);
        edge_pairs_mat=zeros(256,256);
        pp_or=floor(4*rand(256,1))+1;
        for tp1=1:256
            temp=[0 0 0 0];
            for tdir=1:4
                temp(tdir)=piece(tp1,mod(tdir-1-(pp_or(tp1)-1),4)+1);
            end
            features(tp1,:)=temp;
            tind1=(tp1-1)*4+pp_or(tp1);
            for tdir=1:4
                for tpos=1:treesize(tind1,tdir)
                    tind2=tree(tind1,tdir,tpos);
                    if tind2>0
                        tor2=mod((tind2-1),4)+1;
                        tp2=floor((tind2-1)/4)+1;
                        if tor2==pp_or(tp2)
                            edge_pairs_mat(tp1,tp2)=tdir;
                        end
                    end
                end
            end
        end
                
   
    end %if useway==1
    disp(sprintf('Number of connections= %g',sum(edge_pairs_mat(:))));
    if sum(edge_pairs_mat(:))>5000
        disp('thus sending to NN inference.');
        writematrix(features, sprintf('%s//node_features.csv',e2folder));
        writematrix(edge_pairs_mat, sprintf('%s//edge_type_matrix.csv',e2folder));
        %Run python inference code
        commandext=[OceanFolder '/venv/bin/python ' OceanFolder '/RemoveFakeConnections.py'];
        system(commandext);
        pause(3);
        clean_edge_pair_mat=readmatrix([e2folder '/clean_edge_type_matrix.csv']);
    else
        writematrix(features, sprintf('%s//node_features.csv',e2folder));
        writematrix(edge_pairs_mat, sprintf('%s//edge_type_matrix.csv',e2folder));
        clean_edge_pair_mat=edge_pairs_mat;
    end %if sum(edge_pairs_mat(:))>5000

    old_tree=tree;
    old_treesize=treesize;   %will copy from previous tree only parts that do not include the picked orientations
    tree=-1*ones(256*4,4,200); %was 57 % from which object according to position (1 to 256*4) to all possible objects(1 of 256*4) up to 128 (actually there are 49 max),  according to conntection direction (1 of 4).
    treesize=zeros(size(treesize));
    for ind=1:256*4
        for dirct=1:4
            if old_tree(ind,dirct,1)==0 
                tree(ind,dirct,1)=0;   %copy all poisbble connections to margin (0)
                treesize(ind,dirct)=1;
            end
            for pos_t=treesize(ind,dirct)+1:old_treesize(ind,dirct)
                if old_tree(ind,dirct,pos_t)>0
                    tempind=old_tree(ind,dirct,pos_t);
                    if mod(tempind-1,4)+1~=pp_or(floor((tempind-1)/4)+1) || mod(ind-1,4)+1~=pp_or(floor((ind-1)/4)+1)
                        treesize(ind,dirct)=treesize(ind,dirct)+1;
                        tree(ind,dirct,treesize(ind,dirct))=old_tree(ind,dirct,pos_t);
                    end
                end
            end
        end
    end
    for piece1=1:256
        for piece2=1:256
            if piece1==piece2 || clean_edge_pair_mat(piece1,piece2)==0 || pp_or(piece1)==0 || pp_or(piece2)==0
                continue;
            end
            dirct=clean_edge_pair_mat(piece1,piece2);
            tor1=pp_or(piece1);
            tor2=pp_or(piece2);
            ind1=(piece1-1)*4+tor1;
            ind2=(piece2-1)*4+tor2;
            treesize(ind1,dirct)=treesize(ind1,dirct)+1;
            tree(ind1,dirct,treesize(ind1,dirct))=ind2;
        end
    end

    return;


    %recursive function
    function check=recurs_fillchoice(n1,n2)
        p1=pp_list(n1);
        p2=pp_list(n2);
        vectavail=vectavail0(pp_avail(p1,p2,:)==1);
        %identify avaiable slots to choose from
        if pp_choice(p1,p2)>0
            vectavail=pp_choice(p1,p2);
        elseif pp_choice(p1,p2)==-1
            vectavail=[];
        end
        if length(vectavail)>0
            vectavail=vectavail(randperm(length(vectavail)));
            for val=vectavail
                dir=floor((val-1)/16)+1;
                or1=floor((val-1-(dir-1)*16)/4)+1;
                or2=mod(val-1,4)+1;
                if (pp_or(p1)>0 && pp_or(p1)~=or1) || (pp_or(p2)>0 && pp_or(p2)~=or2) %conditions to remove slots
                    vectavail=vectavail(vectavail~=val);
                end
            end
        end
        flag_worked=false; %It works if p1 has connections in all directions, not particularly between p1 and p2, but here we check that one could exist
        flag_set_ppor1=false;
        flag_set_ppor2=false;

        test_pos_ind0=length(vectavail);
        if only_orientation_1
            orrandom=dir0; %if for virtual puzzle use dir0 without shuffeling for faster result 
        else
            orrandom=dir0(randperm(length(dir0))); 
        end

        if pp_choice(p1,p2)<=0
            %vectavail=[-1 vectavail]; %-1=Check that it is also valid that the couple is not connected (erase all possible connections)
            vectavail=[vectavail -1]; %-1=Check that it is also valid that the couple is not connected (erase all possible connections)
            if pp_or(p1)==0 && ~only_orientation_1
                vectavail=[vectavail -1 -1 -1 ]; 
            end
        end
        for test_pos_ind=1:length(vectavail)
            test_pos=vectavail(test_pos_ind);
            reslt=[0 0 0 0];
            pp_choice(p1,p2)=test_pos;
            
            %check with all choice matrix that it is valid
            if test_pos>0  
                p1_chosen_dir=floor((test_pos-1)/16)+1;
                p1_chosen_or1=floor((test_pos-1-(p1_chosen_dir-1)*16)/4)+1;
                reslt(p1_chosen_dir)=1;
            else
                p1_chosen_or1=pp_or(p1);
                if p1_chosen_or1==0
                    p1_chosen_or1=orrandom(test_pos_ind-test_pos_ind0);
                end
            end %if test_pos>0

            ptrandom=pt0(randperm(length(pt0)));
            for pt=ptrandom
                if pt==p1 || pt==p2
                    continue;
                end
                if pp_choice(p1,pt)==0 && pp_or(pt)==0
                    dirrandom=dir0(randperm(length(dir0)));
                    for dir=dirrandom
                        if reslt(dir)==0 && sum(pp_avail(p1,pt,(dir-1)*16+(p1_chosen_or1-1)*4+1:(dir-1)*16+(p1_chosen_or1-1)*4+4))>0
                            reslt(dir)=1;
                            if pt<257
                                break; %mutiple directions for one pair is not a success, unless it is paired with frame= piece257
                            end
                        end
                    end
                elseif pp_choice(p1,pt)==0 && pp_or(pt)>0
                    dirrandom=dir0(randperm(length(dir0)));
                    for dir=dirrandom
                        if reslt(dir)==0 && pp_avail(p1,pt,(dir-1)*16+(p1_chosen_or1-1)*4+pp_or(pt))>0
                            reslt(dir)=1;
                            if pt<257
                                break; %mutiple directions for one pair is not a success
                            end
                        end
                    end
                elseif pp_choice(p1,pt)>0
                    if pp_avail(p1,pt,pp_choice(p1,pt))>0
                        reslt(floor((pp_choice(p1,pt)-1)/16)+1)=1;
                    end
                end
                if sum(reslt)==4
                    break;
                end
            end  %for pt=ptrandom
            %if success then keep the current choice
            if sum(reslt)<4 %piece p1 does not connect at least with one piece or frame in all directions, look for different choice
                continue;
            end
            if pp_or(p1)~=p1_chosen_or1
                pp_or(p1)=p1_chosen_or1; %will be used when asessing search only relevant slots
                flag_set_ppor1=true;
            end
            if n1<n2
                if pp_choice(p1,p2)==-1
                    pp_choice(p2,p1)=-1;  %if p1 breaks from p2, p2 must the same
                end
                if pp_choice(p1,p2)>0
                    if pp_or(p1)~=p1_chosen_or1
                        pp_or(p1)=p1_chosen_or1; %will be used when asessing search only relevant slots
                        flag_set_ppor1=true;
                    end
                    val=pp_choice(p1,p2);
                    dir=floor((val-1)/16)+1;
                    or1=floor((val-1-(dir-1)*16)/4)+1;
                    or2=mod(val-1,4)+1;
                    if pp_or(p1)~=or1
                        disp('mismatch');
                        %pp_or(p1)=or1; %will be used when asessing search only relevant slots
                        %flag_set_ppor1=true;
                    end
                    if pp_or(p2)~=or2
                        pp_or(p2)=or2; %will be overwritten only by alternative slots of this layer, which were already identified
                        flag_set_ppor2=true;
                    end
                    conj_or1=or2;
                    conj_or2=or1;
                    conj_dir=mod((dir-1)+2,4)+1;
                    conj_pos=(conj_dir-1)*16+(conj_or1-1)*4+conj_or2;
                    pp_choice(p2,p1)=conj_pos;  %fill-out the conjugate choice, the back connection between piece2 and piece1
                end
            
                if n1+1<n2
                    next_n1=n1+1;
                    next_n2=n2;
                else
                    next_n1=1;
                    next_n2=n2+1;
                    if next_n2>maxn2
                        maxn2=next_n2;
                    end
                end
                if toc>4*60*60
                    flag_worked=false;
                    break;
                else
                    flag_worked=recurs_fillchoice(n2,n1); %check that p2 is left with connections as well
                end
                if flag_worked %check next couple for chosing one connection
                    if next_n2>256
                        flag_worked=true;
                        break;
                    end
                    flag_worked=recurs_fillchoice(next_n1,next_n2);
                end
            else
                flag_worked=true;
            end %if n1<n2
            if flag_worked
                break;
            else
                pp_choice(p1,p2)=0;
                pp_choice(p2,p1)=0;
                if flag_set_ppor1   %If also the pp_or state was set at this level then erase it
                    pp_or(p1)=0;
                    flag_set_ppor1=false;
                end
                if flag_set_ppor2
                    pp_or(p2)=0;
                    flag_set_ppor2=false;
                end
            end
        end %for test_pos=vectavail
        if ~flag_worked
            pp_choice(p1,p2)=0;
            if flag_set_ppor1   %If also the pp_or state was set at this level then erase it
                pp_or(p1)=0;
                flag_set_ppor1=false;
            end
        end
        check=flag_worked;
        return;
    end %of recursive function



end






%%%%%%  reorder better by exchanging two pieces ######

function nnew_options=byscore_reorder(original,options,show_state,show_index,timlimit)

    nnew_options=options;
    flag_initial=true;
    meshr=16;
    meshc=16;
    strategy=1;  %1 - stop nucleation with the first dead end
    %2 - skip pieces that could not fit the nucleation
    %3 - suggest
    strategyb=0; % 0- regular mode, 1- suggestion mode
    orientprice=6;
    orientlock=0;
    lookconserv=1;
    loopscoremax=0;
    loopcontinue=1;
    user_priceconserv=0.05;%0.5;
    pricecontinuity=0.1;%0.5;
    ncoremax=0;
    scoremax=0;
    ncore=0;
    score=0;
    deltascore=0;
    registerN=0;
    flname=0;
    %symbolname  (num)=' ';
    symbolname{1}=-1;
    ncorev=[0 0 0 0 0 0 0 0];
    phase=0;
    r0v(phase+1)=9;
    c0v(phase+1)=8;
    phase=1;
    r0v(phase+1)=3;
    c0v(phase+1)=14;
    phase=2;
    r0v(phase+1)=14;
    c0v(phase+1)=3;
    phase=3;
    r0v(phase+1)=3;
    c0v(phase+1)=3;
    phase=4;
    r0v(phase+1)=14;
    c0v(phase+1)=14;
    % original {r,c}=[1 2 3 4];
    clue=[1 1];%changed from 1 1 %default, overloaded by origin7.mat or above
    %load('C:\Users\seifer\Documents\Matlab\eternity\origin6.mat','-mat');
    piece=zeros(16*16,4);
    for row=1:16
        for col=1:16
            ob=original{row,col};
            ind=(col-1)+(row-1)*16+1;
            piece(ind,1:4)=ob;
        end
    end
     

    neworiginal=fill_neworiginal(show_state,show_index,options,piece); %fill like before display, according to options
    for row=1:16
        for col=1:16
            if sum(neworiginal{row,col}==[0 0 0 0])==4
                neworiginal{row,col}=original{row,col};
            end
        end
    end

    



    %original will be modified to have the borders on the right sides, for
    %correct names and conditions of pieces on the borders
    rsvn=0;
    for r=1:16
        c=1;
        if neworiginal{r,c}(4)~=0
            rsvn=rsvn+1;
            rsvlist{rsvn}=neworiginal{r,c};
            neworiginal{r,c}=[0 0 0 0];
        end
        c=16;
        if neworiginal{r,c}(2)~=0 
            rsvn=rsvn+1;
            rsvlist{rsvn}=neworiginal{r,c};
            neworiginal{r,c}=[0 0 0 0];
        end
    end
    for c=1:16
        r=1;
        if neworiginal{r,c}(1)~=0
            rsvn=rsvn+1;
            rsvlist{rsvn}=neworiginal{r,c};
            neworiginal{r,c}=[0 0 0 0];
        end
        r=16;
        if neworiginal{r,c}(3)~=0
            rsvn=rsvn+1;
            rsvlist{rsvn}=neworiginal{r,c};
            neworiginal{r,c}=[0 0 0 0];
        end
    end
    for r=1:16
        c=1;
        if sum(neworiginal{r,c}==[0 0 0 0])==4
            for rs=1:rsvn
                if rsvlist{rs}(4)==0 && (rsvlist{rs}(1)==0 || r~=1) && (rsvlist{rs}(3)==0 || r~=16) && (rsvlist{rs}(1)~=0 || r==1) && (rsvlist{rs}(3)~=0 || r==16)
                    neworiginal{r,c}=rsvlist{rs};
                    rsvlist{rs}=[-1 -1 -1 -1];
                    break
                end
            end
        end
        c=16;
        if sum(neworiginal{r,c}==[0 0 0 0])==4
            for rs=1:rsvn
                if rsvlist{rs}(2)==0 && (rsvlist{rs}(1)==0 || r~=1) && (rsvlist{rs}(3)==0 || r~=16) && (rsvlist{rs}(1)~=0 || r==1) && (rsvlist{rs}(3)~=0 || r==16)
                    neworiginal{r,c}=rsvlist{rs};
                    rsvlist{rs}=[-1 -1 -1 -1];
                    break
                end
            end
        end
    end
    for c=1:16
        r=1;
        if sum(neworiginal{r,c}==[0 0 0 0])==4
            for rs=1:rsvn
                if rsvlist{rs}(1)==0 && (rsvlist{rs}(4)==0 || c~=1) && (rsvlist{rs}(2)==0 || c~=16) && (rsvlist{rs}(4)~=0 || c==1) && (rsvlist{rs}(2)~=0 || c==16)
                    neworiginal{r,c}=rsvlist{rs};
                    rsvlist{rs}=[-1 -1 -1 -1];
                    break
                end
            end
        end
        r=16;
        if sum(neworiginal{r,c}==[0 0 0 0])==4
            for rs=1:rsvn
                if rsvlist{rs}(3)==0 && (rsvlist{rs}(4)==0 || c~=1) && (rsvlist{rs}(2)==0 || c~=16) && (rsvlist{rs}(4)==0 || c~=1) && (rsvlist{rs}(2)==0 || c~=16) && (rsvlist{rs}(4)~=0 || c==1) && (rsvlist{rs}(2)~=0 || c==16)
                    neworiginal{r,c}=rsvlist{rs};
                    rsvlist{rs}=[-1 -1 -1 -1];
                    break
                end
            end
        end
    end
    
    tic;
        
    while score<480 && toc<timlimit  
    
    % symbol gray must be zero !!!
    
    lock=zeros(meshr,meshc);
    lock(9,8)=1; %constraint
    if clue(2)==1
       lock(3,14)=1;%clue2
       lock(3,3)=1;%clue4
    end
    if clue(1)==1
        lock(14,3)=1;%clue1
        lock(14,14)=1;%clue3
    end
    % the center of the spiral must be locked
    ncore=1;
    %to read original{1,1}(1) , etc...
    
    current=neworiginal;
    for r=1:meshr
        for c=1:meshc
            current_pos(r,c,1)=r;
            current_pos(r,c,2)=c;
            current_pos(r,c,3)=0; %angle
        end
    end
    
    priceconserv=user_priceconserv;  %defined by the user as he rates his initial setup
    
     score=0;
     for score_r=1:meshr-1
         for score_c=1:meshc-1
            score=score+( current{score_r,score_c}(3)==current{score_r+1,score_c}(1) ) + ( current{score_r,score_c}(2)==current{score_r,score_c+1}(4) ) ;
         end
     end
     score_c=meshc;
     for score_r=1:meshr-1
             score=score+ ( current{score_r,score_c}(3)==current{score_r+1,score_c}(1) ) ;
    end
     score_r=meshr;
     for score_c=1:meshc-1
             score=score+ ( current{score_r,score_c}(2)==current{score_r,score_c+1}(4) )  ;
    end
    scoreprev=score;
    scoremax=0;
    
    if flag_initial
        flag_initial=false;
        disp(sprintf('Initial score=%d   ',score));
    end

    countloop=0;
    
    phase=0;
    
    while (phase>0 || (score>=scoremax)  || (countloop==0) || (loopcontinue==1) ) && toc<timlimit
    current_back=current;
    current_pos_back=current_pos;
    
    conserv=zeros(meshr,meshc);
    for cons_r=1:meshr
         for cons_c=1:meshc
             test=( current{cons_r,cons_c}(3)==(cons_r<meshr)*current{min(cons_r+1,meshr),cons_c}(1) ) + ( current{cons_r,cons_c}(2)==(cons_c<meshc)*current{cons_r,min(cons_c+1,meshc)}(4) ) ;
             test=test+( current{cons_r,cons_c}(1)==(cons_r>1)*current{max(cons_r-1,1),cons_c}(3) ) + ( current{cons_r,cons_c}(4)==(cons_c>1)*current{cons_r,max(cons_c-1,1)}(2) ) ;
              if test==4 && sqrt((cons_r-r0v(phase+1))^2+(cons_c-c0v(phase+1))^2)>=(sqrt(ncorev(phase+1))/2)
                  conserv(cons_r,cons_c)=1;
              end
         end
    end
    
    if clue(1)==1 & clue(2)==1
        if phase==1
            r0=r0v(phase+1);%3, 9;
            c0=c0v(phase+1);%14, 8;
            strategy=1;
        elseif phase==2
            r0=r0v(phase+1);
            c0=c0v(phase+1);
            strategy=1;
        elseif phase==3
            r0=r0v(phase+1);
            c0=c0v(phase+1);
            strategy=1;
        elseif phase==4
            r0=r0v(phase+1);
            c0=c0v(phase+1);
            strategy=1;
        elseif phase==0
            r0=r0v(phase+1);
            c0=c0v(phase+1);
            strategy=1+(countloop>=0);
            %orientlock=(countloop<=6);
        %elseif phase==3
         %   r0=9;
         %   c0=8;
         %   strategy=2;
          %  priceconserv=0.5;
        end
        phase=mod(phase+1,5);
        
    else
         phase=0;
        if phase==0
            r0=r0v(phase+1);
            c0=c0v(phase+1);
            strategy=1+(countloop<2 | countloop>=15);
            %orientlock=(countloop<=10);
        end
        phase=0;%mod(phase+1,2);
    
    end
    
    if strategyb==1
          pricecontinuity=0;
          priceconserv=0;
    end
    
    ncore=1;
    
    countloop=countloop+1 ;   
    temporary=current;
    %stack=''; %temporary replacment track
    stackN=0;
    n=0;
    for r=2:meshr-1
        for c=2:meshc-1
            if lock(r,c)==0  % here is locked pieces that will not be replaced
                n=n+1;
                picklist(n,1:2)=[r c];
            end
        end
    end
    picklistN=n;
    n=0;
    for r=2:meshr-1
        n=n+1;
        picklist_border(n,1:2)=[r 1];
        n=n+1;
        picklist_border(n,1:2)=[r meshc];
    end
    for c=2:meshc-1
        n=n+1;
        picklist_border(n,1:2)=[1 c];
        n=n+1;
        picklist_border(n,1:2)=[meshr c];
    end
    picklist_borderN=n;
     picklist_corner(1,1:2)=[1 1];
     picklist_corner(2,1:2)=[1 meshc];
     picklist_corner(3,1:2)=[meshr 1];
     picklist_corner(4,1:2)=[meshr meshc];
     picklist_cornerN=4;
     
    % main loop: choose no 1 tag along a spiral
    first_r=r0;
    first_c=c0;
    spiralborder_r1=first_r; % the top of filled row in spiral
    spiralborder_r2=first_r;%+1; % the filled bottom line
    spiralborder_c2=first_c; % the rightest filled column
    spiralborder_c1=first_c; % the leftest filled column
    
    choices=0; %reset number of best choices
    
    % anticlockwise
    while first_r<=16 && (first_r>1 || first_c>1)
            stepdir_c=((1)*(first_r>spiralborder_r2)+(-1)*(first_r<spiralborder_r1));
            stepdir_r=((1)*(first_c<spiralborder_c1) +(-1)*(first_c>spiralborder_c2));
            c=first_c+stepdir_c*(first_c>=spiralborder_c1-max(0,stepdir_c) & first_c<=spiralborder_c2-min(0,stepdir_c) )+(first_r==r0 & first_c==c0);
            r=first_r+stepdir_r*(first_r>=spiralborder_r1-max(0,stepdir_r) & first_r<=spiralborder_r2-min(0,stepdir_r)) ;
            if r==first_r+1 && r==spiralborder_r1
                spiralborder_r1=max(1,spiralborder_r1-1);
            end
            if r==first_r-1 && r==spiralborder_r2
                spiralborder_r2=min(meshr,spiralborder_r2+1);
            end
            if c==first_c+1 && c==spiralborder_c1
                spiralborder_c1=max(1,spiralborder_c1-1);
            end
            if c==first_c-1 && c==spiralborder_c2
                spiralborder_c2=min(meshc,spiralborder_c2+1);
            end
            
          if c<1 || c>meshc
            c=1*(c<1)+meshc*(c>meshc);
            if r<r0
                r=spiralborder_r2+1;
            else
                r=spiralborder_r1-1;
            end
            if c==1 && r>meshr
                c=spiralborder_c2+1;
                r=meshr;
                spiralborder_r1=max(1,spiralborder_r1-1);
            end
            if c==meshc && r<1 
                c=spiralborder_c1-1;
                r=1;
                spiralborder_r2=min(meshr,spiralborder_r2+1);
            end
        end
        if r<1
            r=spiralborder_r2+1;
            if r>meshr
                break
            end
        end
        if r>meshr
            r=spiralborder_r1-1;
            if r<1
                break
            end
        end
        if c>meshc || c<1
            break
        end
            
        first_r=r;
        first_c=c;
        value1=current{first_r,first_c};
        
           if lock(first_r,first_c)==0
       
            ncore=ncore+1;   
               
         if (first_r==1 && first_c==1)  ||  (first_r==meshr && first_c==meshc)  ||  (first_r==1 && first_c==meshc)  ||  (first_r==meshr && first_c==1) 
            region=3; % corner (to replace only with corner)
            Npicks=picklist_cornerN;
         elseif first_r==1 || first_r==meshr || first_c==1 || first_c==meshc
            region=2; % border (to replace only with border not corner)
            Npicks=picklist_borderN;
         else
             region=1; % interior (regular, to replace with interior pieces)
            Npicks=picklistN;
         end
         
         look=[0 0 0 0]; %look up, right, down, left to match with existing spiral (first's neighbors)
         if first_r>spiralborder_r2 && first_c<=spiralborder_c2
             look(1)=1;
             if first_c>spiralborder_c1
                  look(4)=1;
             end
         elseif first_r<spiralborder_r1 && first_c>=spiralborder_c1
             look(3)=1;
             if first_c<spiralborder_c2
                 look(2)=1;
             end
         end
         if first_c>spiralborder_c2 && first_r>=spiralborder_r1
             look(4)=1;
             if first_r<spiralborder_r2
                 look(3)=1;
             end
         elseif first_c<spiralborder_c1 && first_r<=spiralborder_r2
             look(2)=1;
               if first_r>spiralborder_r1
                 look(1)=1;
             end
         end
         
         if strategyb==1
            if  ((spiralborder_r1==1 && spiralborder_c2==meshc) ||  (spiralborder_r2==meshr && spiralborder_c1==1))
                if (c==1) || (c==meshc)
                    look(1)=0;
                    look(3)=0;
                elseif (r==1) || (r==meshr)
                    look(2)=0;
                    look(4)=0;
                end
             end
         end
         
         
         if first_r==1
             look(1)=0;
             neighborup=0;   % zero is symbol gray
         else
             neighborup=temporary{first_r-1,first_c}(3);
         end
         if first_r==meshr
             look(3)=0;
              neighbordown=0;
         else
             neighbordown=temporary{first_r+1,first_c}(1);
         end
         if first_c==1
             look(4)=0;
              neighborleft=0;
         else
             neighborleft=temporary{first_r,first_c-1}(2);
         end
         if first_c==meshc
             look(2)=0;
             neighborright=0;
         else
             neighborright=temporary{first_r,first_c+1}(4);
         end
          
       %      sprintf('r=%g,  c=%g,     look=[%g  %g  %g  %g  ]',first_r,first_c,look(1),look(2),look(3),look(4))
         
     profitN=0;
         
     for nn=1:Npicks  %look for replacement to first piece including the same (rotated or not)
         if region==1
             pick=picklist(nn,1:2);
         elseif region==2
             pick=picklist_border(nn,1:2);
         elseif region==3
             pick=picklist_corner(nn,1:2);
         end
        second_r=pick(1);
        second_c=pick(2);
        selfreplace=((second_r==first_r) & (second_c==first_c));  %1- self replace - skip other replacements
             
        
        fslook(1)=first_r>1;
        fslook(3)=first_r<meshr;
        fslook(4)=first_c>1;
        fslook(2)=first_c<meshc;
        slook(1)=second_r>1;
        slook(3)=second_r<meshr;
        slook(4)=second_c>1;
        slook(2)=second_c<meshc;
       
             indrot=0:(3*selfreplace+15*(~selfreplace));
             indrot_mark=zeros(size(indrot));
             if sum(fslook)<4 || sum(slook)<4
                 for nrot=indrot
                     first_rotate=mod(nrot,4);  % 0-no, 1-90, 2-180, 3-270
                     second_rotate=fix(nrot/4);   % 0-no, 1-90, 2-180, 3-270
                     fslooknow=fslook((first_rotate==0)*(1:4)+(first_rotate==1)*[4 1:3]  + (first_rotate==2)*[3:4 1:2]  + (first_rotate==3)*[2:4 1] ) ;    
                     slooknow=slook((second_rotate==0)*(1:4)+(second_rotate==1)*[4 1:3]  + (second_rotate==2)*[3:4 1:2]  + (second_rotate==3)*[2:4 1] ) ;    
                     if sum(fslooknow==slook)==4 &&  sum(slooknow==fslook)==4
                     indrot_mark(nrot+1)=1;
                     else
                     indrot_mark(nrot+1)=0;
                     end
                 end
                 if sum(indrot_mark(1:length(indrot)))>0
                     indrot=indrot(indrot_mark(1:length(indrot))==1);
                 end
             end
        
         for nrot=indrot
         first_rotate=mod(nrot,4);  % 0-no, 1-90, 2-180, 3-270
         second_rotate=fix(nrot/4);   % 0-no, 1-90, 2-180, 3-270
         value2=current{second_r,second_c};
         valuenow1=value1((first_rotate==0)*(1:4)+(first_rotate==1)*[4 1:3]  + (first_rotate==2)*[3:4 1:2]  + (first_rotate==3)*[2:4 1] ) ;    
         valuenow2=value2((second_rotate==0)*(1:4)+(second_rotate==1)*[4 1:3]  + (second_rotate==2)*[3:4 1:2]  + (second_rotate==3)*[2:4 1] ) ;    
         
        pass1=1;
         if look(1)
             pass1=pass1 & valuenow2(1)==neighborup;
         end
         if look(2)
             pass1=pass1 & valuenow2(2)==neighborright;
         end
          if look(3)
             pass1=pass1 & valuenow2(3)==neighbordown;
         end
         if look(4)
             pass1=pass1 & valuenow2(4)==neighborleft;
         end
      
         forth_r=[0 0 0 0];
         forth_c=[0 0 0 0];
         third_rot=[0 0 0 0];
        choice34=0; %number of choices met in the  3-4 permutation
      
        
         % now if pass1==1 it worth to check replacement
         if pass1
            if selfreplace==0
    
             temporary{first_r,first_c}=valuenow2;
             temporary{second_r,second_c}=valuenow1;
             stackN=1;
             stack{stackN}=[first_r  first_c];
             stackN=2;
             stack{stackN}=[second_r  second_c];
             
                
                scheck=[0 0 0 0];
                 if second_r==1
                     %slook(1)=0;
                     neighborup2=0;
                     third_r(1)=0;
                     third_c(1)=0;
                 else
                     neighborup2=temporary{second_r-1,second_c}(3);
                     third_r(1)=second_r-1;
                     third_c(1)=second_c;
                     scheck(1)=( neighborup2==temporary{second_r,second_c}(1) );
                     %slook(1)=1;
                 end
                 if second_r==meshr
                     %slook(3)=0;
                      neighbordown2=0;
                     third_r(3)=0;
                     third_c(3)=0;
                else
                     neighbordown2=temporary{second_r+1,second_c}(1);
                     third_r(3)=second_r+1;
                     third_c(3)=second_c;
                     scheck(3)=( neighbordown2==temporary{second_r,second_c}(3) );
                     %slook(3)=1;
                 end
                 if second_c==1
                     %slook(4)=0;
                     neighborleft2=0;
                     third_r(4)=0;
                     third_c(4)=0;
                 else
                     neighborleft2=temporary{second_r,second_c-1}(2);
                     third_r(4)=second_r;
                     third_c(4)=second_c-1;
                     scheck(4)=( neighborleft2==temporary{second_r,second_c}(4) );
                     %slook(4)=1;
                 end
                 if second_c==meshc
                     %slook(2)=0;
                     neighborright2=0;
                     third_r(2)=0;
                     third_c(2)=0;
                 else
                     neighborright2=temporary{second_r,second_c+1}(4);
                     third_r(2)=second_r;
                     third_c(2)=second_c+1;
                     scheck(2)=( neighborright2==temporary{second_r,second_c}(2) );
                     %slook(2)=1;
                 end
    
                 for dirc=1:4
                     if ((third_r(dirc)==first_r) && (third_c(dirc)==first_c))  || (strategyb==1) 
                         third_r(dirc)=0;
                         third_c(dirc)=0;
                         slook(dirc)=0;
                     end
                 end
                 
                
                 for dirc=1:4
                     profit34N=0;
                     if slook(dirc)==1
                            tlook=[1 1 1 1]; %% third's pointer neighbors
                            if third_r(dirc)==1
                                tlook(1)=0;
                                neighborup3=0;
                            else
                                neighborup3=temporary{third_r(dirc)-1,third_c(dirc)}(3);
                            end
                             if third_r(dirc)==meshr
                                tlook(3)=0;
                                 neighbordown3=0;
                             else
                                neighbordown3=temporary{third_r(dirc)+1,third_c(dirc)}(1);
                          end
                             if third_c(dirc)==meshc
                                tlook(2)=0;
                                  neighborright3=0;
                            else
                                neighborright3=temporary{third_r(dirc),third_c(dirc)+1}(4);
                           end
                             if third_c(dirc)==1
                                tlook(4)=0;
                                 neighborleft3=0;
                              else
                                neighborleft3=temporary{third_r(dirc),third_c(dirc)-1}(2);
                           end
    
                         if (third_r(dirc)==1 && third_c(dirc)==1)  ||  (third_r(dirc)==meshr && third_c(dirc)==meshc)  ||  (third_r(dirc)==1 && third_c(dirc)==meshc)  ||  (third_r(dirc)==meshr && third_c(dirc)==1)
                            region4=3; % corner (to replace only with corner)
                            Npicks4=picklist_cornerN;
                         elseif third_r(dirc)==1 || third_r(dirc)==meshr || third_c(dirc)==1 || third_c(dirc)==meshc
                            region4=2; % border (to replace only with border not corner)
                            Npicks4=picklist_borderN;
                         else
                             region4=1; % interior (regular, to replace with interior pieces)
                            Npicks4=picklistN;
                         end
                         pickable3(dirc)=0;
                         for nnn=1:Npicks4
                              if region4==1
                                 pick4t=picklist(nnn,1:2);
                              elseif region4==2
                                 pick4t=picklist_border(nnn,1:2);
                              elseif region4==3
                                 pick4t=picklist_corner(nnn,1:2);
                              end
                               if sum(pick4t==[third_r(dirc)  third_c(dirc)])==2
                                   pickable3(dirc)=1;
                                   break
                               end
                         end % for nnn
                            
                     if ~scheck(dirc) &&  (pickable3(dirc)==1) && (slook(dirc)==1)
                          for nnn=1:Npicks4  %look for forth piece including the same (rotated or not)
                              if region4==1
                                 pick4=picklist(nnn,1:2);
                              elseif region4==2
                                 pick4=picklist_border(nnn,1:2);
                              elseif region4==3
                                 pick4=picklist_corner(nnn,1:2);
                              end
                              usedotherdir=0;
                              for tempdirc=1:dirc-1
                                    usedotherdir=usedotherdir  ||  (pick4(1)==forth_r(tempdirc) && pick4(2)==forth_c(tempdirc)); 
                              end
                              for tempdirc=1:4
                                    usedotherdir=usedotherdir  ||  (pick4(1)==third_r(tempdirc) && pick4(2)==third_c(tempdirc)); 
                              end
                              if (pick4(1)==first_r && pick4(2)==first_c) ||  (pick4(1)==second_r && pick4(2)==second_c)  || usedotherdir
                                 % skip
                              else
                                %only here check profit34 for each dir indepen
                                %needs to choose best rotation here
                                %rot3 - how to rotate the piece now at pos4
                                %the roatation always refers to the value
                                selfreplace2=(pick4(1)==third_r(dirc) & pick4(2)==third_c(dirc));
                                flook=[1 1 1 1]; %% forth's pointer neighbors
                                if pick4(1)==1
                                    flook(1)=0;
                                    neighborup4=0;
                                else
                                    neighborup4=temporary{pick4(1)-1,pick4(2)}(3);
                                end
                                 if pick4(1)==meshr
                                    flook(3)=0;
                                    neighbordown4=0;
                                else
                                    neighbordown4=temporary{pick4(1)+1,pick4(2)}(1);
                                end
                                 if pick4(2)==meshc
                                    flook(2)=0;
                                     neighborright4=0;
                                else
                                    neighborright4=temporary{pick4(1),pick4(2)+1}(4);
                               end
                                 if pick4(2)==1
                                    flook(4)=0;
                                    neighborleft4=0;
                                else
                                    neighborleft4=temporary{pick4(1),pick4(2)-1}(2);
                                 end
                                 
                                 
                                         ind34rot=0:(3*selfreplace2+15*(~selfreplace2));
                                         ind34rot_mark=zeros(size(ind34rot));
                                         if sum(tlook)<4 || sum(flook)<4
                                             for nrot34=ind34rot
                                                 forth_rotate=mod(nrot34,4);  % 0-no, 1-90, 2-180, 3-270
                                                 third_rotate=fix(nrot34/4);   % 0-no, 1-90, 2-180, 3-270
                                                 tlooknow=tlook((third_rotate==0)*(1:4)+(third_rotate==1)*[4 1:3]  + (third_rotate==2)*[3:4 1:2]  + (third_rotate==3)*[2:4 1] ) ;    
                                                 flooknow=flook((forth_rotate==0)*(1:4)+(forth_rotate==1)*[4 1:3]  + (forth_rotate==2)*[3:4 1:2]  + (forth_rotate==3)*[2:4 1] ) ;    
                                                 if sum(tlooknow==flook)==4 &&  sum(flooknow==tlook)==4
                                                 ind34rot_mark(nrot34+1)=1;
                                                 else
                                                 ind34rot_mark(nrot34+1)=0;
                                                 end
                                             end
                                             if sum(ind34rot_mark(1:length(ind34rot)))>0
                                                 ind34rot=ind34rot(ind34rot_mark(1:length(ind34rot))==1);
                                             end
                                         end
                                         for nrot34=ind34rot
                                         forth_rotate=mod(nrot34,4);  % 0-no, 1-90, 2-180, 3-270
                                         third_rotate=fix(nrot34/4);   % 0-no, 1-90, 2-180, 3-270
                                         %values alway refer to original
                                         %position; neighbors refer to original
                                         %position; contents are switched 3--4
                                         value3=temporary{third_r(dirc),third_c(dirc)};    
                                         value4=temporary{pick4(1),pick4(2)};
                                         valuenow3=value3((third_rotate==0)*(1:4)+(third_rotate==1)*[4 1:3]  + (third_rotate==2)*[3:4 1:2]  + (third_rotate==3)*[2:4 1] ) ;    
                                         valuenow4=value4((forth_rotate==0)*(1:4)+(forth_rotate==1)*[4 1:3]  + (forth_rotate==2)*[3:4 1:2]  + (forth_rotate==3)*[2:4 1] ) ;    
                                          if pick4(1)==third_r(dirc)
                                                if pick4(2)==third_c(dirc)+1
                                                     neighborright3=valuenow3(4);
                                                     neighborleft4=valuenow4(2);
                                               elseif pick4(2)==third_c(dirc)-1
                                                     neighborleft3=valuenow3(2);
                                                     neighborright4=valuenow4(4);
                                               end
                                          elseif pick4(2)==third_c(dirc)
                                                if pick4(1)==third_r(dirc)-1
                                                     neighborup3=valuenow3(3);
                                                     neighbordown4=valuenow4(1);
                                               elseif pick4(1)==third_r(dirc)+1
                                                     neighbordown3=valuenow3(1);
                                                     neighborup4=valuenow4(3);
                                                end
                                          end
                                         
                                         if selfreplace2==1 
                                              profit34rot(nrot34+1)=((valuenow4(1)==neighborup3) - (value3(1)==neighborup3))+((valuenow4(2)==neighborright3) -( value3(2)==neighborright3)) ...
                                                + ((valuenow4(3)==neighbordown3) -( value3(3)==neighbordown3))+((valuenow4(4)==neighborleft3) - (value3(4)==neighborleft3)) ...
                                                +orientlock*(-orientprice*(forth_rotate~=0)  -orientprice*( third_rotate~=0));
                                         else
                                               profit34rot(nrot34+1)=((valuenow4(1)==neighborup3) - (value3(1)==neighborup3))+((valuenow4(2)==neighborright3) - (value3(2)==neighborright3)) ...
                                                + ((valuenow4(3)==neighbordown3) - (value3(3)==neighbordown3))+((valuenow4(4)==neighborleft3) - (value3(4)==neighborleft3));
                                               profit34rot(nrot34+1)=profit34rot(nrot34+1)+((valuenow3(1)==neighborup4) - (value4(1)==neighborup4))+((valuenow3(2)==neighborright4) - (value4(2)==neighborright4)) ...
                                                + ((valuenow3(3)==neighbordown4) - (value4(3)==neighbordown4))+((valuenow3(4)==neighborleft4) - (value4(4)==neighborleft4))  ...
                                                +orientlock*(-orientprice*(forth_rotate~=0) -orientprice*( third_rotate~=0));
                                        end % if selfreplace2
                                     end %for nrot34
                                 choice34rot_compare=ind34rot(profit34rot(1:length(ind34rot))==max(profit34rot(1:length(ind34rot))));
                                 rotchoices=length(choice34rot_compare);
                                 choice34rot=choice34rot_compare(fix(rotchoices*rand(1))+1); % recover nrot34
                                rot4=mod(choice34rot,4);  % 0-no, 1-90, 2-180, 3-270
                                rot3=fix(choice34rot/4);   % 0-no, 1-90, 2-180, 3-270
                                if selfreplace2==1
                                    rot3=5;
                                end
                                
                                profit34N=profit34N+1;
                                profit34_do{profit34N}=[pick4 rot3 rot4 rotchoices];
                                profit34(profit34N)=profit34rot(choice34rot+1);%?
                                
                             end % if (pick
                          end %for nnn
                     end % if  ~scheck(dirc) &  pickable3(dirc)==1
                     if profit34N==0
                            forth_r(dirc)=0;
                            forth_c(dirc)=0;
                            third_rot(dirc)=0;
                            forth_rot(dirc)=0;
                     else
                            indcompare=1:profit34N;
                             prof34compare=indcompare(profit34(1:profit34N)==max(profit34(1:profit34N)));
                             randchoice=prof34compare(fix(length(prof34compare)*rand(1))+1);
                             choice34=choice34+length(prof34compare); % the overalchoices count
                             choice34=choice34+profit34_do{randchoice}(5); % the overalchoices count with rotchoices
                            forth_r(dirc)=profit34_do{randchoice}(1);  %pick4(1)
                             forth_c(dirc)=profit34_do{randchoice}(2); %pick4(2)
                             third_rot(dirc)=profit34_do{randchoice}(3); %rot3
                             forth_rot(dirc)=profit34_do{randchoice}(4); %rot4
                             %enter to <temporary> table 
                                  value3=temporary{third_r(dirc),third_c(dirc)};    
                                  value4=temporary{forth_r(dirc),forth_c(dirc)};
                                  valuenow3=value3((third_rot(dirc)==0 | third_rot(dirc)==5)*(1:4)+(third_rot(dirc)==1)*[4 1:3]  + (third_rot(dirc)==2)*[3:4 1:2]  + (third_rot(dirc)==3)*[2:4 1] ) ;    
                                  valuenow4=value4((forth_rot(dirc)==0)*(1:4)+(forth_rot(dirc)==1)*[4 1:3]  + (forth_rot(dirc)==2)*[3:4 1:2]  + (forth_rot(dirc)==3)*[2:4 1] ) ;    
                                 temporary{third_r(dirc),third_c(dirc)}=valuenow4;  %here is the switch
                                 temporary{forth_r(dirc),forth_c(dirc)}=valuenow3;
                                 stackN=stackN+1;
                                 stack{stackN}=[third_r(dirc)  third_c(dirc)];
                                 stackN=stackN+1;
                                 stack{stackN}=[forth_r(dirc)  forth_c(dirc)];
    
                     end %if profit34N
                     
                     else %of if slook(dirc)==1 (here ==0)
                            forth_r(dirc)=0;
                            forth_c(dirc)=0;
                            third_rot(dirc)=0;
                            forth_rot(dirc)=0;
    
                     end %if slook(dirc)==1
                  end %for dirc
      
                  
            elseif   selfreplace==1
            
                %ido nothing here
                
            end %if ~selfreplace
             
            % calculate profite of permuatations and get temporary table reset
             % keep first_(r,c,rotate), second, third, forth  in the profit table
            profitN=profitN+1;
            if selfreplace==1
                  profit_do{profitN}=[selfreplace first_r first_c first_rotate];
                  profit(profitN)=((valuenow1(1)==neighborup) - (value1(1)==neighborup))+((valuenow1(2)==neighborright) - (value1(2)==neighborright)) ...
                    + ((valuenow1(3)==neighbordown) - (value1(3)==neighbordown))+((valuenow1(4)==neighborleft) - (value1(4)==neighborleft))...
                    +priceconserv*conserv(first_r,first_c) +orientlock*(-orientprice*(first_rotate~=0) -orientprice*(second_rotate~=0));
            else
                 profit_do{profitN}=[selfreplace  first_r  first_c  first_rotate  second_r  second_c  second_rotate  ...
                 forth_r(1)  forth_c(1)  forth_rot(1)  third_rot(1)  forth_r(2)  forth_c(2)  forth_rot(2) third_rot(2) ...
                 forth_r(3)  forth_c(3)  forth_rot(3)  third_rot(3)  forth_r(4)  forth_c(4)  forth_rot(4) third_rot(4) ...
                 choice34 ];
               
            %lostconserv=0;
            %if  lookconserv==1
            %   lostconserv=conserv(second_r,second_c);
            %     if lostconserv==0
            %        for dirctemp=1:4
            %            if forth_r(dirctemp)>0 & forth_c(dirctemp)>0
            %                    lostconserv=lostconserv+(third_rot(dirctemp)<5)*conserv(forth_r(dirctemp),forth_c(dirctemp));
            %            end
            %         end
            %     end
            %end
    
            newconserv=zeros(meshr,meshc);
            scorenow=0;
             scoreprev=0;
             for score_r=1:meshr-1
                 for score_c=1:meshc-1
                    scoreprev=scorenow+( current{score_r,score_c}(3)==current{score_r+1,score_c}(1) ) + ( current{score_r,score_c}(2)==current{score_r,score_c+1}(4) ) ;
                    scorenow=scorenow+( temporary{score_r,score_c}(3)==temporary{score_r+1,score_c}(1) ) + ( temporary{score_r,score_c}(2)==temporary{score_r,score_c+1}(4) ) ;
                    test=( temporary{score_r,score_c}(3)==(score_r<meshr)*temporary{min(score_r+1,meshr),score_c}(1) ) + ( temporary{score_r,score_c}(2)==(score_c<meshc)*temporary{score_r,min(score_c+1,meshc)}(4) ) ;
                    test=test+( temporary{score_r,score_c}(1)==(score_r>1)*temporary{max(score_r-1,1),score_c}(3) ) + ( temporary{score_r,score_c}(4)==(score_c>1)*temporary{score_r,max(score_c-1,1)}(2) ) ;
                    if test==4 && sqrt((score_r-r0)^2+(score_c-c0)^2)>=(sqrt(ncorev(phase+1))/2)
                        newconserv(score_r,score_c)=1;
                    end
                 end
             end
             profit(profitN)=scorenow-scoreprev+pricecontinuity*(sum(sum(newconserv))-sum(sum(conserv))) +orientlock*( -orientprice*(first_rotate~=0) -orientprice*(second_rotate~=0));
             
             
            end %if selfreplace
            profit(profitN)=profit(profitN)+(rand(1)>0.7)*1;   %enetering arbitrariness in selection of moves
            
            
            for n=1:stackN
                temporary{stack{n}(1),stack{n}(2)}=current{stack{n}(1),stack{n}(2)};
            end
            stackN=0;
            
         end % if pass1
         
         end %for nrot
         end % for nn   
     
     % Choose the best permuation (or randomly between equal value) and show
     % number of choices that were.  ---> choiceN,  choices
     skiperfect=0;
     
     if profitN==0
         if strategy==1
            break
         elseif strategy==2
             skiperfect=1;
         end
     end
     
     if skiperfect==0
     
     indcompare=1:profitN;
     profitcompare=indcompare(profit(indcompare)==max(profit(indcompare)));
     choiceN=profitcompare(fix(length(profitcompare)*rand(1))+1);
     if profit_do{choiceN}(1)==0 %not selfreplace
        choices=choices+length(profitcompare)*profit_do{choiceN}(24) ;  %choice12*choice34
     end
    
     % update current and temporary with the accepted changes + remember
      % permutations in current_pos
     if profit_do{choiceN}(1)==1  % selfreplace
         first_r=profit_do{choiceN}(2);
         first_c=profit_do{choiceN}(3);
         a=current{first_r , first_c};
         rot=profit_do{choiceN}(4)*(orientlock==0);
         current{first_r , first_c}=a((rot==0)*(1:4)+(rot==1)*[4 1:3]  + (rot==2)*[3:4 1:2]  + (rot==3)*[2:4 1] ) ;
         current_pos(first_r , first_c,3)=mod(current_pos(first_r , first_c,3)+rot,4);
     else
         first_r=profit_do{choiceN}(2);
         first_c=profit_do{choiceN}(3);
         first_a=current{first_r,first_c};
         first_rot=profit_do{choiceN}(4)*(orientlock==0);
         current_pos(first_r , first_c,3)=mod(current_pos(first_r , first_c,3)+first_rot,4);
         first_track=current_pos(first_r , first_c,1:3);
         second_r=profit_do{choiceN}(5);
         second_c=profit_do{choiceN}(6);
         second_a=current{second_r,second_c};
         second_rot=profit_do{choiceN}(7)*(orientlock==0);
         current_pos(second_r , second_c,3)=mod(current_pos(second_r , second_c,3)+second_rot,4);
         second_track=current_pos(second_r , second_c,1:3);
         slook(1)=second_r>1;
         slook(3)=second_r<meshr;
         slook(4)=second_c>1;
         slook(2)=second_c<meshc;
         for dirc=1:4
            if slook(dirc)==1
             n1=8+(dirc-1)*4;
             n2=9+(dirc-1)*4;
             n3=10+(dirc-1)*4;
             n4=11+(dirc-1)*4;
             if profit_do{choiceN}(n1)~=0
                 forth_r(dirc)=profit_do{choiceN}(n1);
                 forth_c(dirc)=profit_do{choiceN}(n2);
                 if (forth_r(dirc)>0) && (forth_c(dirc)>0)
                     forth_a{dirc}=current{forth_r(dirc),forth_c(dirc)};
                     forth_rot(dirc)=profit_do{choiceN}(n3)*(orientlock==0);
                     current_pos(forth_r(dirc) , forth_c(dirc),3)=mod(current_pos(forth_r(dirc) , forth_c(dirc),3)+forth_rot(dirc),4);
                     forth_track{dirc}=current_pos(forth_r(dirc) , forth_c(dirc),1:3);
                     third_rot(dirc)=profit_do{choiceN}(n4)*(orientlock==0)+5*(profit_do{choiceN}(n4)==5)*(orientlock==1);
                     third_r(dirc)=second_r-(dirc==1)+(dirc==3);
                     third_c(dirc)=second_c-(dirc==4)+(dirc==2);
                     third_a{dirc}=current{third_r(dirc),third_c(dirc)};
                     if third_rot(dirc)~=5 %not selfreplace2
                         current_pos(third_r(dirc) , third_c(dirc),3)=mod(current_pos(third_r(dirc) , third_c(dirc),3)+third_rot(dirc),4);
                         third_track{dirc}=current_pos(third_r(dirc) , third_c(dirc),1:3);
                     end
                 end %if (forth_r(dirc)>0) & (forth_c(dirc)>0)
             else
                 % do not change
             end % if profit_do{choiceN}(n1)~=0
          end % if slook(dirc)==1
         end % for dirc
         if (first_r==1 && second_r~=1) || (first_r==16 && second_r~=16) || (first_c==1 && second_c~=1) || (first_c==16 && second_c~=16) || (first_r~=1 && second_r==1) || (first_r~=16 && second_r==16) || (first_c~=1 && second_c==1) || (first_c~=16 && second_c==16) 
           flagcancel=1;
         else
           flagcancel=0;
         end
         if flagcancel==0 | orientlock==0
             current{first_r,first_c}=second_a((second_rot==0)*(1:4)+(second_rot==1)*[4 1:3]  + (second_rot==2)*[3:4 1:2]  + (second_rot==3)*[2:4 1] ) ;
             current{second_r,second_c}=first_a((first_rot==0)*(1:4)+(first_rot==1)*[4 1:3]  + (first_rot==2)*[3:4 1:2]  + (first_rot==3)*[2:4 1] ) ;
             current_pos(first_r,first_c,1:3)=second_track;
             current_pos(second_r,second_c,1:3)=first_track;
             for dirc=1:4
               if (slook(dirc)==1) && (forth_r(dirc)>0) && (forth_c(dirc)>0)
                 if (third_r(dirc)==1 && forth_r(dirc)~=1) || (third_r(dirc)==16 && forth_r(dirc)~=16) || (third_c(dirc)==1 && forth_c(dirc)~=1) || (third_c(dirc)==16 && forth_c(dirc)~=16) || (third_r(dirc)~=1 && forth_r(dirc)==1) || (third_r(dirc)~=16 && forth_r(dirc)==16) || (third_c(dirc)~=1 && forth_c(dirc)==1) || (third_c(dirc)~=16 && forth_c(dirc)==16) 
                   flagcanceldir=1;
                 else
                   flagcanceldir=0;
                 end
                 if flagcanceldir==0 | orientlock==0
                     n1=8+(dirc-1)*4;
                     n2=9+(dirc-1)*4;
                     n3=10+(dirc-1)*4;
                     n4=11+(dirc-1)*4;
                     if profit_do{choiceN}(n1)~=0    % is zero if third-u/d/r/l is not changed
                             if third_rot(dirc)==5 %selfreplace2:  4 rotated to 3=4
                                current{forth_r(dirc),forth_c(dirc)}=forth_a{dirc}((forth_rot(dirc)==0)*(1:4)+(forth_rot(dirc)==1)*[4 1:3]  + (forth_rot(dirc)==2)*[3:4 1:2]  + (forth_rot(dirc)==3)*[2:4 1] ) ;
                                current_pos(forth_r(dirc),forth_c(dirc),1:3)=forth_track{dirc};
                             else
                                current{forth_r(dirc),forth_c(dirc)}=third_a{dirc}((third_rot(dirc)==0)*(1:4)+(third_rot(dirc)==1)*[4 1:3]  + (third_rot(dirc)==2)*[3:4 1:2]  + (third_rot(dirc)==3)*[2:4 1] ) ;
                                current{third_r(dirc),third_c(dirc)}=forth_a{dirc}((forth_rot(dirc)==0)*(1:4)+(forth_rot(dirc)==1)*[4 1:3]  + (forth_rot(dirc)==2)*[3:4 1:2]  + (forth_rot(dirc)==3)*[2:4 1] ) ;
                                current_pos(forth_r(dirc),forth_c(dirc),1:3)=third_track{dirc};
                                current_pos(third_r(dirc),third_c(dirc),1:3)=forth_track{dirc};
                             end
                     else
                         % do not change
                     end % if profit_do{choiceN}(n1)~=0
                 end %if flagcanceldir==0
              end % if (slook(dirc)==1) & (forth_r(dirc)>0) & (forth_c(dirc)>0)
             end % for dirc
         end %if flagcancel==0
     end % if profit_do{choiceN}(1)==1 
     
     end %if skiperfect==0
     
     temporary=current;
     
     % remove just replaced / rotated/ not replaced piece 1 from the picklist
     if skiperfect==0
         
     indNpicks=1:Npicks;
     for n=1:Npicks  
         if region==1
             pick=picklist(n,1:2);
         elseif region==2
             pick=picklist_border(n,1:2);
         elseif region==3
             pick=picklist_corner(n,1:2);
         end
         if (pick(1)==first_r) && (pick(2)==first_c)
              if region==1
                picklist=picklist(indNpicks(indNpicks~=n),1:2);
                picklistN=picklistN-1;
              elseif region==2
                picklist_border=picklist_border(indNpicks(indNpicks~=n),1:2);
                picklist_borderN=picklist_borderN-1;
              elseif region==3
                picklist_corner=picklist_corner(indNpicks(indNpicks~=n),1:2);
                picklist_cornerN=picklist_cornerN-1;
              end
              break
         end
     end %for n
     
     end % if skiperfect==0
     
     end % if lock(first...)==0
    end %of while
    
    
             score=0;
             for score_r=1:meshr-1
                 for score_c=1:meshc-1
                    score=score+( current{score_r,score_c}(3)==current{score_r+1,score_c}(1) ) + ( current{score_r,score_c}(2)==current{score_r,score_c+1}(4) ) ;
                 end
             end
             score_c=meshc;
             for score_r=1:meshr-1
                     score=score+ ( current{score_r,score_c}(3)==current{score_r+1,score_c}(1) ) ;
            end
             score_r=meshr;
             for score_c=1:meshc-1
                     score=score+ ( current{score_r,score_c}(2)==current{score_r,score_c+1}(4) )  ;
             end
    
    % track of pieces is in table current_pos
    % choices indicated the sum of choices with equal value that were met (not the total number of
    % choices, since we do not multiply the choices).
    
    ncore; %size of core reached
    ncorev(phase+1)=ncore;
    
    registerN=registerN+1;
    register_score(registerN)=score;
    register_ncore(registerN)=ncore;
    register_branches(registerN)=choices;
    
        if ncore>ncoremax
            ncoremax=ncore;
        end
    
        deltascore=score-scoreprev;
        scoreprev=score;
    
    if score>=scoremax
        scoremax=score;
        disp(sprintf('best_score=%d  ',scoremax));
        currentmax=current;
        current_posmax=current_pos;
    
        best_current=current;
    
    
    else
        current=current_back;
        currnet_pos=current_pos_back;
         temporary=current;
    end
    
    
    end % while score>=scoremax
    
    end % while score<480
    
    for ind_location=1:256
        row=floor((ind_location-1)/16)+1;
        col=ind_location-(row-1)*16;
        index=whichisClue(best_current{row,col},original);
        nnew_options(ind_location,:)=0;
        nnew_options(ind_location,index)=1;
    end

end % function




function options=external2native(op,original)

    translatenum=[0 1 9 17 5 13 2 10 18 6 14 3 11 19 7 15 4 12 20 8 16 21 22];
    options=zeros(256,1024);
    for row=1:16
        for col=1:16
            ind_location=(row-1)*16+col;
            obj=op(row,(col-1)*4+1:(col-1)*4+4);
            obj_trans=[0 0 0 0];
            for t=1:4
                obj_trans(t)=translatenum(obj(t)+1);
            end
            indchosen=whichisClue(obj_trans,original);
            if indchosen>0
                indstart=floor((indchosen-1)/4)*4+1;
                options(:,indstart:indstart+3)=0;
                options(ind_location,:)=0;
                options(ind_location,indchosen)=1;
            else
                disp('fault in file');
            end

        end
    end

end

function op=translate2external(piece,show_state,show_index)
    translatenum=[0 1 6 11 16 4 9 14 19 2 7 12 17 5 10 15 20 3 8 13 18 21 22];
    op=zeros(16,64);
    for row=1:16
        for col=1:16
            ind_location=(row-1)*16+col;
            ind=show_index(ind_location);
            piece1=1+floor((ind-1)/4);
            rot1=ind-4*(piece1-1);
            obj=[0 0 0 0];
            obj_trans=[0 0 0 0];
            for direction=1:4
                obj(direction)=piece(piece1,mod(direction-1-(rot1-1),4)+1);
                obj_trans(direction)=translatenum(obj(direction)+1);
            end
            op(row,(col-1)*4+1:(col-1)*4+4)=obj_trans;
        end
    end

end



function vector=any_neighbor_match(options,tree,indlocation)
        indf=1:1024;
        treesz=size(tree,3);
        row=floor((indlocation-1)/16)+1;
        col=indlocation-(row-1)*16;

        vector=zeros(size(indf));
        next_location_vector=zeros(1,4);
        for direction=1:4
            if direction==1
                next_row=row-1;
                next_col=col;
            elseif direction==2
                next_row=row;
                next_col=col+1;
            elseif direction==3
                next_row=row+1;
                next_col=col;
            elseif direction==4
                next_row=row;
                next_col=col-1;
            end
            if next_row==0 || next_row==16+1 || next_col==0 || next_col==16+1
                next_location=0;
                %next_sym=0;
            else
                next_location=(next_col-1)+(next_row-1)*16+1;
            end
            next_location_vector(direction)=next_location;
        end
        
        for direction=1:4
            nextlocation=next_location_vector(direction);
            if nextlocation>0
                rev_direction=mod(direction-1-2,4)+1;
                for ind=1:1024
                    if options(nextlocation,ind)>0
                        for pos=1:treesz
                            a=tree(ind,rev_direction,pos);
                            if a>0
                                vector(a)=1;
                            end
                        end
                    end
                end
            end

        end


end




function startSafeParpool(desiredNumWorkers)
    % Get local cluster
    c = parcluster('local');

    % Determine the maximum allowed by the cluster
    maxClusterWorkers = c.NumWorkers;

    % Get number of physical cores available
    hwInfo = evalc('feature(''numCores'')');
    numCores = str2double(regexp(hwInfo, '\d+', 'match', 'once'));

    % Choose the safest number of workers
    safeNumWorkers = min([desiredNumWorkers, maxClusterWorkers, numCores]);

    % Start or reuse pool
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(c, safeNumWorkers);
    else
        fprintf('Using existing pool with %d workers.\n', poolobj.NumWorkers);
    end
end


function [options,show_state,show_index,failblock]=form_block(threadno,streamkey,options,baseline_options,show_state,show_index,piece,baseline_tree,parallelstate,parallelForce,e2folder,compromize,isprob,h_number,slash_mode,ghost_options,flag_useHamiltonian,filen,pythonserverfile,OceanFolder,flag_snail)
%function [options,show_state,show_index,parallelstate,parallelForce,failblock]=form_block(threadno,options,baseline_options,show_state,show_index,piece,baseline_tree,parallelstate,parallelForce,e2folder,compromize,isprob,h_number,slash_mode,ghost_options)
    
 	bW=16;
	bS=16*16;
	flagloop=false;
    failblock=false;
    prob_options=options;
    extraline=h_number;

    %flag_useHamiltonian=input('Use Ising Hamiltonian? (0-no, 1-yes): ');
    if flag_useHamiltonian
        timestamp = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
        writematrix('stop',[OceanFolder '/command2server' timestamp '.txt']);
        writematrix('0',[OceanFolder '/replyfromserver' timestamp '.txt']);
        outnumbers_count=1;
        load(filen); %load (NodeSpin), mapnodes_ind, mapnodes_loc, options (not needed)
		%consider 
		%pbsScript = 'ocean_job.pbs';
		%system(['qsub -v TIMESTAMP=' timestamp ' ' pbsScript]);

		commandext = "bash -c ""module load mathlib mkl/2023 && module load python/2023 && source " + OceanFolder +"/venv/bin/activate && nohup " + OceanFolder + "/venv/bin/python "+ pythonserverfile + " "+ timestamp+ " > output.log 2>&1 & "" &";
		disp(commandext);
		[status, cmdout] = system(commandext);
		disp(status);
		disp('Output:');
		disp(cmdout);

    end
    options=prob_options;


    sieve_level=2;%input('Sieve level (1-light, 2-serious): ');
    if extraline==8
        flag_skip=true;
    else
        flag_skip=false;
    end

    block_pos=zeros(1,2);
    posind=0;

    if flag_snail==false
        for t=1:extraline
            r=t;
            v_c=1:extraline;
            if mod(r,2)==0
                v_c=v_c(end:-1:1);
            end
            for c=v_c
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
        for t=extraline+1:16-extraline
            r=t;
            v_c=1:extraline;
            %if mod(r,2)==0
            %    v_c=v_c(end:-1:1);
            %end
            for c=v_c
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
        for t=16+1-extraline:16
            r=t;
            v_c=1:extraline;
            if mod(r,2)==0
                v_c=v_c(end:-1:1);
            end
            for c=v_c
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
        for t=extraline+1:16-extraline
            c=t;
            v_r=16+1-extraline:16;
            if true %mod(c,2)==0
                v_r=v_r(end:-1:1);
            end
            for r=v_r
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
        for t=16+1-extraline:16
            c=t;
            v_r=16+1-extraline:16;
            if mod(c,2)==0
                v_r=v_r(end:-1:1);
            end
            for r=v_r
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
        for t=16-extraline:-1:extraline+1
            r=t;
            v_c=16+1-extraline:16;
            if true %mod(r,2)==0
                v_c=v_c(end:-1:1);
            end
            for c=v_c
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
        for t=extraline:-1:1
            r=t;
            v_c=16+1-extraline:16;
            if mod(r,2)==0
                v_c=v_c(end:-1:1);
            end
            for c=v_c
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
        t_vectorcon=[extraline+1:5 (16-extraline):-1:12 6 11 7 10 8 9];
        t_vectorcon(t_vectorcon<=extraline)=[];
        t_vectorcon(t_vectorcon>=17-extraline)=[];
        for t=t_vectorcon   %5:12, but intermittantly to make the stitch working with the most recent choices, which is much more efficient
            c=t;
            v_r=1:extraline;
            %if mod(c,2)==0
            %    v_r=v_r(end:-1:1);
            %end
            for r=v_r
                posind=posind+1;
                block_pos(posind,1:2)=[r c];
            end
        end
    else %now if flag_snail==true
        flag_reverse=(rand>0.5);
        for layerno=1:extraline

            if flag_reverse==false
                c=layerno;
                for r=layerno:17-layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
                r=17-layerno;
                for c=1+layerno:16-layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
                c=17-layerno;
                for r=17-layerno:-1:layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
                r=layerno;
                for c=16-layerno:-1:1+layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
            else
                c=layerno;
                for r=17-layerno:-1:layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
                r=layerno;
                for c=1+layerno:16-layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
                c=17-layerno;
                for r=layerno:17-layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
                r=17-layerno;
                for c=16-layerno:-1:1+layerno
                    posind=posind+1;
                    block_pos(posind,1:2)=[r c];
                end
            end
            
        end

    end %if flag_snail==false

    %skip positions that are empty in ghost, thinking it is a sign of unfavorable place for nucleation 
    pos=0;
    if ~isempty(ghost_options)
        for pos_tag=1:posind
            pos=pos+1;
            row=block_pos(pos,1);
            col=block_pos(pos,2);
            indlocation=(row-1)*16+col;
            if sum(ghost_options(indlocation,:))==0
                block_pos(pos:posind-1,1:2)=block_pos(pos+1:posind,1:2);
                posind=posind-1;
                pos=pos-1;
            end
        end

    end

    indf=1:1024;
    if ~isempty(ghost_options)
        ghost_imglin=sum(ghost_options>0,2);
        tshow_state=1*(ghost_imglin==1);
        tshow_index=zeros(size(tshow_state));
        for tindlocation=1:256
            if tshow_state(tindlocation)==1
                tshow_index(tindlocation)=min(indf(ghost_options(tindlocation,:)>0));
            end
        end
        ghost_tree=update_tree(baseline_tree,tshow_state,tshow_index);
    end



    maxpos=posind;
    disp(sprintf('maxpos=%d',maxpos));
    firstopen_posind=-1;
    shuffled_choices=zeros(maxpos,100);
    count_choices=zeros(maxpos,1);
    posof_choices=ones(maxpos,1);
    for posind=1:maxpos
        row=block_pos(posind,1);
        col=block_pos(posind,2);
        indlocation=(row-1)*16+col;
        if show_state(indlocation)==1
            count_choices(posind)=-1;
            shuffled_choices(posind)=show_index(indlocation);
        elseif firstopen_posind==-1
            firstopen_posind=posind;
        end
    end
    indf=1:1024;
    base_show_state=show_state;
    base_show_index=show_index;
    base_options=options;
    active_posind=1;
    maxfill_posind=0;
    total_time=0;
    flag_ended=false;
    flag_needrefresh=true;
    tic;
    while sum(count_choices~=0)<maxpos && ~(flag_ended && flag_skip) % && (sum(count_choices~=0)<maxpos-2 || toc<3600*6)
        for posind=active_posind:maxpos
            row=block_pos(posind,1);
            col=block_pos(posind,2);
            indlocation=(row-1)*16+col;
            flag_needrefresh=flag_needrefresh || (count_choices(posind)==0);

            score_nostop=1;
            if flag_needrefresh
                %%%%  board to limit further choices
                keep_state5=(show_state==5);
                show_state= 1*(keep_state5 | base_show_state>0); %use only chosen options or baseline options, not including predictions from previous
                show_index(base_show_state>0)=base_show_index(base_show_state>0);
                show_state(indlocation)=0; %so will take from base (state before starting h), fresh
                keep_state5(indlocation)=0;
                for tindlocation=1:256
                    if show_state(tindlocation)==0 
                        options(tindlocation,:)=base_options(tindlocation,:);
                    end
                end
                [options,show_state,show_index]=force_to_options(options,base_options,show_state,show_index);
                [options,show_state,show_index]=force_to_home(options,base_options,show_state,show_index);
                show_state(keep_state5)=5;
                show_index(show_state==0)=0;
                tree=baseline_tree;
                tree=update_tree(tree,show_state,show_index);
                if isprob==0
                    [options,score_nostop]=new_options(compromize,options,tree); %predict (eliminate) choices
                else
                    %[options,score_nostop]=new_options_prob(compromize,options,tree); %predict (eliminate) choices
                end
                [options,show_state,show_index]=force_to_options(options,base_options,show_state,show_index);
                [options,show_state,show_index]=force_to_home(options,base_options,show_state,show_index);
                show_state(keep_state5)=5;
                show_index(show_state==0)=0;
                flag_needrefresh=false;
            end %if flag_needrefresh



            if count_choices(posind)==0
                lenvect=0;
                if all(score_nostop>0 || compromize) %if board correct start to find options, otherwise lenvect=0, go backward
                    indf_possible=indf(options(indlocation,indf)>0);

                    if length(indf_possible)==1
                         index_solved_vector=indf(options(indlocation,:)>0);
                    else

                        for ind=indf_possible
                            if options(indlocation,ind)>0
                                piece1=1+floor((ind-1)/4);
                                rot1=ind-4*(piece1-1);
                                flag=true;
                                for direction=1:4
                                    sym=piece(piece1,mod(direction-1-(rot1-1),4)+1);
                                    if direction==1
                                        next_row=row-1;
                                        next_col=col;
                                    elseif direction==2
                                        next_row=row;
                                        next_col=col+1;
                                    elseif direction==3
                                        next_row=row+1;
                                        next_col=col;
                                    elseif direction==4
                                        next_row=row;
                                        next_col=col-1;
                                    end
                                    if next_row==0 || next_row==17 || next_col==0 || next_col==17
                                        next_sym=0;
                                    else
                                        next_location=(next_col-1)+(next_row-1)*16+1;
                                        indf=1:1024;
                                        next_index_solved_vector=indf(options(next_location,:)>0);
                                        if show_state(next_location)>0 && sum(next_index_solved_vector==show_index(next_location))>0 %user selected one choice
                                            index_piece_v=show_index(next_location);
                                        else
                                            index_piece_v=indf(options(next_location,:)>0); %take all the options of neighbor pieces
                                        end
                                        piece1_v=1+floor((index_piece_v-1)/4);
                                        rot1_v=index_piece_v-4*(piece1_v-1);
                                        next_sym=piece(sub2ind(size(piece),piece1_v,mod(direction-1+2-(rot1_v-1),4)+1));
                                    end
                                    flag=flag && sum(next_sym==sym)>0;
                                end
                                if flag==false
                                    options(indlocation,ind)=0;
                                end
                            end
                        end
                        index_solved_vector=indf(options(indlocation,:)>0);
                        indu=1:256;
                        %remove choices of pieces already shown in other positions
                        for tspos=1:length(index_solved_vector)
                            indt=index_solved_vector(tspos);
                            indt_rot1=floor((indt-1)/4)*4+1;
                            for indt_rot=indt_rot1:indt_rot1+3
                                if sum(show_index(indu~=indlocation)==indt_rot & show_state(indu~=indlocation)>0)>0
                                    options(indlocation,indt)=0;
                                    index_solved_vector(tspos)=-1;
                                end
                            end
                        end
                        index_solved_vector=index_solved_vector(index_solved_vector>0);
                        %%%%%%%%%%%%%%%@@@@@@@@@@@@###########$$$%%%%%%%%%%%%%%%%%%%%%%
                        %index_solved_vector=4*(indlocation-1)+1; %%%% overwrite the known options for virtual puzzles

%Copy to the end of python server C:\Users\seifer\PycharmProjects\Ocean\HamiltonianServer.py
%{
import time
import sys
time_str = sys.argv[1]
command_file = f"/gpfs01/home/seifer/eternity30/Ocean/command2server{time_str}.txt"
reply_file = f"/gpfs01/home/seifer/eternity30/Ocean/replyfromserver{time_str}.txt"
print('Ready to take tasks')
while True:
    try:
        # Read the command file
        with open(command_file, "r") as f:
            content = f.read().strip()
        if content.lower() == "run":
            filename2=f"/gpfs01/home/seifer/eternity30/Ocean/QA/spin{time_str}.csv"
            spin_array2 = np.loadtxt(filename2, dtype=int, delimiter=',')
            spin_dict2 = dict(zip(bqm.variables, spin_array2))
            energy2 = bqm.energy(spin_dict2)
            print("Hamiltonian value (energy) known solution :", energy2)
            with open(command_file, "w") as f:
                f.write(str(energy2))
            with open(reply_file, "a") as f:
                f.write(f"{energy2}\n")
        time.sleep(1)
    except Exception as e:
        print(f"Error: {e}")
        time.sleep(1)
%}

                        toptions=options;
                        timglin=sum(toptions>0,2);
                        tshow_state=1*(timglin==1);
                        tshow_index=zeros(size(tshow_state));
                        for tindlocation=1:bS
                            if tshow_state(tindlocation)==1
                                tshow_index(tindlocation)=min(indf(toptions(tindlocation,:)>0));
                            end
                        end
						
						lenvect=length(index_solved_vector);
                        score_nostop_vector=ones(1,lenvect);
                        if lenvect>=1 && flag_useHamiltonian
                            %Here the python program should be already
                            %running, just communicate with files
                            %use: indchosen, indlocation, options, 
                            flag_startfromDWAVEsolution=false;
                            if flag_startfromDWAVEsolution
                                bestspin=NodeSpin;
                            else
                                bestspin=1*ones(size(mapnodes_loc'));
                                for nd=1:length(mapnodes_loc)
                                    tloc=mapnodes_loc(nd);
                                    tind=mapnodes_ind(nd);
                                    bestspin(nd)=1-2*toptions(tloc,tind);
                                end
                            end
                            indg=1:length(bestspin);
                            for tloc=1:bS
                                if tshow_state(tloc)>0
                                    ind_actual=tshow_index(tloc);
                                    bestspin(indg(mapnodes_ind==ind_actual))=1;
                                    bestspin(indg(mapnodes_loc==tloc))=1;
                                    bestspin(indg(mapnodes_loc==tloc & mapnodes_ind==ind_actual))=-1;
                                end
                            end
                            bestspin0=bestspin;
                            Hamiltonian=0;
                            for testchoice=1:lenvect
                               bestspin=bestspin0;
                               indchosen=index_solved_vector(testchoice);
                               bestspin(indg(mapnodes_ind==indchosen))=1;
                               bestspin(indg(mapnodes_loc==indlocation))=1;
                               bestspin(indg(mapnodes_loc==indlocation & mapnodes_ind==indchosen))=-1;
                               filenamebest=[OceanFolder '/QA/spin' timestamp '.csv'];
                               writematrix(bestspin,filenamebest);
                               writematrix('run',[OceanFolder '/command2server' timestamp '.txt']);
                               while true
                                    outnumbers=readmatrix([OceanFolder '/replyfromserver' timestamp '.txt']);
                                    if length(outnumbers)>outnumbers_count
                                        outnumbers_count=length(outnumbers);
                                        Hamiltonian=outnumbers(outnumbers_count);
                                        break;
                                    end
									pause(3);
                               end
                               score_nostop_vector(testchoice)=-Hamiltonian;  %Use in negative: good if score high
                            end
                        else %if not flag_Hamiltonian
                                %score_nostop_vector=randperm(streamkey,lenvect); %shuffle vector streamkey
                                if lenvect>1
                                    for t=1:round(rand*10)
                                        randperm(lenvect);
                                    end
                                end
                                index_solved_vector=index_solved_vector(randperm(streamkey,lenvect)); %shuffle vector

                        end   %lenvect>=1 && flag_useHamiltonian

                        %sort by score (after vector already shuffled)
                        %score_nostop_vector(score_nostop_vector>0 & score_nostop_vector<=112)=1;
                        %score_nostop_vector=ceil(score_nostop_vector/1)*1;  %make it blind to score differences of up to 4, since the path to solution does not go directly at maximum score 
                        [~, index_vector]=sort(score_nostop_vector,'descend');
                        index_solved_vector=index_solved_vector(index_vector);%apply score order on prioirity of piece-index to choose first
                        index_solved_vector=index_solved_vector(index_solved_vector>0);

                    end %end of selection between one option and more

                    lenvect=length(index_solved_vector);
                    if lenvect>0
                        active_posind=posind;
                        count_choices(posind)=lenvect;
                        posof_choices(posind)=1;
                        if slash_mode && count_choices(posind)>1 %define preferred order according to dwave answers
                            deepfill=false;
                            problemsize=4000;
                            [score_dwave, ~]=dwave_solve(indlocation,options,index_solved_vector, show_state,show_index, piece,deepfill,problemsize,0 );
                            [~, index_vector]=sort(score_dwave,'descend');
                            index_solved_vector=index_solved_vector(index_vector);%apply score order on prioirity of piece-index to choose first
                        elseif ~isempty(ghost_options) && count_choices(posind)>1
                            if sum(ghost_imglin)>100
                                vector_flags_from_neighbors=any_neighbor_match(ghost_options,ghost_tree,indlocation);
                                score_ghost=(ghost_options(indlocation,index_solved_vector)>0) | vector_flags_from_neighbors(index_solved_vector);
                            else
                                score_ghost=(ghost_options(indlocation,index_solved_vector)>0);
                            end
                            [~, index_vector]=sort(score_ghost,'descend'); %reorder already ordered choice list, to have the group from slash mode solution first (ghost_options, as candidate)
                            index_solved_vector=index_solved_vector(index_vector);
                        end
                        shuffled_choices(posind,1:lenvect)=index_solved_vector;
                        show_index(indlocation)=shuffled_choices(posind,posof_choices(posind));
                        show_state(indlocation)=5; %mark as planned choices
                        if posind>maxfill_posind
                            maxfill_posind=posind;
                            %fprintf(' %d[%g] ',posind,threadno);
                            best_options=options;
                            best_show_state=show_state;
                            best_show_index=show_index;
                            total_time=total_time+toc;
                            best_score_show=sum(show_state>0);
                            disp(sprintf(' %d[%g] =%d',posind,threadno,best_score_show));
                            %planv=0; %for virtual puzzle
                            %if planv==0
                            required_options=zeros(256,1024);
                            for t=1:256
                                 required_options(t,(t-1)*4+1)=1;
                            end
                            firstcheck_actual_sol=sum(sum((options>0.5).*(required_options>0.5),2));
                            %    fprintf('(s=%d), ',firstcheck_actual_sol);
                            %end
                            if posind>=maxpos-1 || best_score_show>=160
                                try
                                    filesave=sprintf('%s//goingon39_th%g_%g_%g_known%g.mat',e2folder,threadno,posind,best_score_show,firstcheck_actual_sol);
                                    save(filesave,'best_options','shuffled_choices','posof_choices','count_choices','posind','show_state','show_index','-mat');
                                catch
                                end
                            end
                        end


                        lines=readlines(sprintf('%s//break.txt',e2folder));
                        if contains(lines(1),'end') 
                            %disp('User requested to exit');
                            flagloop=false;
                            failblock=true;
                            disp('break')
                            %options=best_options;
                            %show_state=best_show_state;
                            %show_index=best_show_index;
                            return;
                        end
                    end
                end   %end of test for available choices
                if lenvect==0
                    options(indlocation,:)=prob_options(indlocation,:);%baseline_options(indlocation,:);
                    show_state(indlocation)=0;
                    show_index(indlocation)=0;
                    flag_needrefresh=true;
                     
                    if flag_skip %dont stop, skip
                        continue;
                    else
                        if posind==firstopen_posind
                            disp('No match with initial condition');
                            failblock=true;
                            return;
                        end
                        for t=posind-1:-1:1
                            if count_choices(t)>=0 
                                active_posind=t;
                                break;
                            end
                        end
                        break;
                    end
                end
            elseif posof_choices(posind)<count_choices(posind)
                %advance choice by 1
                active_posind=posind;
                posof_choices(posind)=posof_choices(posind)+1;
                show_index(indlocation)=shuffled_choices(posind,posof_choices(posind));
                show_state(indlocation)=5; %to be sure

            elseif posof_choices(posind)==count_choices(posind)
                show_state(indlocation)=0;
                show_index(indlocation)=0;
                options(indlocation,:)=prob_options(indlocation,:);%baseline_options(indlocation,:);
                posof_choices(posind)=1;
                count_choices(posind)=0;
                flag_needrefresh=true;

               
                for t=posind-1:-1:1
                    if count_choices(t)>=0 
                        active_posind=t;
                        break;
                    elseif t==1
                        disp('No solution');
                        failblock=true;
                        return;
                    end
                end
                break;
            end
            flagloop=false;
        end

        flag_ended=true;

    end
    try
        filesave=sprintf('%s//blockready_th%g_%s.mat',e2folder,threadno,char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss')));
        save(filesave,'options','-mat');
    catch
    end
    disp('Block ready');


end


function [score_dwave, options] = dwave_solve(indlocation,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,filenumber_ready)
    global PaidDwaveService;
	global OceanFolder;
    quantum_computation_paid=PaidDwaveService;
    bW=16; %board width
    bS=bW*bW;%size of board
    b4S=bS*4;


    if length(shuffled_choices_pick)==0
        disp('No fit to this position, try another one');
        indpiece=[];
        return
    end
    max_nodes=problemsize;

    if deepfill==false
        use_hybrid_solver=true;
        %max_couples=24; %Usal maximum 24,  Advantage_system6.4 can use 34
        if use_hybrid_solver
            max_couples=3000; %hybrid case
        else
            max_couples=23; %Usal maximum 24,  Advantage_system6.4 can use 34
            %max_couples=33; %Usal maximum 24,  Advantage_system6.4 can use 34
        end
    elseif deepfill==true
        use_hybrid_solver=true;
        if use_hybrid_solver
            max_couples=3000; 
        else
            max_couples=23; 
        end
    end
    disp(sprintf('Preparing to use dwave with up to %g spins.',max_nodes))

    if filenumber_ready==0

        mapnodes_ind=zeros(1,max_nodes);
        mapnodes_loc=zeros(1,max_nodes);
        mapnodes_notavariable=zeros(1,max_nodes);
        mapnodes_couplers_pos=zeros(1,max_nodes); %slot number
        mapnodes_couplers_count=zeros(1,max_nodes);
        mapnodes_couplers_with=zeros(max_couples,max_nodes);
        mapnodes_couplers_weight=zeros(max_couples,max_nodes);
        mapnodes_magnetic=zeros(1,max_nodes+1);
        missing_coupling_count=zeros(1,max_nodes);

        QAcounter=readmatrix([OceanFolder '/venv/bin/python ' '/QA/counter.txt']); %In windows environment variables define: Ocean = C:\Users\seifer\PycharmProjects\Ocean
        QAcounter=QAcounter+1;
        fileID = fopen([OceanFolder '/' sprintf('ask_d_wave%g.py',QAcounter)],'w');
        fprintf(fileID,'import numpy as np\n');
        fprintf(fileID,'import networkx as nx\n');
        fprintf(fileID,'import matplotlib.pyplot as plt\n');
        fprintf(fileID,'import dimod\n');
        fprintf(fileID,'import csv\n');
        fprintf(fileID,'import os\n');
        fprintf(fileID,'from dwave.system import EmbeddingComposite\n');
        fprintf(fileID,'from dwave.system.samplers import DWaveSampler\n');
        fprintf(fileID,'from dwave.system import LeapHybridSampler\n'); %binary quantum solver: hybrid classical-quantum 
        fprintf(fileID,'from dwave.cloud.exceptions import *\n');


        fprintf(fileID,'G=nx.Graph()\n');

        %add nodes for main selection in question
        fprintf(fileID,sprintf('#Nodes in question: %d \n',length(shuffled_choices_pick)));
        question_size=length(shuffled_choices_pick);
        for node_number=1:min(length(shuffled_choices_pick),max_nodes)
            mapnodes_ind(node_number)=shuffled_choices_pick(node_number);
            mapnodes_loc(node_number)=indlocation;
        end
        CrossCheck_flag=false(bS,bS,4);
        failed_coupling_counter=0;

        indf=1:b4S;
        flag_random=false;

 
        if ~flag_random
            stage_vector=indlocation;
            for row=1:(bW)
                for col=1:(bW)
                    indlocationt=(col-1)+(row-1)*bW+1;
                    if indlocationt~=indlocation
                        ind_vect=indf(options(indlocationt,:)>0);
                        for ind=ind_vect
                            if node_number<max_nodes
                                node_number=node_number+1;
                                mapnodes_ind(node_number)=ind;
                                mapnodes_loc(node_number)=indlocationt;
                            end
                        end
                        stage_vector=[stage_vector indlocationt];
                    end
                end
            end
        else
            stage_vector=1:max_nodes*4;
        end
        
        for stage=stage_vector
            if flag_random
                nodef=1:node_number;
                indlocation_vect=mapnodes_loc(nodef);
                nodef_avail=nodef;%nodef(show_state(indlocation_vect)==0);
                node_randomly=nodef_avail(floor(rand(1)*length(nodef_avail))+1);
                indlocation=mapnodes_loc(node_randomly);
                node_choices=nodef(mapnodes_loc==indlocation);
                ind_choices=mapnodes_ind(node_choices);
            elseif ~flag_random
                indlocation=stage;
                nodef=1:node_number;
                node_choices=nodef(mapnodes_loc==indlocation);
                ind_choices=mapnodes_ind(node_choices);
            end

            group_table_thislocation_nodes=zeros(4,40,200);
            group_table_nextlocation_nodes=zeros(4,40,200);
            group_sym_thislocation=zeros(4,40);
            %group_sym_nextlocation=zeros(4,40);%same slot as thislocation
            group_nodecounter_thislocation=zeros(4,40);
            group_symcounter_thislocation=zeros(4,1);
            group_nodecounter_nextlocation=zeros(4,40);
            %group_symcounter_nextlocation=zeros(4,1); %the same as for thislocation
    
            next_location_vector=zeros(4,1);

            row=floor((indlocation-1)/bW)+1;
            col=indlocation-(row-1)*bW;
            for direction=1:4
                if direction==1
                    next_row=row-1;
                    next_col=col;
                elseif direction==2
                    next_row=row;
                    next_col=col+1;
                elseif direction==3
                    next_row=row+1;
                    next_col=col;
                elseif direction==4
                    next_row=row;
                    next_col=col-1;
                end
                if next_row==0 || next_row==bW+1 || next_col==0 || next_col==bW+1
                    next_location=0;
                    %next_sym=0;
                else
                    next_location=(next_col-1)+(next_row-1)*bW+1;
                end
                next_location_vector(direction)=next_location;
            end


            for direction=1:4
                for place_in_choices=1:length(ind_choices)
                    ind=ind_choices(place_in_choices);
                    node_of_ind=node_choices(place_in_choices);
                    piece1=1+floor((ind-1)/4);
                    rot1=ind-4*(piece1-1);
                    index_piece_next_match=zeros(1,0); %initialize as empty
                    sym=piece(piece1,mod(direction-1-(rot1-1),4)+1);
                    next_location=next_location_vector(direction);
                    if next_location>0 && ~CrossCheck_flag(indlocation,next_location,direction) && sym>0
                        index_piece_v=indf(options(next_location,:)>0); %take all the options of neighbor pieces
                        %%index_piece_v=index_piece_v(index_piece_v~=ind);
                        piece1_v=1+floor((index_piece_v-1)/4);
                        rot1_v=index_piece_v-4*(piece1_v-1);
                        next_sym=piece(sub2ind(size(piece),piece1_v,mod(direction-1+2-(rot1_v-1),4)+1));
                        index_piece_next_match=index_piece_v(next_sym==sym);
                    else
                        continue;
                    end
                    %Important: skip unfruitfull directions that consume the node count
                    %if length(index_piece_next_match)>40
                    %    continue;
                    %end
                    if  sym>0   %~isempty(index_piece_next_match) &&
                        flag_newsym=true;
                        for testslot=1:group_symcounter_thislocation(direction)
                            if sym==group_sym_thislocation(direction,testslot)
                                flag_newsym=false;
                                thisslot=testslot;
                                group_nodecounter_thislocation(direction,testslot)=group_nodecounter_thislocation(direction,testslot)+1;
                                group_table_thislocation_nodes(direction,testslot,group_nodecounter_thislocation(direction,testslot))=node_of_ind;
                                break;
                            end
                        end
                        if flag_newsym
                            thisslot=group_symcounter_thislocation(direction)+1;
                            group_symcounter_thislocation(direction)=thisslot;
                            group_sym_thislocation(direction,thisslot)=sym;
                            group_nodecounter_thislocation(direction,thisslot)=group_nodecounter_thislocation(direction,thisslot)+1;
                            group_table_thislocation_nodes(direction,thisslot,group_nodecounter_thislocation(direction,thisslot))=node_of_ind;
                        end
                        %nxdirection=mod(direction-1+2,4)+1; NOT NEEDED:
                        %THE DIRECTION IS A LABEL OF THE NODES RELATED TO
                        %TWO ADJACENT LOCATIONS, FROM THE FIRST NODE
                        %PESPECTIVE
                        %%nextslot=group_symcounter_nextlocation(direction)+1;
                        %%group_symcounter_nextlocation(direction)=nextslot;
                        %%group_sym_nextlocation(direction,nextslot)=sym;
                        for temp=1:length(index_piece_next_match)
                            temp_index=index_piece_next_match(temp);
                            already_innextlocation=false;
                            for tempos=1:group_nodecounter_nextlocation(direction,thisslot)
                                temp_node=group_table_nextlocation_nodes(direction,thisslot,tempos);
                                if mapnodes_ind(temp_node)==temp_index
                                    already_innextlocation=true;
                                    break;
                                end
                            end
                            if already_innextlocation
                                continue;
                            end
                            do_newnode=true;
                            for temp2=1:node_number
                                if mapnodes_loc(temp2)==next_location
                                    if temp2~=node_of_ind
                                        if mapnodes_ind(temp2)==temp_index
                                            nodepos=group_nodecounter_nextlocation(direction,thisslot)+1;
                                            group_nodecounter_nextlocation(direction,thisslot)=nodepos;
                                            group_table_nextlocation_nodes(direction,thisslot,nodepos)=temp2;
                                            do_newnode=false;
                                        end
                                    end
                                end
                            end
                            if node_number< max_nodes && do_newnode
                                node_number=node_number+1;
                                %fprintf(fileID,'G.add_node(%g)\n',node_number);done after
                                mapnodes_ind(node_number)=index_piece_next_match(temp);
                                mapnodes_loc(node_number)=next_location;

                                nodepos=group_nodecounter_nextlocation(direction,thisslot)+1;
                                group_nodecounter_nextlocation(direction,thisslot)=nodepos;
                                group_table_nextlocation_nodes(direction,thisslot,nodepos)=node_number;
                            elseif node_number== max_nodes && do_newnode
                                failed_coupling_counter=failed_coupling_counter+1;
                                % node_number=prev_node_number;
                                % max_nodes=node_number;  %do not alow to add more, reject adding couplers to the added spins (to ignore them).
                            end
                        end


                    end

                end %for place_in_choices
            end %for direction=1:4

            for direction=1:4
                next_location=next_location_vector(direction);
                if next_location>0
                    CrossCheck_flag(indlocation,next_location,direction)=true; %prevent adding this couple again
                    CrossCheck_flag(next_location,indlocation,mod(direction-1+2,4)+1)=true; %prevent adding this couple again
                end
            end


            for dir=1:4
                for indsym1=1:group_symcounter_thislocation(dir)
                    sym1=group_sym_thislocation(dir,indsym1);
                    sym2=sym1;
                    indsym2=indsym1;
                    %the symbol is stored in the first location only
                    if sym1==sym2
                        A=group_nodecounter_thislocation(dir,indsym1);
                        B=group_nodecounter_nextlocation(dir,indsym2);
                        if A>max_couples || B>max_couples
                            A=0;
                            B=0;
                        end
                        for indnode1=1:A
                            node1=group_table_thislocation_nodes(dir,indsym1,indnode1);
                            if mapnodes_couplers_count(node1)<max_couples
                                mapnodes_magnetic(node1)=mapnodes_magnetic(node1)+0.5*(B-A);
                            end
                        end
                        for indnode2=1:B
                            node2=group_table_nextlocation_nodes(dir,indsym2,indnode2);
                            if mapnodes_couplers_count(node2)<max_couples
                                mapnodes_magnetic(node2)=mapnodes_magnetic(node2)+0.5*(A-B);
                            end
                        end
                        %case one of group A and one from B
                        for indnode1=1:A
                            for indnode2=1:B
                                node1=group_table_thislocation_nodes(dir,indsym1,indnode1);
                                node2=group_table_nextlocation_nodes(dir,indsym2,indnode2);
                                if node2<node1
                                    node2=group_table_thislocation_nodes(dir,indsym1,indnode1);
                                    node1=group_table_nextlocation_nodes(dir,indsym2,indnode2);
                                end
                                if node1<node2

                                    flag_needadd=true;
                                    for cpl_pos=1:mapnodes_couplers_pos(node1)
                                        if mapnodes_couplers_with(cpl_pos,node1)==node2
                                            mapnodes_couplers_weight(cpl_pos,node1)=mapnodes_couplers_weight(cpl_pos,node1)-0.5;
                                            flag_needadd=false;
                                            break;
                                        end
                                    end
                                    if flag_needadd
                                        if mapnodes_couplers_count(node1)<max_couples && mapnodes_couplers_count(node2)<max_couples
                                            mapnodes_couplers_count(node1)=mapnodes_couplers_count(node1)+1;
                                            mapnodes_couplers_count(node2)=mapnodes_couplers_count(node2)+1;
                                            cpos=mapnodes_couplers_pos(node1)+1;
                                            mapnodes_couplers_pos(node1)=cpos;
                                            mapnodes_couplers_with(cpos,node1)=node2;
                                            mapnodes_couplers_weight(cpos,node1)=-0.5;
                                        else
                                            missing_coupling_count(node1)=missing_coupling_count(node1)+1;
                                            missing_coupling_count(node2)=missing_coupling_count(node2)+1;
                                        end
                                    end
                                end
                            end
                        end
                        %case one of group A and one from A
                        for indnode1=1:A
                            for indnode2=1:A
                                node1=group_table_thislocation_nodes(dir,indsym1,indnode1);
                                node2=group_table_thislocation_nodes(dir,indsym1,indnode2);
                                if node1<node2
                                    flag_needadd=true;
                                    for cpl_pos=1:mapnodes_couplers_pos(node1)
                                        if mapnodes_couplers_with(cpl_pos,node1)==node2
                                            mapnodes_couplers_weight(cpl_pos,node1)=mapnodes_couplers_weight(cpl_pos,node1)+0.5;
                                            flag_needadd=false;
                                            break;
                                        end
                                    end
                                    if flag_needadd
                                        if mapnodes_couplers_count(node1)<max_couples && mapnodes_couplers_count(node2)<max_couples
                                            mapnodes_couplers_count(node1)=mapnodes_couplers_count(node1)+1;
                                            mapnodes_couplers_count(node2)=mapnodes_couplers_count(node2)+1;
                                            cpos=mapnodes_couplers_pos(node1)+1;
                                            mapnodes_couplers_pos(node1)=cpos;
                                            mapnodes_couplers_with(cpos,node1)=node2;
                                            mapnodes_couplers_weight(cpos,node1)=+0.5;
                                        else
                                            missing_coupling_count(node1)=missing_coupling_count(node1)+1;
                                            missing_coupling_count(node2)=missing_coupling_count(node2)+1;
                                        end
                                    end
                                end
                            end
                        end
                        %case one of group B and one from B
                        for indnode1=1:B
                            for indnode2=1:B
                                node1=group_table_nextlocation_nodes(dir,indsym2,indnode1);
                                node2=group_table_nextlocation_nodes(dir,indsym2,indnode2);
                                if node1<node2
                                    flag_needadd=true;
                                    for cpl_pos=1:mapnodes_couplers_pos(node1)
                                        if mapnodes_couplers_with(cpl_pos,node1)==node2
                                            mapnodes_couplers_weight(cpl_pos,node1)=mapnodes_couplers_weight(cpl_pos,node1)+0.5;
                                            flag_needadd=false;
                                            break;
                                        end
                                    end
                                    if flag_needadd
                                        if mapnodes_couplers_count(node1)<max_couples && mapnodes_couplers_count(node2)<max_couples
                                            mapnodes_couplers_count(node1)=mapnodes_couplers_count(node1)+1;
                                            mapnodes_couplers_count(node2)=mapnodes_couplers_count(node2)+1;
                                            cpos=mapnodes_couplers_pos(node1)+1;
                                            mapnodes_couplers_pos(node1)=cpos;
                                            mapnodes_couplers_with(cpos,node1)=node2;
                                            mapnodes_couplers_weight(cpos,node1)=+0.5;
                                        else
                                            missing_coupling_count(node1)=missing_coupling_count(node1)+1;
                                            missing_coupling_count(node2)=missing_coupling_count(node2)+1;
                                        end
                                    end
                                end
                            end
                        end


                    end
  
                end
            end


        end % for stage=


        disp(sprintf('Number of rejected coupling due to slot limitation: %g',failed_coupling_counter));


        for temps_node=1:node_number
            temps_loc=mapnodes_loc(temps_node);
            if show_state(temps_loc)>0
                mapnodes_magnetic(temps_node)=max(1,mapnodes_magnetic(temps_node)); %bias to -1 to choose fixed ones naturally
                mapnodes_notavariable(temps_node)=1;

            end
        end

        
        
        %Attraction to one spin down in each location (so the preferred value of sum
        %of spins is the number of spins minus 2). Apply only if coupling
        %is described fully because otherwise I preffer all spins -1.
        ancil_number=0;  %Count of added accessory nodes based on ancilla technique
        Gbase=0.0125;
        nodef=1:node_number;
        for ind_loc_temp=1:bS*(Gbase>0)
            if show_state(ind_loc_temp)==0
                node_list_temp=sort(nodef(mapnodes_loc(nodef)==ind_loc_temp & mapnodes_ind(nodef)>0)) ;
                M=length(node_list_temp);
                if M>0 && M<=600 %upper limit due to memory restriction
                    G=Gbase;
                    incompleteness=sum(missing_coupling_count(node_list_temp));
                    if incompleteness==0 && M>=2  %suitable for reward test
                        K=0.5;

                        for ind_node1=1:M

                            node1=node_list_temp(ind_node1);
                            mapnodes_magnetic(node1)=mapnodes_magnetic(node1)-2*G*K*(M-2);

                            for ind_node2=(ind_node1+1):M
                                node2=node_list_temp(ind_node2);
                                flag_needadd=true;
                                for cpl_pos=1:mapnodes_couplers_pos(node1)
                                    if mapnodes_couplers_with(cpl_pos,node1)==node2
                                        mapnodes_couplers_weight(cpl_pos,node1)=mapnodes_couplers_weight(cpl_pos,node1)+2*K*G;
                                        flag_needadd=false;
                                        break;
                                    end
                                end
                                if flag_needadd

                                    if mapnodes_couplers_count(node1)<max_couples && mapnodes_couplers_count(node2)<max_couples
                                        mapnodes_couplers_count(node1)=mapnodes_couplers_count(node1)+1;
                                        mapnodes_couplers_count(node2)=mapnodes_couplers_count(node2)+1;
                                        cpos=mapnodes_couplers_pos(node1)+1;
                                        mapnodes_couplers_pos(node1)=cpos;
                                        mapnodes_couplers_with(cpos,node1)=node2;
                                        mapnodes_couplers_weight(cpos,node1)=2*K*G;
                                    else
                                        missing_coupling_count(node1)=missing_coupling_count(node1)+1;
                                        missing_coupling_count(node2)=missing_coupling_count(node2)+1;
                                    end
                                    
                                end
                            end

                        end
                    else
                        %mapnodes_magnetic(node1)=max(G,mapnodes_magnetic(node1)); %prefer at least slitely a positive bias to make the spin -1, meaning still potentially a valid option, to avoid rejecting possible solutions without complete balance
                    end
                end
            end
        end


        
        

        disp(sprintf('Number of missing coupling terms: %g',sum(missing_coupling_count)));

        mapnodes_ind=mapnodes_ind(1:node_number+ancil_number); %anc
        mapnodes_loc=mapnodes_loc(1:node_number+ancil_number);%anc


        ancilla_way=false;

        %exclude using the same piece (possibly rotated) at different locations: define inequality in chain
        % For two spins the hamiltonial is (-S1+1)*(-S2+1), so there is both
        % coupling and bias
        if ancilla_way
            excd_penalty_weight=0.015;
            %G=excd_penalty_weight; was 0.02 %0.015;%0.004;%0.25;%0.001;
            donecheck=zeros(1,node_number);
            for test_in=1:node_number*(excd_penalty_weight>0)
                if  donecheck(test_in)==1
                    continue;
                end
                ind_loc_temp=mapnodes_loc(test_in);
                donecheck(test_in)=1;
                test_count=1;
                node_chain=[];
                node_chain(1)=test_in;
                for jj=test_in+1:node_number
                    if floor((mapnodes_ind(test_in)-1)/4)==floor((mapnodes_ind(jj)-1)/4) %&& test_count<4
                        test_count=test_count+1;
                        node_chain(test_count)=jj;
                        donecheck(jj)=1;
                    end
                end
                if test_count>1
                    node_list_temp=sort(node_chain);
                    M=length(node_list_temp);
                    G=excd_penalty_weight;
    
                    isodd=mod(M,2);
                    lastofA=node_list_temp(fix(M/2)+isodd);
                    K=(M+isodd)/4;
                    ancil_number=ancil_number+1; %a=ancil_number
                    mapnodes_loc(node_number+ancil_number)=ind_loc_temp;
                    mapnodes_ind(node_number+ancil_number)=-1;
                    if isodd
                        mapnodes_magnetic(node_number+ancil_number)=G*K*4;
                    end
    
                    for ind_node1=1:M
    
                        node1=node_list_temp(ind_node1);
                        mapnodes_magnetic(node1)=mapnodes_magnetic(node1)-G;
                        if isodd
                            mapnodes_magnetic(node1)=mapnodes_magnetic(node1)+4*G*K*(-1+2*(node1>lastofA));
                        end
                        if mapnodes_couplers_count(node1)<max_couples
                            mapnodes_couplers_count(node1)=mapnodes_couplers_count(node1)+1;
                            cpos=mapnodes_couplers_pos(node1)+1;
                            mapnodes_couplers_pos(node1)=cpos;
                            mapnodes_couplers_with(cpos,node1)=node_number+ancil_number;
                            mapnodes_couplers_weight(cpos,node1)=4*K*G*(-1+2*(node1>lastofA));
                        else
                            missing_coupling_count(node1)=missing_coupling_count(node1)+1;
                        end
    
                        for ind_node2=(ind_node1+1):M
                            node2=node_list_temp(ind_node2);
                            issameside=(node1>lastofA && node2>lastofA) || (node1<=lastofA && node2<=lastofA);
                            flag_needadd=true;
                            for cpl_pos=1:mapnodes_couplers_pos(node1)
                                if mapnodes_couplers_with(cpl_pos,node1)==node2
                                    mapnodes_couplers_weight(cpl_pos,node1)=mapnodes_couplers_weight(cpl_pos,node1)+2*K*G*(-1+2*(issameside));
                                    flag_needadd=false;
                                    break;
                                end
                            end
                            if flag_needadd
    
                                if mapnodes_couplers_count(node1)<max_couples && mapnodes_couplers_count(node2)<max_couples
                                    mapnodes_couplers_count(node1)=mapnodes_couplers_count(node1)+1;
                                    mapnodes_couplers_count(node2)=mapnodes_couplers_count(node2)+1;
                                    cpos=mapnodes_couplers_pos(node1)+1;
                                    mapnodes_couplers_pos(node1)=cpos;
                                    mapnodes_couplers_with(cpos,node1)=node2;
                                    mapnodes_couplers_weight(cpos,node1)=2*K*G*(-1+2*(issameside));
                                else
                                    missing_coupling_count(node1)=missing_coupling_count(node1)+1;
                                    missing_coupling_count(node2)=missing_coupling_count(node2)+1;
                                end
    
                            end
                        end
    
                    end
    
                end
            end

        else 
            %if ancilla_way: false
            excd_penalty_weight=0.25;
            donecheck=zeros(1,node_number);
            for test_in=1:node_number*(excd_penalty_weight>0)
                if  donecheck(test_in)==1
                    continue;
                end
                ind_loc_temp=mapnodes_loc(test_in);
                donecheck(test_in)=1;
                test_count=1;
                node_chain=[];
                node_chain(1)=test_in;
                for jj=test_in+1:node_number
                    if floor((mapnodes_ind(test_in)-1)/4)==floor((mapnodes_ind(jj)-1)/4) %&& test_count<4
                        test_count=test_count+1;
                        node_chain(test_count)=jj;
                        donecheck(jj)=1;
                    end
                end
                if test_count>1
                    node_list_temp=sort(node_chain);
                    M=length(node_list_temp);
                    G=excd_penalty_weight;
                    K=0.5;%(M)/4;
    
                    for ind_node1=1:M
    
                        node1=node_list_temp(ind_node1);
                        mapnodes_magnetic(node1)=mapnodes_magnetic(node1)-2*G*K*(M-2);
                        for ind_node2=(ind_node1+1):M
                            node2=node_list_temp(ind_node2);
                            flag_needadd=true;
                            for cpl_pos=1:mapnodes_couplers_pos(node1)
                                if mapnodes_couplers_with(cpl_pos,node1)==node2
                                    mapnodes_couplers_weight(cpl_pos,node1)=mapnodes_couplers_weight(cpl_pos,node1)+2*K*G;
                                    flag_needadd=false;
                                    break;
                                end
                            end
                            if flag_needadd
    
                                if mapnodes_couplers_count(node1)<max_couples && mapnodes_couplers_count(node2)<max_couples
                                    mapnodes_couplers_count(node1)=mapnodes_couplers_count(node1)+1;
                                    mapnodes_couplers_count(node2)=mapnodes_couplers_count(node2)+1;
                                    cpos=mapnodes_couplers_pos(node1)+1;
                                    mapnodes_couplers_pos(node1)=cpos;
                                    mapnodes_couplers_with(cpos,node1)=node2;
                                    mapnodes_couplers_weight(cpos,node1)=2*K*G;
                                else
                                    missing_coupling_count(node1)=missing_coupling_count(node1)+1;
                                    missing_coupling_count(node2)=missing_coupling_count(node2)+1;
                                end
    
                            end
                        end
    
                    end
    
                end
            end



        end %if ancilla_way



        disp(sprintf('Number of missing coupling terms: %g',sum(missing_coupling_count)));

        %NORMALIZE: make the parameters fit direct quantum processor
        maxJ=max(max(abs(mapnodes_couplers_weight)));
        maxH=max(max(abs(mapnodes_magnetic)));
        %excessive_factor=max(maxJ,maxH/6);%in most advanced QPU
        excessive_factor=max(maxJ,maxH/2); %usual range in QPU is 4 , but we want advantage for already set pieces
        mapnodes_couplers_weight=mapnodes_couplers_weight/excessive_factor;
        mapnodes_magnetic=mapnodes_magnetic/excessive_factor;
        %for nd=hard_node_number_start:hard_node_number_end
        %    mapnodes_magnetic(nd)=4;  %make it maximum despite the normalization
        %end
        %
        %
        %
        %
        %
        %add all the requested nodes in two lines
        fprintf(fileID,'#Definition for nodes \n');
        fprintf(fileID,sprintf('nodes = range(1, 1 + %d)\n',node_number+ancil_number));
        fprintf(fileID,sprintf('G.add_nodes_from(nodes)\n'));
        fprintf(fileID,'#Definition for edges (couplers) \n');
        for nd=1:node_number
            if mapnodes_couplers_pos(nd)>0
                csvString ='G.add_edges_from([';
                for pos=1:mapnodes_couplers_pos(nd)
                    if ~any(mapnodes_couplers_with(pos,nd)==mapnodes_couplers_with(1:pos-1,nd))
                        csvString =[csvString  sprintf('(%d,%d,{"obj":%d}),',nd,mapnodes_couplers_with(pos,nd),mapnodes_couplers_weight(pos,nd))];
                    end
                end
                fprintf(fileID,[csvString ' ])\n']);
            end
        end
        fprintf(fileID,'#Definition for magnetic field (bias) \n');
        fprintf(fileID,sprintf('magnetic_field=np.zeros(%d)\n',node_number+ancil_number));
        prev_start_nd=1;
        for nd=2:node_number+ancil_number+1
            if mapnodes_magnetic(nd)~=mapnodes_magnetic(nd-1)
                if nd>prev_start_nd+1
                    fprintf(fileID,sprintf('magnetic_field[%d:%d]=%d\n',prev_start_nd-1,nd-1,mapnodes_magnetic(nd-1)));
                else
                    fprintf(fileID,sprintf('magnetic_field[%d]=%d\n',prev_start_nd-1,mapnodes_magnetic(nd-1)));
                end
                prev_start_nd=nd;
            end
        end
        fprintf(fileID,'h = {i + 1: magnetic_field[i] for i in range(len(magnetic_field))}\n');
        fprintf(fileID,'bqm = dimod.BinaryQuadraticModel(h,{},0.0,''SPIN'')\n');
        fprintf(fileID,'for u, v, data in G.edges(data=True):\n');
        fprintf(fileID,'    bqm.add_quadratic(u, v, data[''obj''])\n');
        %fprintf(fileID,'#plot_bqm(bqm)\n');
        %fprintf(fileID,'variables = list(bqm.linear.keys())\n');
        %fprintf(fileID,'for var, bias in zip(variables, magnetic_field):\n');
        %fprintf(fileID,'    bqm.linear[var] = bias\n');

        for nd=1:node_number
            if mapnodes_notavariable(nd)==1
                fprintf(fileID,sprintf('bqm.fix_variable(%d,-1)\n',nd));
            end
        end

        %CHECK CORRECTNESS BEFORE DEPLOYING

        if deepfill==true
            checktable=zeros(bS,b4S);
            checkloc=zeros(1,bS);
            countfault=0;
            for indnode=1:node_number
                indloc=mapnodes_loc(indnode);
                checkloc(indloc)=1;
                indoption=mapnodes_ind(indnode);
                if checktable(indloc,indoption)==0
                    checktable(indloc,indoption)=1;
                else
                    countfault=countfault+1;
                end
            end
            disp(sprintf('Count of nodes: %d, board coverage= %d,  Count of faults (duplicate nodes per option): %d ',node_number,sum(checkloc)/bS,countfault));
        end



        ext_node_number=node_number+ancil_number;
        if ext_node_number<=1024
            time_hybrid=3;
        end
        if ext_node_number>1024 && ext_node_number<4096
            time_hybrid=3+(10-3)*(ext_node_number-1024)/(4096-1024)+0.5;
        end
        if ext_node_number>=4096 && ext_node_number<=10000
            time_hybrid=10+(40-10)*(ext_node_number-4096)/(10000-4096)+1;
        end
        if ext_node_number>10000 && ext_node_number<=30000
            time_hybrid=40+(200-40)*(ext_node_number-10000)/(30000-10000)+1;
        end
        if ext_node_number>30000 && ext_node_number<=100000
            time_hybrid=600;%200+(600-200)*(node_number-30000)/(100000-30000)+1;
        end
        if ext_node_number>100000 && ext_node_number<=1000000
            time_hybrid=600;
        end

        %here is altenative of using the hybrid solver, to aim at much larger graphs 
        %fprintf(fileID,'iteration = hybrid.RacingBranches(hybrid.InterruptableTabuSampler(),hybrid.EnergyImpactDecomposer(size=2) | hybrid.QPUSubproblemAutoEmbeddingSampler() | hybrid.SplatComposer()) | hybrid.ArgMin()\n');
        %fprintf(fileID,'workflow = hybrid.LoopUntilNoImprovement(iteration, convergence=3)\n');
        %fprintf(fileID,'init_state = hybrid.State.from_problem(bqm)\n');
        %fprintf(fileID,'data=workflow.run(init_state).result()\n');
        %fprintf(fileID,'    sampler_qpu = DWaveSampler(solver={''lower_noise'': False, ''qpu'': True})\n');
        if quantum_computation_paid
            fprintf(fileID,'try:\n');
            if use_hybrid_solver
                %hybrid machine:
                fprintf(fileID,'    sampler_qpu = LeapHybridSampler(solver=''hybrid_binary_quadratic_model_version2p'')\n');
            else
                %request for best available
                fprintf(fileID,'    sampler_qpu = DWaveSampler()\n');
                %newest machine:
                %fprintf(fileID,'    sampler_qpu = DWaveSampler(solver=''Advantage2_system1.4'')\n'); %most advanced, no reply
                %fprintf(fileID,'    sampler_qpu = DWaveSampler(solver=''Advantage_system4.1'')\n');
                %largest coupling size:
                %fprintf(fileID,'    sampler_qpu = DWaveSampler(solver=''Advantage2_prototype2.6'')\n');
            end
            fprintf(fileID,'    print("Connected to QPU {} .".format(sampler_qpu.solver.id))\n');
            fprintf(fileID,'except SolverNotFoundError:\n');
            fprintf(fileID,'    print("The solver is not available.")\n');
            fprintf(fileID,'    exit()\n');
        else
            if node_number<10
                fprintf(fileID,'sampler_qpu = dimod.ExactSolver()\n');
                fprintf(fileID,'print("Exact solver on the local computer, no charge")\n');
            else
                fprintf(fileID,'sampler_qpu = dimod.SimulatedAnnealingSampler()\n');
                fprintf(fileID,'print("Simulated annealing on the local computer, no charge")\n');
            end
        end
        
        if quantum_computation_paid
            if use_hybrid_solver
                %for hybrid sampler use the following
                fprintf(fileID,sprintf('result = sampler_qpu.sample(bqm,time_limit=%d)\n',time_hybrid));
                fprintf(fileID,'data=result.first.sample\n');
            else
                %for d-wave direct sampler use the following:
                fprintf(fileID,'num_reads = 100\n');
                fprintf(fileID,'result = {}\n');
                fprintf(fileID,'result[''QPU''] = EmbeddingComposite(sampler_qpu).sample(bqm, num_reads=num_reads, answer_mode=''raw'',annealing_time=20)\n');
                fprintf(fileID,'data=result[''QPU''].first.sample\n');
            end
        else
            fprintf(fileID,'result = sampler_qpu.sample(bqm)\n');%,num_reads=100
            fprintf(fileID,'print("Minimum energy= {}".format(result.first.energy))\n');
            fprintf(fileID,'data=result.first.sample\n');
        end

        fprintf(fileID,'print("Best solution found: {}".format(data))\n');
        fprintf(fileID,sprintf('for key in range(1,%d):\n',node_number+1));
        fprintf(fileID,'     data.setdefault(key, 0)\n');
        fprintf(fileID,'data_sorted=dict(sorted(data.items()))\n');
        fprintf(fileID,'array = np.array(list(data_sorted.values()))\n');
        fprintf(fileID,sprintf('filename=os.environ.get(''Ocean'')+''/QA/d_wave_answer%g.csv''\n',QAcounter));
        fprintf(fileID,'np.savetxt(filename, array,fmt=''%%d'',  delimiter='','')\n');
        fclose(fileID);

        writematrix(QAcounter,[OceanFolder '/QA/counter.txt']);
        %Here send execution sign for d-wave , click manually to continue
        disp(sprintf(' task for QM (%d,%d) ',ext_node_number,question_size));
        tempmatfilename=[OceanFolder '/QA/' sprintf('framework%g.mat',QAcounter)];
        save(tempmatfilename,'shuffled_choices_pick','options','mapnodes_loc','mapnodes_ind','-v7.3');

    else
        QAcounter=filenumber_ready;
        tempmatfilename=[OceanFolder '/QA/' sprintf('framework%g.mat',QAcounter)];
        load(tempmatfilename,'shuffled_choices_pick','options','mapnodes_loc','mapnodes_ind');
        node_number=length(mapnodes_loc(mapnodes_ind>0));
    end %if filenumber_ready==0
    
    try
        commandext=[OceanFolder '/venv/bin/python ' OceanFolder '/' sprintf('ask_d_wave%d.py',QAcounter)];
        system(commandext);
        pause(3);

        NodeSpin=readmatrix([OceanFolder '/QA/' sprintf('d_wave_answer%g.csv',QAcounter)], "Delimiter",",");
        selected_spin=NodeSpin(1:length(shuffled_choices_pick));
        %if sum(selected_spin)<0
        %    selected_spin=-selected_spin; %flip sign, so -1 is the minority
        %    NodeSpin=-NodeSpin;
        %elseif sum(selected_spin)==0 && sum(NodeSpin)<0
        %    selected_spin=-selected_spin; %flip sign, so -1 is the minority
        %    NodeSpin=-NodeSpin;
        %end
        selected_nodes=mapnodes_ind(1:length(shuffled_choices_pick));
        negative_selected_nodes=selected_nodes(selected_spin<0);
        size_result=length(negative_selected_nodes);
        if size_result==0
            fprintf(' C ')
        elseif size_result==1
            fprintf(' A ')
        else
            fprintf(' B ')
        end
        score_dwave=1*(selected_spin<0);  %1-recommended, 0-otherwise


        if deepfill==true
            supressed_count=0;
            for indnode=1:node_number
                indloc=mapnodes_loc(indnode);
                indoption=mapnodes_ind(indnode);
                if show_state(indloc)==0
                    if NodeSpin(indnode)==-1 %selected
                        if options(indloc,indoption)==0
                            %do not agree to select this state, but count as faults
                            supressed_count=supressed_count+1;
                            %options(indloc,indoption)=1; %may return states that were already eliminated by selection
                        end
                    elseif  NodeSpin(indnode)==1 %else if zero, change nothing, this spin is not reported by dwave
                        options(indloc,indoption)=0;
                    end
                end
            end
            %disp(sprintf('Number of states selected by QM but disapproved in previous runs: %d ',supressed_count));
        end
        tempmatfilename=[OceanFolder '/QA/' sprintf('framework%g.mat',QAcounter)];
        save(tempmatfilename,'shuffled_choices_pick','NodeSpin','options','mapnodes_loc','mapnodes_ind','-v7.3');
        disp(sprintf('Sum(options)=%g ',sum(options(:))));
        for t=1:256
            qv(t)=sum(options(t,:));
            if qv(t)>0
                fprintf('%g:(%g), ',t,qv(t));
            end
        end
        disp('');
        disp('QM calculation finished');

    catch
        score_dwave=ones(size(shuffled_choices_pick));
        disp(sprintf('DWAVE not available. Program prepared: Filenumber=%g',QAcounter));
    end

end


function [new_options]=relaxed_solver(threadno,options,tree,treesize)
    kseed = RandStream.create('mrg32k3a','Seed',10*threadno+floor(10*rand(1)));
    RandStream.setGlobalStream(kseed);
    new_options=zeros(size(options));
    imglin=sum(options>0,2);
    show_state=1*(imglin==1);
    for loc=1:256
        if show_state(loc)==1
            new_options(loc,:)=options(loc,:);
        end
    end
    step_dir_forward=zeros(1,256);
    step_pos_left_neighbor=zeros(1,256);
    step_dir_left=4;%zeros(1,256);
    dir_sideway=step_dir_left;
    step_indlocation=zeros(1,256);
    pos=0;
    rand_pattern=rand;
    if rand_pattern<=0.25
        dir_sideway=4;
        r=16;
        c=1;
        for t=1:8
            if c>1
                pos=pos+1;
                r=r+1;
                step_dir_forward(pos)=2;
                step_indlocation(pos)=(r-1)*16+c;
                c=c+1;
            end
            r=r+1;
            for s=1:15
               pos=pos+1;
               r=r-1;
               step_dir_forward(pos)=1;
               step_indlocation(pos)=(r-1)*16+c;
               if c>1 && r<16
                    step_pos_left_neighbor(pos)=pos-1-2*(16-r);
               end
            end
            pos=pos+1;
            r=r-1;
            step_dir_forward(pos)=2;
            step_indlocation(pos)=(r-1)*16+c;
            c=c+1;
            r=r-1;
            for s=1:15
               pos=pos+1;
               r=r+1;
               step_dir_forward(pos)=3;
               step_indlocation(pos)=(r-1)*16+c;
               if c>1 && r>1
                    step_pos_left_neighbor(pos)=pos-1-2*(r-1);
               end
            end
        end
        pos=256;
        r=r+1;
        step_indlocation(pos)=(r-1)*16+c;
        step_dir_forward(pos)=1;
        step_pos_left_neighbor(pos)=pos-1-2*(r-1);
    elseif rand_pattern>0.25 && rand_pattern<0.5
        dir_sideway=2;
        r=16;
        c=16;
        for t=1:8
            if c<16
                pos=pos+1;
                r=r+1;
                step_dir_forward(pos)=4;
                step_indlocation(pos)=(r-1)*16+c;
                c=c-1;
            end
            r=r+1;
            for s=1:15
               pos=pos+1;
               r=r-1;
               step_dir_forward(pos)=1;
               step_indlocation(pos)=(r-1)*16+c;
               if c<16 && r<16
                    step_pos_left_neighbor(pos)=pos-1-2*(16-r);
               end
            end
            pos=pos+1;
            r=r-1;
            step_dir_forward(pos)=4;
            step_indlocation(pos)=(r-1)*16+c;
            c=c-1;
            r=r-1;
            for s=1:15
               pos=pos+1;
               r=r+1;
               step_dir_forward(pos)=3;
               step_indlocation(pos)=(r-1)*16+c;
               if c<16 && r>1
                    step_pos_left_neighbor(pos)=pos-1-2*(r-1); %igonre the name, it is not left here
               end
            end
        end
        pos=256;
        r=r+1;
        step_indlocation(pos)=(r-1)*16+c;
        step_dir_forward(pos)=1;
        step_pos_left_neighbor(pos)=pos-1-2*(r-1);

    elseif rand_pattern>=0.5 && rand_pattern<=0.75
        dir_sideway=1;
        c=16;
        r=1;
        for t=1:8
            if r>1
                pos=pos+1;
                c=c+1;
                step_dir_forward(pos)=3;
                step_indlocation(pos)=(r-1)*16+c;
                r=r+1;
            end
            c=c+1;
            for s=1:15
               pos=pos+1;
               c=c-1;
               step_dir_forward(pos)=4;
               step_indlocation(pos)=(r-1)*16+c;
               if r>1 && c<16
                    step_pos_left_neighbor(pos)=pos-1-2*(16-c); %igonre the name, it is not left here
               end
            end
            pos=pos+1;
            c=c-1;
            step_dir_forward(pos)=3;
            step_indlocation(pos)=(r-1)*16+c;
            r=r+1;
            c=c-1;
            for s=1:15
               pos=pos+1;
               c=c+1;
               step_dir_forward(pos)=2;
               step_indlocation(pos)=(r-1)*16+c;
               if r>1 && c>1
                    step_pos_left_neighbor(pos)=pos-1-2*(c-1);
               end
            end
        end
        pos=256;
        c=c+1;
        step_indlocation(pos)=(r-1)*16+c;
        step_dir_forward(pos)=4;
        step_pos_left_neighbor(pos)=pos-1-2*(c-1);
    elseif rand_pattern>0.75 
        dir_sideway=3;
        c=16;
        r=16;
        for t=1:8
            if r<16
                pos=pos+1;
                c=c+1;
                step_dir_forward(pos)=1;
                step_indlocation(pos)=(r-1)*16+c;
                r=r-1;
            end
            c=c+1;
            for s=1:15
               pos=pos+1;
               c=c-1;
               step_dir_forward(pos)=4;
               step_indlocation(pos)=(r-1)*16+c;
               if r<16 && c<16
                    step_pos_left_neighbor(pos)=pos-1-2*(16-c); %igonre the name, it is not left here
               end
            end
            pos=pos+1;
            c=c-1;
            step_dir_forward(pos)=1;
            step_indlocation(pos)=(r-1)*16+c;
            r=r-1;
            c=c-1;
            for s=1:15
               pos=pos+1;
               c=c+1;
               step_dir_forward(pos)=2;
               step_indlocation(pos)=(r-1)*16+c;
               if r<16 && c>1
                    step_pos_left_neighbor(pos)=pos-1-2*(c-1);
               end
            end
        end
        pos=256;
        c=c+1;
        step_indlocation(pos)=(r-1)*16+c;
        step_dir_forward(pos)=4;
        step_pos_left_neighbor(pos)=pos-1-2*(c-1);

    end %if randpattern


    forward_connections=false(256,1024,1024);  %step number, piece index, connected piece indeX (in forward direction)
    sideway_connections=false(256,1024,1024);  %step number, piece index, connected piece indeX  (in left side direction)
    step_choices=zeros(256,700);
    step_choices_size=zeros(1,256);
    
    indp=1:256;
    indf=1:1024;
    for pos=1:255
        indloc=step_indlocation(pos);
        
        next_indloc=step_indlocation(pos+1);
        dir_forward=step_dir_forward(pos);
        vect_choice=indf(options(indloc,:)>0);
        step_choices(pos,1:length(vect_choice))=vect_choice;
        step_choices_size(pos)=length(vect_choice);
        for choice=vect_choice
            for t=1:treesize(choice,dir_forward)
                next=tree(choice,dir_forward,t);
                if next>0
                    if options(next_indloc,next)>0
                        forward_connections(pos,choice,next)=true;
                    end
                end
            end
            for t=1:treesize(choice,dir_sideway)
                next2=tree(choice,dir_sideway,t);
                if next2>0 && step_pos_left_neighbor(pos)>0
                    side_indloc=step_indlocation(step_pos_left_neighbor(pos));
                    if options(side_indloc,next2)>0
                        sideway_connections(pos,choice,next2)=true;
                    end
                end
            end
        end
    end
    racepath=zeros(256,1024,256); % node position=(step number, piece index choice), and slots that Contains chosen piece index along the path
    racepath_timer=-1*ones(256,1024);
    pos=1;
    racepath(pos,step_choices(pos),1)=step_choices(pos); %store the first position in the paths with the index number
    if step_choices_size(pos)==1
        racepath_timer(pos,step_choices(pos))=1;
    else
        racepath_timer(pos,step_choices(pos))=100;
    end
    for clk=1:256*300+40*50000
        racepath_timer=racepath_timer-1; %time works for everybody
        pos=256;
        compete_inds=indf(racepath_timer(pos,:)==0);
        if ~isempty(compete_inds)
            final_path=compete_inds(1);
            for tpos=1:256
                indloc=step_indlocation(tpos);
                ind=racepath(256,final_path,tpos);
                if ind>0
                    new_options(indloc,ind)=1;
                end
            end
            disp('Found relaxed solution');
            return;
        end
        for pos=1:255
            compete_inds=indf(racepath_timer(pos,:)==0); %this contains poistions in slot of racepath (piece indices at the current pos that their time has come to work)
            compete_size=length(compete_inds);
            if compete_size>1
                compete_selected=compete_inds(floor(rand*compete_size)+1); %piece index selected
                racepath_timer(pos,compete_inds(compete_inds~=compete_selected))=racepath_timer(pos,compete_inds(compete_inds~=compete_selected))+1; %postpone the others for next clock tick
                compete_size=1;
            elseif compete_size==1
                compete_selected=compete_inds(1);
            end
            if compete_size==1
                % initiate new message on the next position that 
                % contains the the history of the path. We deal only with one index
                % in the current position, the winning compete_selected
                next_options=step_choices(pos+1,1:step_choices_size(pos+1));
                flag_possible_vect=false(1,length(next_options));
                sideway_pos=step_pos_left_neighbor(pos+1);

                for n=1:length(next_options)
                    next_ind=next_options(n);
                    if sideway_pos>0
                        sideway_ind=racepath(pos,compete_selected,sideway_pos);
                        flag_possible_vect(n)=forward_connections(pos,compete_selected,next_ind) && sideway_connections(pos+1,next_ind,sideway_ind);
                    else
                        flag_possible_vect(n)=forward_connections(pos,compete_selected,next_ind);
                    end
                end
                badscore=50000;
                goodscore=3+min(sum(flag_possible_vect)*3,295);
                for n=1:length(next_options)
                    next_ind=next_options(n);
                    if sum(next_ind== racepath(pos,compete_selected,1:pos))>0
                        flag_repeated_use=true;
                    else
                        flag_repeated_use=false;
                    end
                    if flag_possible_vect(n) && ~flag_repeated_use
                        score=goodscore;
                    else
                        score=badscore;
                    end
                    racepath(pos+1,next_ind,1:pos)=racepath(pos,compete_selected,1:pos);
                    racepath(pos+1,next_ind,pos+1)=next_ind;
                    racepath_timer(pos+1,next_ind)=score; %set delay for action
                end
            end
        end

    end % for clk
    maxpos=0;
    for tpos=1:256
        if sum(racepath_timer(tpos,:)>=0)>0
            maxpos=tpos;
        end
    end
    compete_inds=indf(racepath_timer(maxpos,:)>=0);
    final_path=compete_inds(1);
    for tpos=1:maxpos
        indloc=step_indlocation(tpos);
        ind=racepath(maxpos,final_path,tpos);
        if ind>0
            new_options(indloc,ind)=1;
        end
    end
    disp('Found relaxed solution, but not for entire board');



end