function miniSolver10()
global screenCoord;
global mousepress;
global CH;
global e2folder;
global flagloop;
global extraline;
global bW; %board width
global bS;%size of board
global b4S;
global quantum_computation_paid;
global use_only_QPU;

bW=input('Board dimension: '); %board width
bS=bW*bW;%size of board
b4S=bS*4;
payquest=input('Use DWAVE paid subscription (1) or simulation (0): ');
quantum_computation_paid=logical(payquest);
if payquest==1
    payquestQPU=input('Fit to quantum processor unit (1) or use hybrid solver (0): ');
    use_only_QPU=logical(payquestQPU);
end


quantum_d_wave=false; %this is the slash option to use dwave while  h function

compromize=false;
extraline=0;
flagloop=false;
e2folder=input('Eternity folder: ','s');
screenCoord=[0 0];
%parallelForce=input('Number of CPU cores? ');%19;
disp('Mini puzzle solver written by Shahar Seifer (C) 2025');
disp('Directions:  Hover with mouse over tiles and press key. Do not move the board window.');
disp('w- deduction, h- nucleation, f- fix faults and show condition ');
disp('c- remove, v- match tiles, up/down arrows- select different match, f- fix one selection per shown tile, g-propagate one step');
disp('t- enter tile as 4-letter text, z- show a list of all open choices, space- fit all general choices, l- load board state file');
disp('s- save board state,  e- export txt file compatible with eternityII editor, esc- end program.');
disp('/ - toggle for light assistance from dwave during nucleation (frame block). ');
disp('1- run full caclulation in D-Wave');
disp(' ');

original=[];

askifnew=input('Build a new game? 1-yes: ');
if askifnew==1
   for r=1:bW
       for c=1:bW
          vect=-1*ones(1,4);
          if r==1
              vect(1)=0;
          elseif r==bW
              vect(3)=0;
          end
          if c==1
              vect(4)=0;
          elseif c==bW
              vect(2)=0;
          end
          if c>1
              ob=original{r,c-1};
              vect(4)=ob(2);
          end
          if r>1
              ob=original{r-1,c};
              vect(1)=ob(3);
          end
          for t=1:4
              if vect(t)==-1
                 vect(t)=floor(8*rand(1))+1;
              end
          end
          original{r,c}=vect;
       end
   end
   filename=0;
else
    % load(sprintf('%s\\origin_mini_virtual.mat',e2folder),'-mat'); %load originial cell array
    disp('Original table can be loaded with the mat file of the board state.');
    [filename,path] = uigetfile(sprintf('%s\\*.mat',e2folder),'Fetch mat of saved board');
end



planv=0;

keep_options=[];
machine_state=0;


isprob=0;


if sum(filename>0)>0
    filen=[path filename];
    disp(['loading ' filen]);
    load(filen,'-mat');
    keep_options=options;
    if isempty(isprob)
        isprob=0;
    end
else
    if planv==0
        flaggen=input('0- empty board,  1- generate full board,  2- generate frame ? ');
        if flaggen==1
            options=zeros(bS,b4S);
            for ind_location=1:bS
                options(ind_location,4*(ind_location-1)+1)=1;
            end
            keep_options=options;
        elseif flaggen==0
            keep_options=zeros(bS,b4S);
        elseif flaggen==2
            options=zeros(bS,b4S);
            for ind_location=1:bW
                options(ind_location,4*(ind_location-1)+1)=1;
            end
            for ind_location=bS-(bW-1):bS
                options(ind_location,4*(ind_location-1)+1)=1;
            end
            for ind_location=bW+1:bW:bS-(2*bW)+1
                options(ind_location,4*(ind_location-1)+1)=1;
            end
            for ind_location=(2*bW):bW:bS-bW
                options(ind_location,4*(ind_location-1)+1)=1;
            end
            keep_options=options;
        end
    else
        keep_options=zeros(bS,b4S);
    end
end



%load oringal{r:1-bW,c=1-bW}= [sym_up, sym_right, sym_down, sym_left]
piece=zeros(bW*bW,4);
for row=1:bW
    for col=1:bW
        ob=original{row,col};
        ind=(col-1)+(row-1)*bW+1;
        piece(ind,1:4)=ob;
    end
end




%build tree:  from obj1(piece1,rot1) via leg (1:up, 2:right, 3:down,4:left) to obj2(piece2,rot2)
% if edge symbol write to zero object
% The rotations are all anticlockwise !
%Initial state, all connections included. [piece2]=floor((obj2-1)/4)+1, [rot2]=(obj2-1)%4
tree=-1*ones(bS*4,4,57);  % from which object according to position (1 to bS*4) to all possible objects(1 of bS*4) up to 128 (actually there are 49 max),  according to conntection direction (1 of 4).
treesize=zeros(bS*4,4);
for direct1=1:4
    for piece1=1:bS
        for rot1=1:4   %rotation(1,2,3,4)=(no rot, 1 CW, 2 CW, 3 CW),  direction(1,2,3,4)=up, right,down, left
            
            treeind1=(piece1-1)*4+rot1;
            sym1=piece(piece1,mod(direct1-1-(rot1-1),4)+1);
            if sym1==0
                tree(treeind1,direct1,1)=0;
                treesize(treeind1,direct1)=1;
                continue;
            end
            for piece2=1:bS
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



options=double(ones(bS,b4S)); %location=(row-1)*bW+col,  index valid table


tempv=1:bS;


for row=2:(bW-1)
    for col=2:(bW-1)
        options(location(row,col),1:bW*4)=0;
        options(location(row,col),(bS-bW)*4+1:bS*4)=0;
        for row2=2:(bW-1)  
            options(location(row,col),(location(row2,1)-1)*4+1:(location(row2,1)-1)*4+4)=0;
            options(location(row,col),(location(row2,bW)-1)*4+1:(location(row2,bW)-1)*4+4)=0;
        end
    end
    
end
%define options to margin locations
for r=1:bW
    options(location(r,1),:)=0;
    options(location(r,bW),:)=0;
end
for c=2:(bW-1)
    options(location(1,c),:)=0;
    options(location(bW,c),:)=0;
end
options(:,1:4:(bW-1)*4+1)=0;
options(:,4*bW*(bW-1)+3: 4 :4*bW*(bW-1)+(bW-1)*4+3)=0;
options(:,2: 4*bW :4*bW*(bW-1)+2)=0;
options(:,4*(bW-1)+4: 4*bW :4*bW*(bW-1)+4*(bW-1)+4)=0;

%corners may be exchanged after rotation
options(location(1,1),1)=1; %use left up
options(location(1,1),4*(bW-1)+4)=1; %use right up
options(location(1,1),4*bW*(bW-1)+2)=1; %use left down
options(location(1,1),4*bW*(bW-1)+4*(bW-1)+3)=1; %use right down
options(location(1,bW),2)=1; %use left up
options(location(1,bW),4*(bW-1)+1)=1; %use right up
options(location(1,bW),4*bW*(bW-1)+3)=1; %use left down
options(location(1,bW),4*bW*(bW-1)+4*(bW-1)+4)=1; %use right down
options(location(bW,1),4)=1; %use left up
options(location(bW,1),4*(bW-1)+3)=1; %use right up
options(location(bW,1),4*bW*(bW-1)+1)=1; %use left down
options(location(bW,1),4*bW*(bW-1)+4*(bW-1)+2)=1; %use right down
options(location(bW,bW),3)=1; %use left up
options(location(bW,bW),4*(bW-1)+2)=1; %use right up
options(location(bW,bW),4*bW*(bW-1)+4)=1; %use left down
options(location(bW,bW),4*bW*(bW-1)+4*(bW-1)+1)=1; %use right down
%margin-non-corner pieces may be excahnge after rotation
for c=2:(bW-1)
    options(location(1,c),4+1:4:(bW-2)*4+1)=1;
    options(location(1,c),4*bW*(bW-1)+4+3: 4 :4*bW*(bW-1)+(bW-2)*4+3)=1;
    options(location(1,c),4*bW+2: 4*bW :4*bW*(bW-2)+2)=1;
    options(location(1,c),4*bW+4*(bW-1)+4: 4*bW :4*bW*(bW-2)+4*(bW-1)+4)=1;
    
    options(location(bW,c),4*bW*(bW-1)+4+1:4:4*bW*(bW-1)+(bW-2)*4+1)=1;
    options(location(bW,c),4+3:4:(bW-2)*4+3)=1;
    options(location(bW,c),4*bW+4: 4*bW :4*bW*(bW-2)+4)=1;
    options(location(bW,c),4*bW+4*(bW-1)+2: 4*bW :4*bW*(bW-2)+4*(bW-1)+2)=1;
end
for r=2:(bW-1)
    options(location(r,1),4+4:4:(bW-2)*4+4)=1;
    options(location(r,1),4*bW*(bW-1)+4+2: 4 :4*bW*(bW-1)+(bW-2)*4+2)=1;
    options(location(r,1),4*bW+1: 4*bW :4*bW*(bW-2)+1)=1;
    options(location(r,1),4*bW+4*(bW-1)+3: 4*bW :4*bW*(bW-2)+4*(bW-1)+3)=1;
    
    options(location(r,bW),4*bW*(bW-1)+4+4:4:4*bW*(bW-1)+(bW-2)*4+4)=1;
    options(location(r,bW),4+2:4:(bW-2)*4+2)=1;
    options(location(r,bW),4*bW+3: 4*bW :4*bW*(bW-2)+3)=1;
    options(location(r,bW),4*bW+4*(bW-1)+1: 4*bW :4*bW*(bW-2)+4*(bW-1)+1)=1;
end

baseline_options=options;
baseline_tree=tree;

if ~isempty(keep_options)
    options=keep_options;
end

required_options=zeros(bS,b4S);
for t=1:bS
    required_options(t,(t-1)*4+1)=1;
end
RRR=lsqminnorm(double(options),ones(bS,1));
chosenRRR=RRR;
vec4=1:4;
for ind=1:4:b4S
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



%options=options.*(chosenRRR>0)';



for ind_location=1:bS
    if sum(options(ind_location,:)>0)==1
        indf=1:b4S;
        indu=1:bS;
        indchoice=indf(options(ind_location,:)>0);
        indstart=floor((indchoice-1)/4)*4+1;
        options(indu~=ind_location,indstart:indstart+3)=0;
    end
end



if planv==0
    actual_sol=sum(sum((options>0.5).*(required_options>0.5),2))

end
sol=sum(sum(options>0.5,2))

indf=1:b4S;
imglin=sum(options>0,2);
show_state=1*(imglin==1);
show_index=zeros(size(show_state));
for ind_location=1:bS
    if show_state(ind_location)==1
       show_index(ind_location)=min(indf(options(ind_location,:)>0));
    end
end
[options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);

keep_options=options;
keep_show_index=show_index;
keep_show_state=show_state;

dwave_options=[];


parallelstate=false;
delete(gcp('nocreate'));

symbollabel=['@' 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];

%RRR=lsqminnorm(double(options),ones(bS,1));
%agebraic_solved=sum(RRR>0.5)

%h.HTMLSource=htmltxt;
fig = uifigure;
fig.WindowState = 'maximized';
fig.KeyPressFcn = @onKeyPress;           % prefer this in uifigure
fig.WindowKeyPressFcn  = @onKeyPress;
fig.Interruptible = 'on';
fig.BusyAction = 'cancel';
fig.Position=[38   198   787   818];

fig.WindowButtonDownFcn = @(src,evt) onMouseDownRefresh(src, evt);

% Create uiaxes
ax = uiaxes(fig);
ax.PickableParts = 'visible';
ax.HitTest = 'on';
h = uihtml(fig);
%h.Position = [ 40 40 fig.Position(3) fig.Position(4)]; 



flagloop=false;

for testno=1:10000


    %SHOWING SESSION
    neworiginal=fill_neworiginal(show_state,show_index,options,piece);
    show_original(h,fig,neworiginal,testno);

    mousepress=0;
    CH='';
    screenCoord=[0 0];
    if flagloop==false && machine_state==0
        while (sum(screenCoord)==0 && mousepress==0 && isempty(CH))
            pause(0.5);
        end
        %disp(screenCoord);
    
        col=floor(16*(screenCoord(1)-91)/(1160-491))+1;
        row=floor(16*(989-screenCoord(2))/(989-270))+1;
        if col>bW
            col=bW;
        end
        if row>bW
            row=bW;
        end
        if col<1 || row<1
            col=1;
            row=1;
        end
    else
        col=1;
        row=1;
        CH='';
    end

    h_number=-1;

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
        elseif machine_state==1
            show_index=keep_show_index;
            show_state=keep_show_state;

            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            machine_state=machine_state+1;
        end
    end

    if strcmp(CH,'escape')
        close(fig);
        return;
    end


    if col>=1 && col<=bW && row>=1 && row<=bW
        ind_location=(col-1)+(row-1)*bW+1;
        if strcmp(CH,'uparrow') && show_state(ind_location)>0
            for t=show_index(ind_location)+1:b4S
                if options(ind_location,t)>0
                    show_index(ind_location)=t;
                    break;
                end
            end
        end
        if strcmp(CH,'downarrow') && show_state(ind_location)>0
            for t=show_index(ind_location)-1:-1:1
                if options(ind_location,t)>0
                    show_index(ind_location)=t;
                    break;
                end
            end
        end
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
            quantum_d_wave=~quantum_d_wave;
            disp(sprintf('Using d-wave=  %d',quantum_d_wave))
            pause(1);
        end
        if strcmp(CH,'return') && sum(options(ind_location,:)>0)>0
            show_state(ind_location)=1; %open recommeded option
            options(ind_location,:)=keep_options(ind_location,:);
            show_index(ind_location)=min(indf(options(ind_location,:)>0));
        end
        if strcmp(CH,'s') 
            %options=force_to_options(options,baseline_options);
            filesave=sprintf('%s\\savestate%d_8x8.mat',e2folder,testno);
            save(filesave,'options','original','-mat');
            disp(sprintf('saving savestate%d',testno));
            pause(1);
        end
        if strcmp(CH,'hyphen') || strcmp(CH,'1') || strcmp(CH,'2') || strcmp(CH,'3') || strcmp(CH,'4')  || strcmp(CH,'5')  || strcmp(CH,'6') || strcmp(CH,'7') || strcmp(CH,'8') || strcmp(CH,'9')
            %options=force_to_options(options,baseline_options);
            disp('Effecting options table of entire board..');
            deepfill=true;
            filenumber_ready=0;
            
            problemsize=input('Enter number of spins to invest in the problem: ')
            use_ancilla=input('Choose  1-soft constraint with ancilla method,  0-hard constraint: ')
            ancilla_way=logical(use_ancilla);
        
            if strcmp(CH,'hyphen')
                filenumber_ready=input('Enter QA file number to retrieve:  ');
                problemsize=100000; %dummy
            end
            shuffled_choices_pick=indf(options(ind_location,:)>0);
            [~, dwave_options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,filenumber_ready,ancilla_way);
            if ~quantum_d_wave %if not during slash mode
                options=dwave_options;
                imglin=sum(options>0,2);
                show_state=1*(imglin==1);
                show_index=zeros(size(show_state));
                for indlocation=1:bS
                    if show_state(indlocation)==1
                       show_index(indlocation)=min(indf(options(indlocation,:)>0));
                    end
                end
                tree=update_tree(tree,show_state,show_index);
                disp('options table updated by dwave computation.Use f to convert empty positions to base options.');
            else %slash mode, so use dwave to prepare for h function
                disp('dwave_options updated by dwave calculation, which can be used in h function.')
                filesave=sprintf('%s\\dwave_latest.mat',e2folder);
                save(filesave,'dwave_options','-mat');
                disp('dwave_options saved.');
                quantum_d_wave=false;
                dwimglin=sum(dwave_options>0,2);
                dwshow_state=1*(dwimglin==1);
                dwshow_index=zeros(size(dwshow_state));
                for dwindlocation=1:bS
                    if dwshow_state(dwindlocation)==1
                       dwshow_index(dwindlocation)=min(indf(dwave_options(dwindlocation,:)>0));
                    end
                end
                [dwave_options,dwshow_state,dwshow_index]=force_to_options(dwave_options,options,dwshow_state,dwshow_index);
                disp('dwave_options regularized');        
                disp('slash function deactivated:  quantum_d_wave=false.');
            end
            pause(1);
        end
        if strcmp(CH,'0')
            disp('Effecting only the one piece, and check for consistency');
            deepfill=true;
            bakup_options=options;
            row_in_options0=options(ind_location,:);
            shuffled_choices_pick=indf(options(ind_location,:)>0);
            problemsize=input('How many qubits to invest in the start run? ');
            problemsize0=problemsize;
            [~, options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,0,false);
            row_in_options1=options(ind_location,:);
            options=bakup_options;
            if sum(row_in_options1==row_in_options0)<b4S
                problemsize=problemsize0*2;
                [~, options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,0,false);
                row_in_options2=options(ind_location,:);
                options=bakup_options;
                if sum(row_in_options1==row_in_options2)==b4S
                    options(ind_location,:)=row_in_options2; %accept verified update in selections of the pointed piece
                    disp('Options table has been updated in the location of the piece according to verified dwave computation.');
                else
                    disp('No update')
                    do_more=input('Check with more qutbits? (0-no, 1-yes): ');
                    if do_more==1
                        problemsize=problemsize0*4;
                        [~, options]=dwave_solve(ind_location,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,0,false);
                        row_in_options3=options(ind_location,:);
                        options=bakup_options;
                        if sum(row_in_options2==row_in_options3)==b4S
                            options(ind_location,:)=row_in_options3; %accept verified update in selections of the pointed piece
                            disp('Options table has been updated in the location of the piece according to verified dwave computation.');
                        else
                            disp('No update also in further trial')
                        end
                    end
                end
            else %NOT if sum(row_in_options1==row_in_options0)<b4S, meaning that QM did not eliminate any option
                disp('Trial shows no benefit from QM.')
            end

            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:bS
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                end
            end
            tree=update_tree(tree,show_state,show_index);
            pause(1);
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

            for ind=1:b4S
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
                        if all(next_row==0 || next_row==(bW+1) || next_col==0 || next_col==(bW+1))
                            next_sym=0;
                        else
                            next_location=(next_col-1)+(next_row-1)*bW+1;
                            indf=1:b4S;
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
            indu=1:bS;
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
                for testchoice=1:maxchoices
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
                        if abs(score-bS)<abs(best_score-bS)
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
            [filename,path] = uigetfile(sprintf('%s\\*.mat',e2folder),'Fetch mat file with options');
            loadsave=[path filename];
            disp(['loading ' loadsave]);
            if ~quantum_d_wave
                load(loadsave,'-mat');
                keep_options=options;
                imglin=sum(options>0,2);
                show_state=1*(imglin==1);
                show_index=zeros(size(show_state));
                for indlocation=1:bS
                    if show_state(indlocation)==1
                       show_index(indlocation)=min(indf(options(indlocation,:)>0));
                    end
                end
            elseif quantum_d_wave
                remoptions=options;
                remshow_state=show_state;
                remshow_index=show_index;
                load(loadsave,'-mat');
                imglin=sum(options>0,2);
                show_state=1*(imglin==1);
                show_index=zeros(size(show_state));
                for indlocation=1:bS
                    if show_state(indlocation)==1
                       show_index(indlocation)=min(indf(options(indlocation,:)>0));
                    end
                end
                [options,~,~]=force_to_options(options,remoptions,show_state,show_index);
                dwave_options=options;
                options=remoptions;
                show_state=remshow_state;
                show_index=remshow_index;
                disp('File loaded to dwave_options because this is a slash mode. Use h function next.')
                quantum_d_wave=false;
                disp('Slash mode deactivated.');
            end

        end
        if strcmp(CH,'comma')

            [options,show_state,show_index]=choose_CornersNN(options,show_state,show_index,tree,treesize);

            pause(1);
        end

        if strcmp(CH,'d')
            [filename,path] = uigetfile(sprintf('%s\\*.mat',e2folder),'Fetch mat file with options');
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
            for indlocation=1:bS
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
            indf=1:b4S;
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
        if strcmp(CH,'f')

            if planv==0
                required_options=zeros(bS,b4S);
                for t=1:bS
                    required_options(t,(t-1)*4+1)=1;
                end
                firstcheck_actual_sol=sum(sum((options>0.5).*(required_options>0.5),2));
                disp(sprintf('Still could be solved= %d  ',firstcheck_actual_sol));
            end
            [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
            disp('operation: f');
            if planv==0
                required_options=zeros(bS,b4S);
                for t=1:bS
                    required_options(t,(t-1)*4+1)=1;
                end
                firstcheck_actual_sol=sum(sum((options>0.5).*(required_options>0.5),2));
                disp(sprintf('Still could be solved= %d  ',firstcheck_actual_sol));
            end
            tree=update_tree(baseline_tree,show_state,show_index);

        end
        if strcmp(CH,'g')
            [options]=correct_options(options,tree);
            imglin=sum(options>0,2);
            show_state=1*(imglin==1);
            show_index=zeros(size(show_state));
            for indlocation=1:bS
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                else
                    show_index(indlocation)=0;
                end
            end

            tree=update_tree(tree,show_state,show_index);

        end

         if strcmp(CH,'h') || flagloop || machine_state>=4 
            disp(sprintf('H=%d,  repeat=%d',h_number,flagloop));
            if ~flagloop
                if isprob==0
                    [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index);
                end
                loop_options=options;
                loop_show_state=show_state;
                loop_show_index=show_index;
            else
                options=loop_options;
                show_state=loop_show_state;
                show_index=loop_show_index;
            end
           
            [options,show_state,show_index,failblock]=form_block(options,baseline_options,show_state,show_index,piece,baseline_tree,e2folder,compromize,isprob,h_number,quantum_d_wave,dwave_options);
            if sum(show_state>0)>=230 || failblock
                flagloop=false;
            end
            disp(sprintf('Solved= %d',sum(show_state>0)));
            if sum(show_state>0)>150
                filesave=sprintf('%s\\savesolve%d_C%s.mat',e2folder,testno,clue_string);
                save(filesave,'options','-mat');
                disp(sprintf('saving savestate%d',testno));
                pause(1);
            end

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
            for indlocation=1:bS
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
            RRR=lsqminnorm(double(options),ones(bS,1));
            extra=max(sum(RRR),0);
            algb_score=abs(bS-extra)+1;
            rank_score=rank(double(options));
            disp(sprintf('algebraic score=%d,  Rank=%d',algb_score,rank_score));
            if planv==0
                required_options=zeros(bS,b4S);
                for t=1:bS
                    required_options(t,(t-1)*4+1)=1;
                end
                firstcheck_actual_sol=sum(sum((options>0.5).*(required_options>0.5),2));
                disp(sprintf('Still could be solved= %d  ',firstcheck_actual_sol));
            end

        end
    end
    screenCoord=[0 0];
    pause(1);
end

end %main function


function [options,show_state,show_index]=force_to_options(options,baseline_options,show_state,show_index)
    global bW; %board width
    global bS;%size of board
    global b4S;
            %force the displayed pieces chosen by user to be one option
            %only. Meaning they are now clues.
            indf=1:b4S;
            for indlocation=1:bS
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
            for indlocation=1:bS
                if show_state(indlocation)==1
                   show_index(indlocation)=min(indf(options(indlocation,:)>0));
                else
                    show_index(indlocation)=0;
                end
            end

            for indlocation=1:bS
                if show_state(indlocation)>0 
                    indf=1:b4S;
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
                for indlocation=1:bS
                    if sum(options(indlocation,:)>0,2)==0 %needed after dwave to fill empty lines
                        options(indlocation,:)=baseline_options(indlocation,:);
                        show_state(indlocation)=0;
                    end
                end
                for indlocation=1:bS
                    if sum(options(indlocation,:)>0)==1 && show_state(indlocation)>0
                        indf=1:b4S;
                        indu=1:bS;
                        indchoice=indf(options(indlocation,:)>0);
                        show_index(indlocation)=indchoice;
                        indstart=floor((indchoice-1)/4)*4+1;
                        options(indu~=indlocation,indstart:indstart+3)=0;
                    end
                end
            end

end

function [options,show_state,show_index]=force_to_home(options,baseline_options,show_state,show_index)
    %necessary for h function, after force_to_options
    global bW; %board width
    global bS;%size of board
    global b4S;
    %for indlocation=1:bS
    %    if show_state(indlocation)==0
    %        options(indlocation,:)=baseline_options(indlocation,:);
    %    end
    %end

    for rpt=1:5
        for indlocation=1:bS
            if sum(options(indlocation,:)>0,2)~=1
                options(indlocation,:)=baseline_options(indlocation,:);
                show_state(indlocation)=0;
            end
        end
        for indlocation=1:bS
            if sum(options(indlocation,:)>0)==1 && show_state(indlocation)>0
                indf=1:b4S;
                indu=1:bS;
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




function locationr=location(r,c)
    global bW; %board width
    global bS;%size of board
    global b4S;
    locationr=bW*(r-1)+c;
end



function neworiginal=fill_neworiginal(show_state,show_index,options,piece)
    global bW; %board width
    global bS;%size of board
    global b4S;

    indf=1:b4S;
    for row=1:bW
        for col=1:bW
            ind_location=(col-1)+(row-1)*bW+1;
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
                    neworiginal{row,col}=[-1*(sum(options(ind_location,:)>0)==0) 0 0 0];
                end
            elseif sum(options(ind_location,:)>0)>0
                neworiginal{row,col}=[0 0 0 0];
            else
                neworiginal{row,col}=[-1 0 0 0];
            end
        end
    end

end




function [new_optionsr,score_nostop]=new_options(compromize,options,tree)
global bW; %board width
global bS;%size of board
global b4S;

    new_optionsr=options>0;
    score_nostop=1;
    imglin=sum(new_optionsr,2);
    score_solved_prev=sum(imglin==1);

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

    for clk=2:255   
       %%prev_sumoptions=sum(new_optionsr(:));
       prev_options=new_optionsr;
       vect_ind=1:bS*4;
       %prev_options=options;
       for location_ind=1:bS
            %if sum(new_optionsr(location_ind,1:4:end) | new_optionsr(location_ind,2:4:end) | new_optionsr(location_ind,3:4:end) | new_optionsr(location_ind,4:4:end) )==1 %certinty in one location
            %    ind=max(vect_ind(new_optionsr(location_ind,:)));
            %    remain_loc=[1:location_ind-1 location_ind+1:bS];
            %    ind_p=floor((ind-1)/4)*4+1:floor((ind-1)/4)*4+4;
            %    new_optionsr(remain_loc,ind_p)=0; %remove inds for the piece already determined from all the other locations 
            %end
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
                
                %$if next_location==43
                 %   dummy=0;
                %end

                ind_flags=false(1,bS*4);
                v_t=vect_ind(prev_options(location_ind,:));
                if ~isempty(v_t) %|| ~compromize
                    for jj=1:length(tree(1,1,:))
                        indnext=tree(v_t,direction,jj)';
                        indnextn0=indnext(indnext>0);
                        tempflag=false(1,bS*4);
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
           if score_nostop>0 || compromize
               score_nostop=score_solved;
           end
           break;
       end
       if score_solved>score_solved_prev  %update tree if number of definite solved pieces increased
            score_solved_prev=score_solved;
            show_state=1*(imglin==1);
            for indlocation=1:bS
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
    global bW; %board width
    global bS;%size of board
    global b4S;

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
global bW; %board width
global bS;%size of board
global b4S;
       new_tree=tree;
       slen=length(tree(1,1,:));
       vect_ind=1:bS*4;
       for location_ind=1:bS
            if show_state(location_ind)>0 && show_index(location_ind)>0
                ind=show_index(location_ind);
                indpiecelist=floor((ind-1)/4)*4+1:floor((ind-1)/4)*4+4;
                indremovelist=indpiecelist(indpiecelist~=ind);
                new_tree(indremovelist,:,:)=-1;
            end
       end

       for location_ind=1:bS
            %if sum(new_optionsr(location_ind,1:4:end) | new_optionsr(location_ind,2:4:end) | new_optionsr(location_ind,3:4:end) | new_optionsr(location_ind,4:4:end) )==1 %certinty in one location
            %    ind=max(vect_ind(new_optionsr(location_ind,:)));
            %    remain_loc=[1:location_ind-1 location_ind+1:bS];
            %    ind_p=floor((ind-1)/4)*4+1:floor((ind-1)/4)*4+4;
            %    new_optionsr(remain_loc,ind_p)=0; %remove inds for the piece already determined from all the other locations 
            %end
            if all(show_state(location_ind)==0 ||show_index(location_ind)==0)
                continue;
            end
            here_index=show_index(location_ind);
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


function show_original(h,fig,original,testno)
global bW; %board width
global bS;%size of board
global b4S;
	global e2folder;
    warning('off','all');
    meshr=bW;
    meshc=bW;
    symbolname{1}=-1;
    current=original;
      rowname={'a' ' ' ' ' 'b' ' ' ' ' 'c' ' ' ' ' 'd' ' ' ' ' 'e'  ' ' ' ' 'f' ' ' ' '  'g'  ' ' ' ' 'h' ' ' ' '  'i' ' ' ' '  'j' ' ' ' '  'k' ' ' ' '  'l' ' ' ' '  'm' ' ' ' '  'n' ' ' ' '  'o' ' ' ' '  'p' ' ' ' ' };
        colname={' ' '1' ' ' ' '  '2' ' ' ' '  '3' ' ' ' '  '4' ' ' ' '  '5' ' ' ' '  '6' ' ' ' '  '7' ' ' ' '  '8' ' ' ' '  '9' ' ' ' '  '10' ' ' ' '  '11' ' ' ' ' '12'  ' ' ' '   '13'  ' ' ' '   '14'  ' ' ' '   '15'   ' ' ' '   '16' ' ' ' ' ' '  };
        if symbolname{1}==-1
            symbolname={'@' 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};
        end
        if mod(meshc,3)>0
            colorlist={'yellow' 'aqua' '#F2DBDB'};
        elseif mod(meshc,4)>0
            colorlist={'yellow' 'aqua' '#F2DBDB' '#DAEEF3'};
        else
            colorlist={'yellow' 'aqua' '#F2DBDB' '#DAEEF3' '#EAF1DD'};
        end        
        readoutflname=sprintf('%s\\readout%d.html',e2folder,testno);
        %readoutflnamemax='C:\\Users\\seifer\\Documents\\Matlab\\eternity\\readoutmax.html';
        %readoutflnamew=readoutflname;
        try
            delete (readoutflname);
        catch
        end
        if testno>1
            try
                readoutflname_prev=sprintf('%s\\readout%d.html',e2folder,testno-1);
                delete (readoutflname_prev);
            catch
            end
        end
        fid2=fopen(readoutflname,'w');
        fprintf(fid2,['<html>' char(10)]);
        fprintf(fid2,['<head>' char(10)]);
        fprintf(fid2,['<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/></head>' char(10)]);
        fprintf(fid2,['<body style="tab-interval:8.0pt">' char(10)]);
        fprintf(fid2,['<div class="Section1">' char(10)]);
        fprintf(fid2,['<table align="left">' char(10)]);
        fprintf(fid2,['<tr>' char(10)]);
        for cc=1:meshc*3+2
            fprintf(fid2,['<td bgcolor="white"><p ><span style="font-size:8.0pt">' colname{cc} '</span></p></td>' char(10)]);
        end
        for r=1:meshr
            for subline=1:3
                fprintf(fid2,['<tr>' char(10)]);
                fprintf(fid2,['<td bgcolor="white"><p ><span style="font-size:8.0pt">' rowname{(r-1)*3+subline} '</span></p></td>' char(10)]);
                for c=1:meshc
                  colorsel=mod((r-1)*meshc+(c-1),length(colorlist))+1;
                  flag_empty=false;
                  if current{r,c}(1)==-1
                      current{r,c}(1)=0;
                      flag_empty=true;
                  end
                
                  if subline==1 || subline==3
                      fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt"></span></p></td>' char(10)]);
                      if flag_empty
                          fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt;color:'  'white'   '">'  '@'  '</span></p></td>' char(10)]);
                      else
                          fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt;color:'  'black'   '">'  symbolname{current{r,c}(subline)+1}  '</span></p></td>' char(10)]);
                      end
                      fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt"></span></p></td>' char(10)]);
                  elseif subline==2
                      if flag_empty
                          fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt;color:'  'white'   '">'  '@'  '</span></p></td>' char(10)]);
                      else
                          fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt;color:'  'black'   '">'  symbolname{current{r,c}(4)+1}  '</span></p></td>' char(10)]);
                      end
                      fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt"></span></p></td>' char(10)]);
                      if flag_empty
                          fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt;color:'  'white'   '">'  '@'  '</span></p></td>' char(10)]);
                      else
                          fprintf(fid2,['<td bgcolor="' colorlist{colorsel}  '"><p ><span style="font-size:8.0pt;color:'  'black'   '">'  symbolname{current{r,c}(2)+1}  '</span></p></td>' char(10)]);
                      end
                  end
                  if flag_empty
                      current{r,c}(1)=-1;
                  end
                end %for c
                fprintf(fid2,['<td bgcolor="white"><p ><span style="font-size:8.0pt">' rowname{(r-1)*3+subline} '</span></p></td>' char(10)]);
                fprintf(fid2,['</tr>' char(10)]);
            end %for subline=1:3
        end %for r
        fprintf(fid2,['<tr>' char(10)]);
        for cc=1:meshc*3+2
            fprintf(fid2,['<td bgcolor="white"><p ><span style="font-size:8.0pt">' colname{cc} '</span></p></td>' char(10)]);
        end
        fprintf(fid2,['</tr>' char(10)]);
    
        fprintf(fid2,['</table>' char(10)]);
        fprintf(fid2,['</div>' char(10)]);
        fprintf(fid2,['</body>' char(10)]);
        fprintf(fid2,['</html>' char(10)]);
        fclose(fid2);
        %[stat, h]=web (readoutflnamew,'-browser');
        %htmltxt=fileread(readoutflname);
        h.HTMLSource = readoutflname;
        h.Position = [ 40    40   738   769]; 
        focus(fig);
        %h.ButtonDownFcn = @onClickImage;
        %gcoor=ginput(1);
        %close(fig);


end


% Capture clicks on uiaxes
% Capture key presses
function onKeyPress(~, event)
    global screenCoord;
    global CH;
    % Which kwy was pressed?
    clickType = event.Key;
    screenCoord = get(0,'PointerLocation');
    CH=clickType;
    %disp(CH);
    
    % remove datatip
    if strcmp(clickType,'escape')
        try
            dt = findall(fig,'Type','datatip');
            if ~isempty(dt)
                delete(dt)
            end
        catch
            %
        end
    end

end

function onMouseDownRefresh(fig, ~)
    global mousepress;
    mousepress=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%    ANALYZE    %%%%%%%%%%%
function [score1,options]=analyze(level,tree,treesize,options)
global bW; %board width
global bS;%size of board
global b4S;
   compromize=false;
   indf=1:b4S;
   last_options=options;
   backup_tree=tree;

   for level=0:level     
    

        %Perturbations
        %perturb_vectors
        perturb_vectors=zeros(bS*4,1000,bS);
        perturb_param=zeros(bS*4,1000,5);
        perturb_score=zeros(bS*4,1000);
        prev_min_score=zeros(bS*4,1);
        
        if level>=1

    
            for ind=1:bS*4
                small0big1=1;
                [temp_score,temp_index_vectors,temp_perturb_vectors,temp_perturb_param]=learn(ind,tree,treesize,options,small0big1);
                score(ind,:,:)=temp_score;
                min_score(ind)=min(temp_score(temp_score(:)>0));
                index_vectors(ind)=temp_index_vectors;
                perturb_vectors(ind,:,:)=temp_perturb_vectors;
                perturb_param(ind,:,:)=temp_perturb_param;
            end
            
            tempv=1:bS*4;
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
        for indlocation=1:bS
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
global bW; %board width
global bS;%size of board
global b4S;
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
                        %imshow(reshape(abs(imglin),[bW bW])*bS/max(abs(imglin)));
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
                            %imshow(reshape(abs(imglin),[bW bW])*bS/max(abs(imglin)));
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
                %imshow(reshape(abs(imglin),[bW bW])*256/max(abs(imglin)));
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
            %imshow(reshape(abs(imglin),[bW bW])*256/max(abs(imglin)));
        end

    end %if large

    if index_vectors>1000
        perturb_vectors=perturb_vectors(1:1000,:);
        perturb_param=perturb_param(1:1000,:);
        temp_score=temp_score(1:1000,1);
    elseif index_vectors<1000 && index_vectors>0
        perturb_vectors=[perturb_vectors ; zeros(1000-index_vectors,bS)];
        perturb_param=[perturb_param ; zeros(1000-index_vectors,5)];
        temp_score=[temp_score; zeros(1000-index_vectors,1)];
        if max(temp_score)==0
            temp_score(1,1)=256000;
        end
    else
        perturb_vectors= zeros(1000,bS);
        perturb_param=zeros(1000-index_vectors,5);
        temp_score= zeros(1000-index_vectors,1);
        temp_score(1,1)=256000;
    end
    


end


function [options,show_state,show_index]=choose_CornersNN(options,show_state,show_index,tree,treesize )
    %% NN here stands for nearest neighbors
    global bW; %board width
    
    block_pos=zeros(1,2);
    block_pos(1,1:2)=[1 1];    
    block_pos(2,1:2)=[1 bW];
    block_pos(3,1:2)=[bW 1];    
    block_pos(4,1:2)=[bW bW];    
    block_pos(5,1:2)=[1 2];    
    block_pos(6,1:2)=[2 1];    
    block_pos(7,1:2)=[bW 2];    
    block_pos(8,1:2)=[bW-1 1];    
    block_pos(9,1:2)=[1 bW-1];    
    block_pos(10,1:2)=[2 bW];    
    block_pos(11,1:2)=[bW bW-1];    
    block_pos(12,1:2)=[bW-1 bW];    
    maxpos=12;
    shuffled_choices=zeros(maxpos,100);
    count_choices=zeros(maxpos,1);
    indf=1:bW*bW*4;
    base_show_state=show_state;
    base_show_index=show_index;
    base_options=options;
    base_tree=tree;
    for posind=1:maxpos
        row=block_pos(posind,1);
        col=block_pos(posind,2);
        indlocation=(row-1)*bW+col;
        indf_possible=indf(options(indlocation,indf)>0);
        lenvect=length(indf_possible);
        count_choices(posind)=lenvect;
        shuffled_choices(posind,1:lenvect)=indf_possible;
    end
    number_of_config=1;
    for t=1:maxpos
        number_of_config=number_of_config*count_choices(t);
    end
    parallelForce=8;
    best_score=zeros(1,parallelForce);
    best_configno=zeros(1,parallelForce);
    for stepfor=1:parallelForce
        score=0;
        posof_choices=ones(maxpos,1);

        for configno=stepfor:parallelForce:number_of_config
            for topos=1:maxpos
                modulus=count_choices(topos);
                divideby=1;
                for tt=1:topos-1
                    divideby=divideby*count_choices(tt);
                end
                posof_choices(topos)=mod(floor((configno-1)/divideby),modulus)+1;
            end
            toptions=base_options;
            ttree=base_tree;
            for posind=1:maxpos
                row=block_pos(posind,1);
                col=block_pos(posind,2);
                ttindlocation=(row-1)*bW+col;
                toptions(ttindlocation,:)=0;
                toptions(:,shuffled_choices(posind,posof_choices(posind)))=0;
                toptions(ttindlocation,shuffled_choices(posind,posof_choices(posind)))=1;
            end
            timglin=sum(toptions>0,2);
            tshow_state=1*(timglin==1);
            tshow_index=zeros(size(tshow_state));
            for tttindlocation=1:bW*bW
                if tshow_state(tttindlocation)==1
                    tshow_index(tttindlocation)=min(indf(toptions(tttindlocation,:)>0));
                end
            end
    
            ttree=update_tree(base_tree,tshow_state,tshow_index);
            [score,~]=analyze(0,ttree,treesize,toptions); % W
            if score>best_score(stepfor)
                best_score(stepfor)=score;
                best_configno(stepfor)=configno;
            end
        end %for configno

    end %parfor step
    index_scor_vect=1:length(best_score);
    index_scor=min(index_scor_vect(best_score==max(best_score)));
    best_score_final=best_score(index_scor);
    best_configno_final=best_configno(index_scor);

    configno=best_configno_final;

    for topos=1:maxpos
        modulus=count_choices(topos);
        divideby=1;
        for tt=1:topos-1
            divideby=divideby*count_choices(tt);
        end
        posof_choices(topos)=mod(floor((configno-1)/divideby),modulus)+1;
    end
    options=base_options;
    for posind=1:maxpos
        row=block_pos(posind,1);
        col=block_pos(posind,2);
        indlocation=(row-1)*bW+col;
        options(indlocation,:)=0;
        options(:,shuffled_choices(posind,posof_choices(posind)))=0;
        options(indlocation,shuffled_choices(posind,posof_choices(posind)))=1;
    end
    disp(['Showing result of best configuration with score=' string(best_score_final)]);
    imglin=sum(options>0,2);
    show_state=1*(imglin==1);
    show_index=zeros(size(show_state));
    for indlocation=1:bW*bW
        if show_state(indlocation)==1
            show_index(indlocation)=min(indf(options(indlocation,:)>0));
        end
    end

end





function [options,show_state,show_index,failblock]=form_block(options,baseline_options,show_state,show_index,piece,baseline_tree,e2folder,compromize,isprob,h_number,quantum_d_wave,dwave_options)
global bW; %board width
global bS;%size of board
global b4S;
    
    global flagloop;
    global extraline;
    failblock=false;
    prob_options=options;
    if flagloop==false || extraline==0
        if h_number<0
            extraline=input('How many lines in frame ?  ');
        else
            extraline=h_number;
        end
    end

    flag_useHamiltonian=input('Use Hamiltonian server? (0-no,1-yes): ');
    if flag_useHamiltonian
        writematrix('stop',[getenv('Ocean') '\command2server.txt']);
        writematrix('0',[getenv('Ocean') '\replyfromserver.txt']);
        outnumbers_count=1;
        [filename,path] = uigetfile([getenv('Ocean') '\QA\framework*.mat'],'Now run python code and Fetch frameowrk.mat file with node mapping');
        filen=[path filename];
        load(filen); %load NodeSpin, mapnodes_ind, mapnodes_loc, options (not needed)
    end
    options=prob_options;




    sieve_level=2;%input('Sieve level (1-light, 2-serious): ');
    if extraline==8
        flag_skip=true;
    else
        flag_skip=false;
    end

    writelines('start',sprintf('%s\\break.txt',e2folder));
    disp('Starting to solve for frame block ...');
    block_pos=zeros(1,2);
    posind=0;
    for t=1:extraline
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
    for t=extraline+1:bW-extraline
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
    for t=bW+1-extraline:bW
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
    for t=extraline+1:bW-extraline
        c=t;
        v_r=bW+1-extraline:bW;
        if true %mod(c,2)==0
            v_r=v_r(end:-1:1);
        end
        for r=v_r
            posind=posind+1;
            block_pos(posind,1:2)=[r c];
        end
    end
    for t=bW+1-extraline:bW
        c=t;
        v_r=bW+1-extraline:bW;
        if true %mod(c,2)==0
            v_r=v_r(end:-1:1);
        end
        for r=v_r
            posind=posind+1;
            block_pos(posind,1:2)=[r c];
        end
    end
    for t=bW-extraline:-1:extraline+1
        r=t;
        v_c=bW+1-extraline:bW;
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
        v_c=bW+1-extraline:bW;
        if true%mod(r,2)==0
            v_c=v_c(end:-1:1);
        end
        for c=v_c
            posind=posind+1;
            block_pos(posind,1:2)=[r c];
        end
    end
    t_vectorcon=(bW-extraline):-1:(extraline+1);%  [extraline+1:5 (bW-extraline):-1:12 6 11 7 10 8 9];
    t_vectorcon(t_vectorcon<=extraline)=[];
    t_vectorcon(t_vectorcon>=(bW+1)-extraline)=[];
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
    maxpos=posind;
    disp(sprintf('maxpos=%d',maxpos));
    firstopen_posind=-1;
    shuffled_choices=zeros(maxpos,100);
    count_choices=zeros(maxpos,1);
    posof_choices=ones(maxpos,1);
    for posind=1:maxpos
        row=block_pos(posind,1);
        col=block_pos(posind,2);
        indlocation=(row-1)*bW+col;
        if show_state(indlocation)==1
            count_choices(posind)=-1;
            shuffled_choices(posind)=show_index(indlocation);
        elseif firstopen_posind==-1
            firstopen_posind=posind;
        end
    end
    indf=1:b4S;
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
            indlocation=(row-1)*bW+col;
            flag_needrefresh=flag_needrefresh || (count_choices(posind)==0);

            score_nostop=1;
            if flag_needrefresh
                %%%%  board to limit further choices
                keep_state5=(show_state==5);
                show_state= 1*(keep_state5 | base_show_state>0); %use only chosen options or baseline options, not including predictions from previous
                show_index(base_show_state>0)=base_show_index(base_show_state>0);
                show_state(indlocation)=0; %so will take from base (state before starting h), fresh
                keep_state5(indlocation)=0;
                for tindlocation=1:bS
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
                                    if next_row==0 || next_row==(bW+1) || next_col==0 || next_col==(bW+1)
                                        next_sym=0;
                                    else
                                        next_location=(next_col-1)+(next_row-1)*bW+1;
                                        indf=1:b4S;
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
                        indu=1:bS;
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

%Copy to python server
%{
import time
command_file = "C:/Users/seifer/PycharmProjects/Ocean/command2server.txt"
reply_file = "C:/Users/seifer/PycharmProjects/Ocean/replyfromserver.txt"
print('Ready to take tasks')
while True:
    try:
        # Read the command file
        with open(command_file, "r") as f:
            content = f.read().strip()
        if content.lower() == "run":
            filename2='C:/Users/seifer/PycharmProjects/Ocean/QA/bestspin.csv'
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




                        lenvect=length(index_solved_vector);
                        %shuffle vector
                        if lenvect>1
                            for t=1:round(rand*10)
                                randperm(lenvect);
                            end
                        end
                        index_solved_vector=index_solved_vector(randperm(lenvect)); %shuffle vector
        
                        toptions=options;
                        timglin=sum(toptions>0,2);
                        tshow_state=1*(timglin==1);
                        tshow_index=zeros(size(tshow_state));
                        for tindlocation=1:bS
                            if tshow_state(tindlocation)==1
                                tshow_index(tindlocation)=min(indf(toptions(tindlocation,:)>0));
                            end
                        end

                        score_nostop_vector=ones(1,lenvect);
                        if lenvect>=1 && flag_useHamiltonian
                            %Here the python program should be already
                            %running, just communicate with files
                            %use: indchosen, indlocation, options, 
                            flag_startfromDWAVEsolution=false;
                            if flag_startfromDWAVEsolution
                                bestspin=NodeSpin;
                            else
                                bestspin=1*ones(size(NodeSpin));
                                for nd=1:length(NodeSpin)
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
                               filenamebest=[getenv('Ocean') '\QA\bestspin.csv'];
                               writematrix(bestspin,filenamebest);
                               writematrix('run',[getenv('Ocean') '\command2server.txt']);
                               while true
                                    outnumbers=readmatrix([getenv('Ocean') '\replyfromserver.txt']);
                                    if length(outnumbers)>outnumbers_count
                                        outnumbers_count=length(outnumbers);
                                        Hamiltonian=outnumbers(outnumbers_count);
                                        break;
                                    end
                                    pause(0.1);
                               end
                               score_nostop_vector(testchoice)=-Hamiltonian;  %Use in negative: good if score high
                            end
                        end
                        %sort by score (after vector already shuffled)
                        %%score_nostop_vector(score_nostop_vector>0 & score_nostop_vector<=112)=1;
                        %score_nostop_vector=(score_nostop_vector);  
                        [~, index_vector]=sort(score_nostop_vector,'descend');
                        index_solved_vector=index_solved_vector(index_vector);%apply score order on prioirity of piece-index to choose first
                        index_solved_vector=index_solved_vector(index_solved_vector>0);
                        

                    end %end of selection between one option and more

                    lenvect=length(index_solved_vector);
                    if lenvect>0
                        active_posind=posind;
                        count_choices(posind)=lenvect;
                        posof_choices(posind)=1;
                        if quantum_d_wave && count_choices(posind)>1 %define preferred order according to dwave answers
                            deepfill=false;
                            problemsize=2000;
                            [score_dwave, ~]=dwave_solve(indlocation,options,index_solved_vector, show_state,show_index, piece,deepfill,problemsize,0,false );
                            [~, index_vector]=sort(score_dwave,'descend');
                            index_solved_vector=index_solved_vector(index_vector);%apply score order on prioirity of piece-index to choose first
                        elseif ~isempty(dwave_options) && count_choices(posind)>1
                            score_dwave=dwave_options(indlocation,index_solved_vector);
                            [~, index_vector]=sort(score_dwave,'descend');
                            index_solved_vector=index_solved_vector(index_vector);
                        end
                        shuffled_choices(posind,1:lenvect)=index_solved_vector;
                        show_index(indlocation)=shuffled_choices(posind,posof_choices(posind));
                        show_state(indlocation)=5; %mark as planned choices
                        if posind>maxfill_posind
                            maxfill_posind=posind;
                            fprintf(' %d  ',posind);
                            best_options=options;
                            best_show_state=show_state;
                            best_show_index=show_index;
                            total_time=total_time+toc;
                            planv=0; %for virtual puzzle
                            if planv==0
                                required_options=zeros(bS,b4S);
                                for t=1:bS
                                    required_options(t,(t-1)*4+1)=1;
                                end
                                firstcheck_actual_sol=sum(sum((options>0.5).*(required_options>0.5),2));
                                fprintf('(s=%d), ',firstcheck_actual_sol);
                            end

                            tic;
                        end

                        lines=readlines(sprintf('%s\\break.txt',e2folder));
                        if contains(lines(1),'end') 
                            disp('User requested to exit');
                            flagloop=false;
                            failblock=true;
                            options=best_options;
                            show_state=best_show_state;
                            show_index=best_show_index;
                            return;
                        end
                    end
                end   %end of test for available choices
                if lenvect==0 %failed to generate options for new posind, so go back to previous posind and advance choice
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
                %at this point all is find, continue for next posind 
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

    disp('Block ready');


end


function [score_dwave, options]=dwave_solve(indlocation,options,shuffled_choices_pick, show_state, show_index,piece,deepfill,problemsize,filenumber_ready,ancilla_way)
global bW; %board width
global bS;%size of board
global b4S;
global quantum_computation_paid;
global use_only_QPU;

    if length(shuffled_choices_pick)==0
        disp('No fit to this position, try another one');
        indpiece=[];
        return
    end
    max_nodes=problemsize;

    use_hybrid_solver=~use_only_QPU;
    if deepfill==false
        %max_couples=24; %Usal maximum 24,  Advantage_system6.4 can use 34
        if use_hybrid_solver
            max_couples=4000; %hybrid case
        else
            max_couples=23; %Usal maximum 24,  Advantage_system6.4 can use 34
            %max_couples=33; %Usal maximum 24,  Advantage_system6.4 can use 34
        end
    elseif deepfill==true
        if use_hybrid_solver
            max_couples=4000; 
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

        QAcounter=readmatrix([getenv('Ocean') '\QA\counter.txt']);
        QAcounter=QAcounter+1;
        fileID = fopen([getenv('Ocean') sprintf('\\ask_d_wave%g.py',QAcounter)],'w');
        fprintf(fileID,'import numpy as np\n');
        fprintf(fileID,'import networkx as nx\n');
        fprintf(fileID,'import matplotlib.pyplot as plt\n');
        fprintf(fileID,'import dimod\n');
        fprintf(fileID,'import csv\n');
        %fprintf(fileID,'from helpers.draw import plot_bqm \n');
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
                if M>0 && M<=60 %upper limit due to memory restriction
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


        %exclude using the same piece (possibly rotated) at different locations: define inequality in chain
        % For two spins the hamiltonial is (-S1+1)*(-S2+1), so there is both
        % coupling and bias
        if ancilla_way
            excd_penalty_weight=0.0075;
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

        else %if ancilla_way: false


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
    
                end
            end



        end %if ancilla_way


        %NORMALIZE: make the parameters fit direct quantum processor
        %maxJ_vect=sort(abs(mapnodes_couplers_weight(:)));
        %maxH_vect=sort(abs(mapnodes_magnetic(:)));
        %maxJ=maxJ_vect(floor(end*1.000));%+std(maxJ_vect);
        %maxH=maxH_vect(floor(end*1.000));%+std(maxH_vect);
        maxJ=max(abs(mapnodes_couplers_weight(:)));
        maxH=max(abs(mapnodes_magnetic(:)));
        %excessive_factor=max(maxJ,maxH/6);%in most advanced QPU
        excessive_factor=max(maxJ,maxH/2); %usual range in QPU is 4 , but we want advantage for already set pieces
        mapnodes_couplers_weight=sign(mapnodes_couplers_weight).*min(abs(mapnodes_couplers_weight)/excessive_factor,1);
        mapnodes_magnetic=sign(mapnodes_magnetic).*min(abs(mapnodes_magnetic)/excessive_factor,2);
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
            disp(sprintf('Count of nodes: %d (+ %d accessory), board coverage= %d,  Count of faults (duplicate nodes per option): %d ',node_number,ancil_number,sum(checkloc)/bS,countfault));
        end


        ext_node_number=node_number+ancil_number;

        if ext_node_number<=1024
            time_hybrid=5;
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
            time_hybrid=200+(600-200)*(ext_node_number-30000)/(100000-30000)+1;
        end
        if ext_node_number>100000 && ext_node_number<=1000000
            time_hybrid=600;
        end
        %time_hybrid=time_hybrid*2; %!!!!!!!!!!!!!!!!! 
        disp(time_hybrid);
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
            fprintf(fileID,'result = sampler_qpu.sample(bqm)\n');
            fprintf(fileID,'print("Minimum energy= {}".format(result.first.energy))\n');
            fprintf(fileID,'data=result.first.sample\n');
        end

        fprintf(fileID,'print("Best solution found: {}".format(data))\n');
        fprintf(fileID,sprintf('for key in range(1,%d):\n',node_number+ancil_number+1));
        %fprintf(fileID,'     data.setdefault(key, 0)\n');
        fprintf(fileID,'     data.setdefault(key, -1)\n');
        fprintf(fileID,'data_sorted=dict(sorted(data.items()))\n');
        fprintf(fileID,'array = np.array(list(data_sorted.values()))\n');
        fprintf(fileID,sprintf('filename=''%s/QA/d_wave_answer%g.csv''\n',strrep(getenv('Ocean'),'\','/'),QAcounter));
        fprintf(fileID,'np.savetxt(filename, array,fmt=''%%d'',  delimiter='','')\n');
        fclose(fileID);

        writematrix(QAcounter,[getenv('Ocean') '\QA\counter.txt']);
        %Here send execution sign for d-wave , click manually to continue
        disp(sprintf(' task for QM (%d,%d) ',ext_node_number,question_size));
        tempmatfilename=[getenv('Ocean') sprintf('\\QA\\framework%g.mat',QAcounter)];
        save(tempmatfilename,'shuffled_choices_pick','options','mapnodes_loc','mapnodes_ind','-v7.3');

    else
        QAcounter=filenumber_ready;
        tempmatfilename=[getenv('Ocean') sprintf('\\QA\\framework%g.mat',QAcounter)];
        load(tempmatfilename,'shuffled_choices_pick','options','mapnodes_loc','mapnodes_ind');
        node_number=length(mapnodes_loc(mapnodes_ind>0));
    end %if filenumber_ready==0
    
    try
        commandext=[getenv('Ocean') sprintf('\\.venv\\Scripts\\python.exe ') getenv('Ocean') sprintf('\\ask_d_wave%d.py',QAcounter)];
        system(commandext);
        pause(3);

        NodeSpin=readmatrix([getenv('Ocean') sprintf('\\QA\\d_wave_answer%g.csv',QAcounter)], "Delimiter",",");
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
        tempmatfilename=[getenv('Ocean') sprintf('\\QA\\framework%g.mat',QAcounter)];
        save(tempmatfilename,'shuffled_choices_pick','NodeSpin','options','mapnodes_loc','mapnodes_ind','-v7.3');
        disp(sprintf('Sum(options)=%g ',sum(options(:))));
        for t=1:bS
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

