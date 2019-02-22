%Output to the Command Window
prompt='Enter the name of file ';
%Asking for user input as a string and storing it
str=input(prompt,'s');
%Set the file name

% Check if  File exists.
if isfile(str)
   %load the file
    file =(str);  
else
    % File does not exist.
    fprintf('Invalid: File does not exist'); 
    return;
end;

%Delimiter used for data reading
delimiter='\t';
%Delimiter used for heading reading
delim='\n';

%flags used to determine if a heading is already read or not
flag_J=0.0;
flag_M=0.0;
flag_R=0.0;
flag_E=0.0;
flag_U=0.0;

%Position on the 1st heading
headerlinesIn = 1.0;

%Read the first heading
f=fopen(file);
Heading=textscan(f,'%s',1,'delimiter','\n', 'headerlines',headerlinesIn-1);
fclose(f);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%READING DATA%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Loop to go through all 5 headers
for i=1:5
    if(Heading{1,1}=="[JOINT COORDINATES]")
        %Read after the first heading till the second heading
        Joint_Data= importdata(file,delimiter,headerlinesIn);
        flag_J=1;
    end;
    if(Heading{1,1}=="[MEMBER JOINT CONNECTIVITY]")
        %Read after the second heading till the third heading
        Member_Data = importdata(file,delimiter,headerlinesIn);
        flag_M=1;
    end;
    if(Heading{1,1}=="[REACTIONS AT NODES]")
        %Re-open file to use the textscan command
        f=fopen(file);
        %Due to different data type present,each column is stored in a cell
        Reaction_Data = textscan(f, '%d %d %s','HeaderLines' ,headerlinesIn);
        %close file
        fclose(f);
        flag_R=1;
    end;
    if(Heading{1,1}=="[EXTERNAL FORCES]")
        %Re-open file to use the textscan command
        f=fopen(file);
        %All four column stored in a different cell
        ExternalForce_Data=textscan(f, '%d %d %f %s' ,'HeaderLines' ,headerlinesIn);
        %close file
        fclose(f);
        flag_E=1;
    end;
    if(Heading{1,1}=="[FORCE UNITS]")
        %Open file
        f=fopen(file);
        %Read the units of the force
        Force_Units = textscan(f, '%s\n %*[^\n]','HeaderLines',headerlinesIn);
        %close the file
        fclose(f);
        flag_U=1;
    end;
    
    if(flag_J==1.0)
        %Compute the position of the heading after [JOINTS COORDINATES]
        headerlinesIn=headerlinesIn+size(Joint_Data.data,1)+2.0;
        flag_J=0.0;
    end;
    if(flag_M==1.0)
        %Compute the position of the heading after [MEMBER JOINT CONNECTIVITY]
        headerlinesIn=headerlinesIn+size(Member_Data.data,1)+2.0;
        flag_M=0.0;
    end;
    if(flag_R==1.0)
       %Compute the position of the heading after [REACTION FORCES]
        headerlinesIn=headerlinesIn+size(Reaction_Data{1,3},1)+2.0;
        flag_R=0.0;
    end;
    if(flag_E==1.0)
        %Compute the position of the heading after [EXTERNAL FORCES]
        headerlinesIn=headerlinesIn+size(ExternalForce_Data{1,3},1)+2.0;
        flag_E=0.0;
    end;
    if(flag_U==1.0)
        %Compute the position of the heading after [FORCE UNITS]
        headerlinesIn=headerlinesIn+size(Force_Units{1,1},1)+2.0;
        flag_U=0.0;
    end;
    %Reopen file and read the next heading
    f=fopen(file);
    Heading=textscan(f,'%s',1,'delimiter','\n', 'headerlines',headerlinesIn-1);
    fclose(f);
end;

%number of joints
NJ=size(Joint_Data.data,1);
%number of members
NM=size(Member_Data.data,1);
%number of reactions forces
NR=size(Reaction_Data{1,3},1);
%number of external forces
NX=size(ExternalForce_Data{1,3},1);

if(2.0*NJ==NM+NR)
     % DO NOTHING but continue.
else
    %DATA in file is not correct
    fprintf('Incorrect data detected');
    return;
end;
%make a matrix M and fill it with zeroes
M=zeros(2.0*NJ, 2.0*NJ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%POPULATING THE M MATRIX%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Go through all the members and calulate values needed
for i=1:size(Member_Data.data,1)
    %Retrieves joints to which each member is connected
    joint_from=Member_Data.data(i,2)+1.0;
    joint_to=Member_Data.data(i,3)+1.0;
    %retrieves the joints position in X/Y coordinates
    x1=Joint_Data.data(joint_from,2);
    y1=Joint_Data.data(joint_from,3);
    x2=Joint_Data.data(joint_to,2);
    y2=Joint_Data.data(joint_to,3);
    
    %Calculates length of each member
    Length=sqrt((x2-x1)^2+(y2-y1)^2);
    
    %Populates Member colums of M matrix with fractional values
    M(2.0*joint_from-1.0,i) = (x2-x1)/Length;
    M(2.0*joint_to-1.0,  i) = -(x2-x1)/Length;
    M(2.0*joint_from,  i) = (y2-y1)/Length;
    M(2.0*joint_to,    i) = -(y2-y1)/Length;
end;

%Go through reaction forces to populate second half of M matrix
for j=1:NR
    %Retrieves reaction data the node and the direction of force
    Node=Reaction_Data{1,2}(j,1);
    Direction = Reaction_Data{1,3}(j,1);
    
    %Checks direction and puts it in appropriate M matrix spot
    if((Direction == "y") || (Direction == "Y"))
        M(2.0*Node+2, NM+j)=M(2.0*Node+2.0, NM+j)+1.0;
    elseif ((Direction == "x")||(Direction == "X"))
        M(2.0*Node+1, NM+j)=M(2.0*Node+1.0, NM+j)+1.0;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%POPULATING THE External MATRIX%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creates external forces matrix
External=zeros(2.0*NJ,1);

%Goes through external forces data to poulare External matrix
for k=1:NX
    %Retrieves External Forces data
    Node_j =ExternalForce_Data{1,2}(k,1);
    Force=ExternalForce_Data{1,3}(k,1);
    Direction_1 = ExternalForce_Data{1,4}(k,1);
    
    %Checks if forces given in X Y or with an angle
    if(Direction_1~="X"&&Direction_1~="x"&&Direction_1~="Y"&&Direction_1~="y")
        %If angle is given calculate X and Y values and puts them in
        %External matrix accordingly
        Theta = str2double(Direction_1{1,1});
        External((2.0*Node_j+1.0))=External((2.0*Node_j+1.0))-Force*cos(Theta*(pi/180.0));
        External((2.0*Node_j+2.0))=External((2.0*Node_j+2.0))-Force*sin(Theta*(pi/180.0));
    else
        %If given in X Y checks direction and populates External matrix
        %accordingly
        if(Direction_1=="X"||Direction_1=="x")
            External(2.0*Node_j+1.0)=Force;
        else
            External(2.0*Node_j+2.0)=Force;
        end;
    end;
end;

%Solve the system of equations
A=-M\External;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Creeating Output File%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Output to the Command Window
prompt='Enter the name of file to save to: ';
%Asking for user input as a string and storing it
str_1=input(prompt,'s');
%Set the file name

%Open output file and write the matrix header
file=fopen(str_1,'w');
fprintf(file,'M MATRIX\n');
fprintf(file,'--------\n');

%Sets precision of output file to 2 decimal places
digits(2);

%Cycles through M matrix and writes each value to output file
for a=1:2*NJ
    %First for loop cycles through rows and nested cycles through colums
    for b=1:2*NJ
        fprintf(file,'%+.2f ', M(a,b));
    end
    %Goes to new line at the end of each row of the matrix in the output
    %file
    fprintf(file,'\n');
end

%Writes solution header to output file
fprintf(file, '\n\nSOLUTION\n');
fprintf(file,'--------\n');

%converts units taken from input file to print in output file
unit = string(Force_Units{1,1}(1,1));

%Cycles through the Result matrix and writes member forces to output file
for c=1:NM
    fprintf(file, 'Member %02d: F= %8.2f %s ', (c-1), abs(A(c)), unit);
    
    %Checks if force is positive or negative and chooses [C] ot [T]
    %accordingly
    if(A(c)<0)
        fprintf(file, '[C]\n');
    else
        fprintf(file, '[T]\n');
    end
end

%Writes reaction forces to output file
%For loops cycle tthrgough M matrix starting with the first reaction forces
%colum checking for a 1 or a 0
for d=1:2*NJ
    for e=NM+1:2*NJ
        %If M m atrix element is not zero writes reaction force to out file
        %from result matrix
        if(M(d,e) ~= 0)
            %Checks row index, if odd Xdirection, if even Y direction
            %((d-x)/2), x being 1 or 2 to get joint number depending on X
            %or Y direction, x = 1 for X, x = 2 for Y
            if(mod(d,2))
                fprintf(file, 'X-direction reaction on joint %d = %+9.2f %s\n', ((d-1.0)/2.0), A(e), unit);
            else
                fprintf(file, 'Y-direction reaction on joint %d = %+9.2f %s\n', ((d-2.0)/2.0), A(e), unit);
            end
        end
    end
end
   
   
   
   
   
   
   
   
   
   
   
   
   
  
