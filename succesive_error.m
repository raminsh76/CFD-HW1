close all
clear
clc

L=7.5;                                           % Fin Length %
h=0.375;                                         % Grid Size %
a=zeros((L/h)-1); 
for i=1:(L/h)-1
    for j=1:(L/h)-1
if i==j
    a(i,j)=-6.6406;
    a(i,j+1)=1;                                  % Sparse coefficients Matrix %
    a(i+1,j)=1;
end
   end
end
a(:,(L/h))=[];
a((L/h),:)=[];
a((L/h)-1,(L/h)-1)=a((L/h)-1,(L/h)-1)+0.9302;    % Neumann Boundary condition for Coefficient Matrix %
a=sparse(a);
b=zeros((L/h)-1,1)-4.6406;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0697;                % Boundary Condition Effect on L-1th nod %
y1=cgs(a,b,1e-6,50);                             % Conjugate Gradients Squared Method %
y1((L/h),:)=0;                                
for l=1:(L/h)-1
    y1((L/h)+1-l,1)=y1((L/h)-l,1);               % Ordering Output %  
end
y1(1,:)=1.2;                                     % Boundary Condition on 1st nod %
y1((L/h)+1,:)=0.9302*y1((L/h),1)+0.0697;         % Boundary Condition on Last nod %
x1=linspace(0,L,(L/h)+1);
%%
h=0.15;                                          % Grid Size %
a=zeros((L/h)-1);
for i=1:(L/h)-1
    for j=1:(L/h)-1
if i==j
    a(i,j)=-2.7425;
    a(i,j+1)=1;                                  % Sparse coefficients Matrix %
    a(i+1,j)=1;
end
   end
end
a(:,(L/h))=[];
a((L/h),:)=[];
a((L/h)-1,(L/h)-1)=a((L/h)-1,(L/h)-1)+0.9709;    % Neumann Boundary condition for Coefficient Matrix %
b=zeros((L/h)-1,1)-0.7425;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0291;                % Boundary Condition Effect on L-1th nod %
y2=cgs(a,b,1e-6,50);                             % Conjugate Gradients Squared Method %
y2((L/h),:)=0;                                 
for l=1:(L/h)-1
    y2((L/h)+1-l,1)=y2((L/h)-l,1);               % Ordering Output %  
end
y2(1,:)=1.2;                                     % Boundary Condition on 1st nod %
y2((L/h)+1,:)=0.9709*y2((L/h),1)+0.0291;         % Boundary Condition on Last nod %
x2=linspace(0,L,(L/h)+1);
%%
h=0.075;                                         % Grid Size %
a=zeros((L/h)-1);
for i=1:(L/h)-1
    for j=1:(L/h)-1
if i==j
    a(i,j)=-2.1856;
    a(i,j+1)=1;                                  % Sparse coefficients Matrix %
    a(i+1,j)=1;
end
   end
end
a(:,(L/h))=[];
a((L/h),:)=[];
a((L/h)-1,(L/h)-1)=a((L/h)-1,(L/h)-1)+0.9852;    % Neumann Boundary condition for Coefficient Matrix %
b=zeros((L/h)-1,1)-0.1856;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0148;                % Boundary Condition Effect on L-1th nod %
y3=cgs(a,b,1e-6,50);                             % Conjugate Gradients Squared Method %
y3((L/h),:)=0;                                 
for l=1:(L/h)-1
    y3((L/h)+1-l,1)=y3((L/h)-l,1);               % Ordering Output %  
end
y3(1,:)=1.2;                                     % Boundary Condition on 1st nod %
y3((L/h)+1,:)=0.9852*y3((L/h),1)+0.0148;         % Boundary Condition on Last nod %
x3=linspace(0,L,(L/h)+1);
%%
h=0.0375;                                        % Grid Size %
a=zeros((L/h)-1);
for i=1:(L/h)-1
    for j=1:(L/h)-1
if i==j
    a(i,j)=-2.0464;
    a(i,j+1)=1;                                  % Sparse coefficients Matrix %
    a(i+1,j)=1;
end
   end
end
a(:,(L/h))=[];
a((L/h),:)=[];
a((L/h)-1,(L/h)-1)=a((L/h)-1,(L/h)-1)+0.9926;    % Neumann Boundary condition for Coefficient Matrix %
b=zeros((L/h)-1,1)-0.0464;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0074;                % Boundary Condition Effect on L-1th nod %
y4=cgs(a,b,1e-6,50);                             % Conjugate Gradients Squared Method %
y4((L/h),:)=0;                                 
for l=1:(L/h)-1
    y4((L/h)+1-l,1)=y4((L/h)-l,1);               % Ordering Output %  
end
y4(1,:)=1.2;                                     % Boundary Condition on 1st nod %
y4((L/h)+1,:)=0.9926*y4((L/h),1)+0.0074;         % Boundary Condition on Last nod %
x4=linspace(0,L,(L/h)+1);
%%
y120=y2(1,1)-y1(1,1);
y230=y3(1,1)-y2(1,1);       % for x=0 %
y340=y4(1,1)-y3(1,1);
%%
y121=y2(6,1)-y1(3,1);
y231=y3(11,1)-y2(6,1);       % for x=0.75 %
y341=y4(21,1)-y3(11,1);
%%
y122=y2(11,1)-y1(5,1);
y232=y3(21,1)-y2(11,1);       % for x=1.5 %
y342=y4(41,1)-y3(21,1);
%%
y123=y2(16,1)-y1(7,1);
y233=y3(31,1)-y2(16,1);       % for x=2.25 %
y343=y4(61,1)-y3(31,1);
%%
y124=y2(21,1)-y1(9,1);
y234=y3(41,1)-y2(21,1);       % for x=3 %
y344=y4(81,1)-y3(41,1);
%%
y125=y2(26,1)-y1(11,1);
y235=y3(51,1)-y2(26,1);       % for x=3.75 %
y345=y4(101,1)-y3(51,1);
%%
y126=y2(31,1)-y1(13,1);
y236=y3(61,1)-y2(31,1);       % for x=4.5 %
y346=y4(121,1)-y3(61,1);
%%
y127=y2(36,1)-y1(15,1);
y237=y3(71,1)-y2(36,1);       % for x=5.25 %
y347=y4(141,1)-y3(71,1);
%%
y128=y2(41,1)-y1(17,1);
y238=y3(81,1)-y2(41,1);       % for x=6 %
y348=y4(161,1)-y3(81,1);
%%
y129=y2(46,1)-y1(19,1);
y239=y3(91,1)-y2(46,1);       % for x=6.75 %
y349=y4(181,1)-y3(91,1);
%%
y1210=y2(51,1)-y1(21,1);
y2310=y3(101,1)-y2(51,1);       % for x=6.75 %
y3410=y4(201,1)-y3(101,1);
%%
y_12=[y120 y121 y122 y123 y124 y125 y126 y127 y128 y129 y1210];
y_23=[y230 y231 y232 y233 y234 y235 y236 y237 y238 y239 y2310];
y_34=[y340 y341 y342 y343 y344 y345 y346 y347 y348 y349 y3410];
x_err=[0 1 2 3 4 5 6 7 8 9 10];
plot(x_err,y_12);
hold on
plot(x_err,y_23);
hold on
plot(x_err,y_34);
legend(' succesive error between 20 & 50 nods ',' succesive error between 50 & 100 nods ',' succesive error between 100 & 200 nods ')