close all
clear
clc

L=7.5;                                           % Fin Length %
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
a=sparse(a);
b=zeros((L/h)-1,1)-0.1856;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0148;                % Boundary Condition Effect on L-1th nod %
y=cgs(a,b,1e-6,50);                              % Conjugate Gradients Squared Method %
y((L/h),:)=0;                                 
for l=1:(L/h)-1
    y((L/h)+1-l,1)=y((L/h)-l,1);                 % Ordering Output %  
end
y(1,:)=1.2;                                      % Boundary Condition on 1st nod %
y((L/h)+1,:)=0.9852*y((L/h),1)+0.0148;           % Boundary Condition on Last nod %
x=linspace(0,L,(L/h)+1);
plot(x,y,'r');
hold on
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
a=sparse(a);
b=zeros((L/h)-1,1)-0.1856;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0148;                % Boundary Condition Effect on L-1th nod %
y=gmres(a,b);                                    % Generalized Minimum Residual Method %
y((L/h),:)=0;                                 
for l=1:(L/h)-1
    y((L/h)+1-l,1)=y((L/h)-l,1);                 % Ordering Output %  
end
y(1,:)=1.2;                                      % Boundary Condition on 1st nod %
y((L/h)+1,:)=0.9852*y((L/h),1)+0.0148;           % Boundary Condition on Last nod %
x=linspace(0,L,(L/h)+1);
plot(x,y,'b');
hold on
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
a=sparse(a);
b=zeros((L/h)-1,1)-0.1856;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0148;                % Boundary Condition Effect on L-1th nod %
y=lsqr(a,b,1e-6,50);                             % Least Squares Method %
y((L/h),:)=0;                                 
for l=1:(L/h)-1
    y((L/h)+1-l,1)=y((L/h)-l,1);                 % Ordering Output %  
end
y(1,:)=1.2;                                      % Boundary Condition on 1st nod %
y((L/h)+1,:)=0.9852*y((L/h),1)+0.0148;           % Boundary Condition on Last nod %
x=linspace(0,L,(L/h)+1);
plot(x,y,'k');
hold on
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
a=sparse(a);
b=zeros((L/h)-1,1)-0.1856;                       % Sourse Terms %
b(1,1)=b(1,1)-1.2;                               % Boundary Condition Effect on 2nd nod %
b((L/h)-1,1)=b((L/h)-1,1)-0.0148;                % Boundary Condition Effect on L-1th nod %
y=minres(a,b,1e-6,50);                           % Minimum Residual Method %
y((L/h),:)=0;                                 
for l=1:(L/h)-1
    y((L/h)+1-l,1)=y((L/h)-l,1);                 % Ordering Output %  
end
y(1,:)=1.2;                                      % Boundary Condition on 1st nod %
y((L/h)+1,:)=0.9852*y((L/h),1)+0.0148;           % Boundary Condition on Last nod %
x=linspace(0,L,(L/h)+1);
plot(x,y,'c');
xlabel(' Fin Length ');
ylabel(' Temprature ');
legend(' cgs  ',' gmres ',' lsqr ',' minres ');
%%
