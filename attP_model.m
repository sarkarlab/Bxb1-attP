% "Input_reading": Convert the input sequences from letters to number matrixs
filename = '1 nM.xlsx'; % select the sequence samples
l = 200; % select the numbers of random sequences
L = 3000; % select the number of top sequences considered for data analysis

% Function of "Input_reading"
[M,X,y] = Input_reading(filename,l,L); % M = sequence in matrix format; X= sequences in expanded matrix format; y=experiments freqeuncy result


%randomly spilt data to training set and validating set and repeat the same
%thing 'perm' imes
perm =1; 
W = zeros (3,9,perm);
w = zeros (27, perm);
m=0.7*l; % training data set
n=l-m; % validating data set

for j=1:perm
    [Xtr, ytr, Xcv, ycv]=Data_spliting(X, y, m, n, l);
    
    % try different inital guess values to see if the calculate weigth score converge
    rep = 1; % number of initial value guess
    A = zeros (3,9,rep); 
    a = zeros (27,rep);
    for i=1: rep
    [A(:,:,i),a(:,i)] = Weight_score(Xtr, ytr);   
    end


% validate the derived weight score 'w' and 'W" using the 'cv' data set
[hcv,Dev_cv] = model_validation (Xcv, ycv, a(:,i), n);
[htr,Dev_tr] = model_validation (Xtr, ytr, a(:,i), m);
figure
scatter (ycv,hcv);
hold on
scatter (ytr,htr);

W (:,:,j) = A (:,:,i);
w (:,j)= a (:,i);
end


%% Functions
function [M,X,y] = Input_reading(filename,l,L)
 M = zeros (3,9,l);
 X = zeros (l,27);
 y = zeros (l,1);
 
% Use WT sequence as the template and tranform the sequence to a matrix
ACGTtemplate=xlsread('ACGT_tranformation.xlsx', 'number code of the matrix'); % the number version of the letter matrix based on the parent sequence reference

% randomly read l sequences 'seq' and frequency 'y' from the excel file
rawdata = xlsread(filename, 'Matlab_output');

% overview of the data input
figure
plot (1:L,rawdata(1:L,13)); 

%randomly pick l sequences from top L sequences
k = randperm(L);
seq=rawdata(k(1:l),2:11);
y=rawdata(k(1:l),13);
    
% convert letter sequences to number matrix M and expanded X
   
    for s=1:l
     for i=1:9
            for j=2:4
            
             if seq(s,i)== ACGTtemplate(j,i)
                M (j,i,s)=1; % 0/1 according to the new reference ACGT template
             end
    
            end
     end
    end
    
% convert (4-1)X10 ACGT coding to 1X30 vector  
    for s=1:l
    
        p=1;
        for i=1:9 
         for j=2:4
            X(s,p)= M(j,i,s);
            p=p+1;
         end
        end
    
    end
end
 

% split data set into training data set 'tr' and cross validating data set
% 'cv'
function [Xtr, ytr, Xcv, ycv]=Data_spliting(X, y, m, n, l)
k = randperm(l);
% training data set with m sequences 
Xtr = X(k(1:m),:); 
ytr = y(k(1:m),:);
% valifating data set with m sequences 
Xcv = X(k(m+1:m+n),:); 
ycv = y(k(m+1:m+n),:);
end

% least square fitting to calculate the weight score aij

function [W,w]= Weight_score(X, y)
W = zeros (3,9);

C=X;
d=y;
A=ones(1,27);
b=10;
Aeq=[];
beq=[];
lb=-ones(27,1)*2;
ub=ones(27,1)*2;
x0=randn(27,1); % initial guess
x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0); % calculated weight score (extended)
w = x;
% fold the extended weight score x to a folded matrix Wij
q=1;
for j= 1:9
    for i=1:3
        W(i,j)=x(q);
        q=q+1;
    end
end

end

% model validation

function [h,Dev] = model_validation (X, y, w, n)
h = X*w;
Dev = 0;
    for i = 1: n
    Dev = Dev + (h(i,1)-y(i,1))^2/n;
    end
end
