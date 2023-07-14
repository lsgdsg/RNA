clear;
clc;
warning('off');
%

currentFolder = pwd;              
addpath(genpath(currentFolder));  
n=1;
m=1;
  y=zeros(10,10);
  w=0;
% Wdd = load('diseasesimilarity.txt');
% A = load('known_lncRNA_disease_interaction.txt');
% Wrr_v=lncRNAfunsim(Wdd,A);
 K=15;
    for arf2=0.1:0.1:1
        w=w+1;
        
        
        %% 
 Wdd = load('DiseaseSimilarityModel1.csv');
    for i=1:length(Wdd)
           Wdd(i,i)=1;
    end
 A = load('DiseaseAndRNABinary.csv');
 A = A';
 Wrr_v = load('B.mat');
 Wrr_v = Wrr_v.B;   
 %% 

 
 
 k1=0.7;
 k2=1;

 L1=10;
 L2=3;
 arf=0.3;
 
 
 t1=15;
 t2=30;
 arf1=0.1;
 
%  K=10;
%  arf2=0.6;
 
A_ori=A;
[nl,nd] = size(A);
 Wrr_v3=lncRNAfunsim(Wdd,A);
[similairty_matrix1,similairty_matrix2] = gaussiansimilarity(A,nl,nd);  
    K1 = [];
    K1(:,:,1)=similairty_matrix1;
    K1(:,:,2)=Wrr_v; 
    K1(:,:,3)=Wrr_v3;
    Wrr_v1=SKF({K1(:,:,1),K1(:,:,2),K1(:,:,3)},t1,t2,arf1);
    
    K2 = [];
    K2(:,:,1)=similairty_matrix2;
    K2(:,:,2)=Wdd;
    Wdd1=SKF({K2(:,:,1),K2(:,:,2)},t1,t2,arf1);
    
y_train=WKNKN( A, Wrr_v1, Wdd1,  K, arf2  ); 

similairty_matrix22=GetImprovedLNSimilarity(A',Wdd1,k1,k2);
similairty_matrix11=GetImprovedLNSimilarity(A,Wrr_v1,k1,k2);


F_1_ori=BR5(y_train,similairty_matrix11,similairty_matrix22,L1,L2,arf);

index=find(A_ori==1);
for u=1:length(index)
    w
    u
     A=A_ori ;
    A(index(u))=0;
  
     Wrr_v4=lncRNAfunsim(Wdd,A);
     [similairty_matrix1,similairty_matrix2] = gaussiansimilarity(A,nl,nd);  
    K1 = [];
    K1(:,:,1)=similairty_matrix1;
    K1(:,:,2)=Wrr_v; 
    K1(:,:,3)=Wrr_v4; 
    Wrr_v2=SKF({K1(:,:,1),K1(:,:,2),K1(:,:,3)},t1,t2,arf1);
    
    K2 = [];
    K2(:,:,1)=similairty_matrix2;
    K2(:,:,2)=Wdd;
    Wdd2=SKF({K2(:,:,1),K2(:,:,2)},t1,t2,arf1);

y_train=WKNKN( A, Wrr_v2, Wdd2, K, arf2 ); 

similairty_matrix222=GetImprovedLNSimilarity(A',Wdd2,k1,k2);
similairty_matrix111=GetImprovedLNSimilarity(A,Wrr_v2,k1,k2);


F_1=BR5(y_train,similairty_matrix111,similairty_matrix222,L1,L2,arf);

    F_1_ori(index(u))=F_1(index(u));
    
    %A=A_ori ;
end
    pre_label_score = F_1_ori(:);
  
    label_y = A_ori(:);
    auc=roc_1(pre_label_score,label_y,'red');
    y(n,m)=auc;
    save loocv_result/y10 y;
    m=m+1;
    end

save loocv_result/y10 y;
  