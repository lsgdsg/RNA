function [Rt]=BR5(A,circSimi,disSimi,L1,L2,alpha)
normWrr = normFun1(circSimi,0.5);
normWdd = normFun1(disSimi,0.5);
% normWrr= matrix_normalize(circSimi);
% normWdd= matrix_normalize(disSimi);

   [rows,cols]=size(A);

A_ori=A;
%% 

for i=1:length(circSimi)
    circSimi(i,i)=1;
end
for j=1:length(disSimi)
    disSimi(j,j)=1;
end

 
     % 避免假阴性数据对矩阵分解的影响，因此对邻接矩阵进行重构
   

   A1=circSimi*A;
    A2=A*disSimi;
    for i=1:rows
        for j=1:cols
            A(i,j)=max(A1(i,j),A2(i,j));
        end
    end
%     A=(A1+A2)/2;
    % 归一化
    AA=zeros(rows,cols);
    for j =1:cols
        colList=A(:,j);
        for i =1:rows
            if max(colList)==0
                AA(i,j)=0;
             else
        AA(i,j)=(A(i,j)-min(colList))/(max(colList)-min(colList));
            end
        end
    end
    A=AA;
    for j =1:cols
        for i =1:rows
            A(i,j)=max(A_ori(i,j),A(i,j));
        end
    end
  

    %% 
   

Rtr0=A/sum(A(:));
Rtd0=Rtr0;
Rtr=Rtr0;
Rtd=Rtr0;
R0=Rtr;

% %bi-random walk on the heterogeneous network
for t=1:max(L1,L2)
    
    ftl = 0;
    ftr = 0;
    
    %random walk on the lncRNA similarity network
    if(t<=L1)
        nRtleft =(1- alpha) * normWrr * Rtr + alpha*Rtr0;
        ftl = 1;
    end
     Rtr = nRtleft;
    %random walk on the disease similarity network
    if(t<=L2)
        nRtright = (1-alpha) *  Rtd * normWdd + alpha*Rtd0;
        ftr = 1;
    end
    Rtd = nRtright;
    %Rt: predictive association scores between each lncRNA-disease pair
   
  
end
% Rtd=norm1(Rtd);
% Rtr=norm1(Rtr);
   Rt = (ftl*Rtr + ftr*Rtd)/(ftl + ftr);
   

% for t=1:max(L1,L2)
%     
%     ftl = 0;
%     ftr = 0;
% 
%     if(t<=L1)
%         Rt1 = Rtr;
%          nRtleft1 =(1-alpha) * (normWrr * Rt1) + alpha*R0;
%         for i=1:2
%             nRtleft0 = (normWrr * Rt1)  ;
%             ftl = 1;
%             Rt1 = nRtleft0;
%         end
%         nRtleft2 =(1-alpha) * nRtleft0 + alpha*R0;
%         nRtleft = (nRtleft1+nRtleft2)/2;
%     end
%     Rtr=nRtleft;
% 
%     if(t<=L2)
%         Rt2 = Rtd;
%           nRtright1 = (1-alpha)* (Rt2 * normWdd) + alpha*R0;
%         for i=1:2
%           nRtright0 =  (Rt2 * normWdd) ;
%           ftr = 1;
%           Rt2=nRtright0;
%         end
%         nRtright2 = (1-alpha)* nRtright0 + alpha*R0;
%         nRtright = (nRtright1 + nRtright2)/2;
%     end
%   Rtd=  nRtright;
%   
%     
% end
%    Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);
end

