function similarity_matrix=GetImprovedLNSimilarity(train_interaction_matrix1,feature_matrix,neighbor_alpha1,neighbor_alpha2)
%判断x的元素是不是0，如果是标记为1，否则标记为0.由于参数是2，表示从行来看，如果全行元素都是1，也就是原矩阵全行都是0，该行放进index里面
train_interaction_matrix=train_interaction_matrix1;
index_0_row=find(all(train_interaction_matrix==0,2)==1);

%--------------------------------初始化相似度矩阵-----------------------------------------------------------------
similarity_matrix=ones(size(train_interaction_matrix, 1));   %查询第一个维度的长度size，并且建立一个size*size全1的方阵。
%因为这些行全部元素都是0，所以在初始化矩阵时要把已经初始化的1改为0
for i=1:size(index_0_row, 1) %第一个维度的长度
    similarity_matrix(index_0_row(i, 1), :)=zeros(1, size(train_interaction_matrix, 1));  %生成 1* size的全零矩阵
    similarity_matrix(:, index_0_row(i, 1))=zeros(size(train_interaction_matrix, 1), 1);  %生成size * 1的全零矩阵
end

%填充交互谱相似度
temp_interaction=train_interaction_matrix;     % 这里用相互作用矩阵  %7963*1340
temp_interaction(index_0_row, :)=[];           % 原矩阵全行都是0的行设为空，也就是去掉了   size(temp_interaction, 1) % 5720 
%round(size(temp_interaction, 1)*neighbor_alpha)  % size(A,1) 求矩阵的行数，neighbor_alpha=0.8  这个式子是假定邻居的数量    4576
temp_similarity_matrix=GetLNSimilarity(temp_interaction, round(size(temp_interaction, 1)*neighbor_alpha1),round(size(temp_interaction, 1)*neighbor_alpha1*neighbor_alpha2));
count = 1;
for j=1:size(similarity_matrix(:), 1)   %A(:)  矩阵A中的每列合并成一个长的列向量  然后又求了一共多少行
    if(similarity_matrix(j)==1)
        similarity_matrix(j)=temp_similarity_matrix(count);
        count=count+1;
    end
end

%填充特征相似度   % 这里用特征矩阵
feature_similarity_matrix=GetLNSimilarity(feature_matrix, round(size(feature_matrix, 1)*neighbor_alpha1),round(size(feature_matrix, 1)*neighbor_alpha1*neighbor_alpha2));  %7963*0.8=6370

for k=1:size(index_0_row, 1)
    similarity_matrix(index_0_row(k, 1), :)=feature_similarity_matrix(index_0_row(k, 1), :);
end
end

% 得到线性邻近度
function W=GetLNSimilarity(feature_matrix,neighbor_num1,neighbor_num2)
%线性邻域相似度的快速计算方法(W==similarity_matrix) 得到相似度矩阵
iteration_max=40;
mu=3;
X=feature_matrix;
size(X);    

%=======================================================================
%用余弦距离计算邻居 求C1
row_num=size(X,1);
distance_matrix=pdist2(X,X,'cosine');
e=ones(row_num,1);
distance_matrix=distance_matrix+diag(e*inf);
[~, si]=sort(distance_matrix,2,'ascend');
nearst_neighbor_matrix=zeros(row_num,row_num);
index=si(:,1:neighbor_num1);
for i=1:row_num
    nearst_neighbor_matrix(i,index(i,:))=1;
end

%迭代求相似矩阵  （余弦相似性）
C1=nearst_neighbor_matrix;
% rand('state',1);
G=rand(row_num,row_num);
%  rand('seed',sum(100*clock));
G=(C1.*G);
lamda=mu*e;
P=X*X'+lamda*e';
for i=1:iteration_max
    Q=(C1.*G)*P;
    G=(C1.*G).*P./Q;
    G(isnan(G))=0;   %NULL值全部变为0
end
% size(G)   %5720*5720

%----------------------------------------------------------------------
%得到了G,作为初始值咯   用欧氏距离计算邻居
X=G*X;
row_num=size(X,1);    %5720
distance_matrix=pdist2(X,X,'euclidean');
e=ones(row_num,1);
distance_matrix=distance_matrix+diag(e*inf);
[~, si]=sort(distance_matrix,2,'ascend');
nearst_neighbor_matrix2=zeros(row_num,row_num);
index=si(:,1:neighbor_num2);
for i=1:row_num
    nearst_neighbor_matrix2(i,index(i,:))=1;
end

%迭代求相似矩阵
C2=nearst_neighbor_matrix2;
W=G;
W=(C2.*W);
lamda=mu*e;
P=X*X'+lamda*e';
for i=1:iteration_max
    Q=(C2.*W)*P;
    W=(C2.*W).*P./Q;
    W(isnan(W))=0;  %NULL值全部变成0
end
end