function similarity_matrix=GetImprovedLNSimilarity(train_interaction_matrix1,feature_matrix,neighbor_alpha1,neighbor_alpha2)
%�ж�x��Ԫ���ǲ���0������Ǳ��Ϊ1��������Ϊ0.���ڲ�����2����ʾ�������������ȫ��Ԫ�ض���1��Ҳ����ԭ����ȫ�ж���0�����зŽ�index����
train_interaction_matrix=train_interaction_matrix1;
index_0_row=find(all(train_interaction_matrix==0,2)==1);

%--------------------------------��ʼ�����ƶȾ���-----------------------------------------------------------------
similarity_matrix=ones(size(train_interaction_matrix, 1));   %��ѯ��һ��ά�ȵĳ���size�����ҽ���һ��size*sizeȫ1�ķ���
%��Ϊ��Щ��ȫ��Ԫ�ض���0�������ڳ�ʼ������ʱҪ���Ѿ���ʼ����1��Ϊ0
for i=1:size(index_0_row, 1) %��һ��ά�ȵĳ���
    similarity_matrix(index_0_row(i, 1), :)=zeros(1, size(train_interaction_matrix, 1));  %���� 1* size��ȫ�����
    similarity_matrix(:, index_0_row(i, 1))=zeros(size(train_interaction_matrix, 1), 1);  %����size * 1��ȫ�����
end

%��佻�������ƶ�
temp_interaction=train_interaction_matrix;     % �������໥���þ���  %7963*1340
temp_interaction(index_0_row, :)=[];           % ԭ����ȫ�ж���0������Ϊ�գ�Ҳ����ȥ����   size(temp_interaction, 1) % 5720 
%round(size(temp_interaction, 1)*neighbor_alpha)  % size(A,1) ������������neighbor_alpha=0.8  ���ʽ���Ǽٶ��ھӵ�����    4576
temp_similarity_matrix=GetLNSimilarity(temp_interaction, round(size(temp_interaction, 1)*neighbor_alpha1),round(size(temp_interaction, 1)*neighbor_alpha1*neighbor_alpha2));
count = 1;
for j=1:size(similarity_matrix(:), 1)   %A(:)  ����A�е�ÿ�кϲ���һ������������  Ȼ��������һ��������
    if(similarity_matrix(j)==1)
        similarity_matrix(j)=temp_similarity_matrix(count);
        count=count+1;
    end
end

%����������ƶ�   % ��������������
feature_similarity_matrix=GetLNSimilarity(feature_matrix, round(size(feature_matrix, 1)*neighbor_alpha1),round(size(feature_matrix, 1)*neighbor_alpha1*neighbor_alpha2));  %7963*0.8=6370

for k=1:size(index_0_row, 1)
    similarity_matrix(index_0_row(k, 1), :)=feature_similarity_matrix(index_0_row(k, 1), :);
end
end

% �õ������ڽ���
function W=GetLNSimilarity(feature_matrix,neighbor_num1,neighbor_num2)
%�����������ƶȵĿ��ټ��㷽��(W==similarity_matrix) �õ����ƶȾ���
iteration_max=40;
mu=3;
X=feature_matrix;
size(X);    

%=======================================================================
%�����Ҿ�������ھ� ��C1
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

%���������ƾ���  �����������ԣ�
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
    G(isnan(G))=0;   %NULLֵȫ����Ϊ0
end
% size(G)   %5720*5720

%----------------------------------------------------------------------
%�õ���G,��Ϊ��ʼֵ��   ��ŷ�Ͼ�������ھ�
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

%���������ƾ���
C2=nearst_neighbor_matrix2;
W=G;
W=(C2.*W);
lamda=mu*e;
P=X*X'+lamda*e';
for i=1:iteration_max
    Q=(C2.*W)*P;
    W=(C2.*W).*P./Q;
    W(isnan(W))=0;  %NULLֵȫ�����0
end
end