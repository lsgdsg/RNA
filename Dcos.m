function cosd = Dcos(interaction)
[m,n] = size(interaction);
% cosd =zeros(m,m);
cosd0 = zeros(m,m);
interaction_copy = interaction;
for i=1:m
    for j=i:m
        b=interaction_copy(i,:).*interaction_copy(j,:);
        b1=sum(b);
        c=sqrt(sum(interaction_copy(i,:).^2));
        c1=sqrt(sum(interaction_copy(j,:).^2));
        c2=c*c1;
        cosd0(i,j)=b1/c2;
    end
end
cosd = cosd0 + cosd0' - diag(diag(cosd0));
end


