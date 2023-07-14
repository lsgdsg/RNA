function cosl = Lcos(interaction)
[m,n] = size(interaction);
% cosl1 =zeros(n,n);
cosl0 = zeros(n,n);
interaction_copy = interaction;
for i=1:n
    for j=i:n
        b=interaction_copy(i,:).*interaction_copy(j,:);
        b1=sum(b);
        c=sqrt(sum(interaction_copy(i,:).^2));
        c1=sqrt(sum(interaction_copy(j,:).^2));
        c2=c*c1;
        cosl0(i,j)=b1/c2;
    end
end
cosl = cosl0 + cosl0' - diag(diag(cosl0));
end


