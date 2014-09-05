function [ corW, Z ] = clustering_W( inW, mLen, hei, wid)
%clustering_W hierarchical clustering on w

%%Computing the correlation between W
corW = zeros( mLen, mLen );
for i = 1:mLen
    for j = 1:mLen
        temp = corrcoef( inW(i,:), inW(j,:));
        corW(i,j) = temp(1,2);
    end
end
figure(); imagesc(corW);

Z = linkage(inW(:,:),'weighted','correlation');
figure();dendrogram(Z);
c = cluster(Z,'maxclust',4);
id = [2 15 26 5 14 18 23 29];
for i=1:length(id)
    temp1 = reshape(inW(:,:),hei, wid);
    subplot(2,4,i);imagesc(temp1);colorbar;colormap(bone);title(['Distribution of Element' num2str(id(i))]);
end
% id = [10 13 18 26];
% for i=1:length(id)
%     temp1 = reshape(curW(id(i),:,:),hei, wid);
%     subplot(2,2,i);imagesc(temp1);colorbar;colormap(bone);title(['Distribution of Element' num2str(id(i))]);
% end
% id = [11 12 3 16 27];
% for i=1:length(id)
%     temp1 = reshape(curW(id(i),:,:),hei, wid);
%     subplot(2,3,i);imagesc(temp1);colorbar;colormap(bone);title(['Distribution of Element' num2str(id(i))]);
% end
% id = [1 5 2 30 4 25];
% for i=1:length(id)
%     temp1 = reshape(curW(id(i),:,:),hei, wid);
%     subplot(2,3,i);imagesc(temp1);colorbar;colormap(bone);title(['Distribution of Element' num2str(id(i))]);
% end
% id = [21 28];
% for i=1:length(id)
%     temp1 = reshape(curW(id(i),:,:),hei, wid);
%     subplot(1, 2,i);imagesc(temp1);colorbar;colormap(bone);title(['Distribution of Element' num2str(id(i))]);
% end
% id = [23 17];
% for i=1:length(id)
%     temp1 = reshape(curW(id(i),:,:),hei, wid);
%     subplot(1, 2,i);imagesc(temp1);colorbar;colormap(bone);title(['Distribution of Element' num2str(id(i))]);
% end
% 
% id = [15 20 29];
% for i=1:length(id)
%     temp1 = reshape(curW(id(i),:,:),hei, wid);
%     subplot(2,2,i);imagesc(temp1);colorbar;colormap(bone);title(['Distribution of Element' num2str(id(i))]);
% end


end

