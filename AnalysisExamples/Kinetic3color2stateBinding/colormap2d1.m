function map2d=colormap2d1(nrb)

rbmap = zeros(nrb*7-1,3);
rbmap(1:nrb,3) = (0.5+1/2/nrb:1/nrb/2:1)';
rbmap(nrb+1:3*nrb,3) = 1;
rbmap(nrb+1:3*nrb,2) = (1/2/nrb:1/nrb/2:1)';
rbmap(3*nrb+1:4*nrb,3) = flipud(rbmap(1:nrb,3));
rbmap(3*nrb+1:4*nrb,2) = 1;
rbmap(3*nrb+1:4*nrb,1) = rbmap(1:nrb,3);
rbmap(4*nrb+1:7*nrb-1,:)=fliplr(flipud(rbmap(1:3*nrb-1,:)));

rbmap(1:nrb*3,:)=[];
map2d=rbmap;