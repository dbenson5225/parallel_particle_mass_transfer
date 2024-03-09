%  starting from scratch to verify some things ...
clear all;

%gpu=gpuDeviceCount;      %  See if there's a gpu
gpu=0;                   %  Turn off gpu
%define numerical and physical parameters
np=40000; ens=1;  % # A points, #particles per point, # runs, point radius
nspec=3; ndim=2; eps=1e-10; ds=1;
for i=1:ndim
   xmin(i)=0; xmax(i)=10; 
   ds=ds*(xmax(i)-xmin(i));
end
ds=ds/np;
ttot=100; dt=10; ntsteps=ttot/dt;
nsource=1;
Dtot=1e-3;  kappa=1;  DRW=(1-kappa)*Dtot; DMT=kappa*Dtot;  % Diffusion stufff
beta = 1;    % SPH beta
sqnp=floor(sqrt(np)); np=sqnp*sqnp;  % For debug only
Sx=[1 5 9 15]
%Sx=[ 37 31 25 19 13 9 5 3 2 ]  %  # subdomains in one direction
%Sx=31
for nn=1:length(Sx);       % Go through a bunch of discretizations
constnx=Sx(nn);
A=zeros(np,6+nspec);   % x,y,z,i,j,k,A,B,C ...  
                       % Hate to hard-wire dimensions but it is much easier
A(:,4:6)=1;            % start out the i,j,k, so chem species are 7:end                     
if(gpu>0) 
    A=gpuArray(A);
    xmin=gpuArray(xmin); xmax=gpuArray(xmax);
end

% First do this all in one script ...

% Define ICs %%%%%%%%%%%
for i=1:ndim
  A(1:np,i)=(xmax(i)-xmin(i))*rand(np,1)+xmin(i);
end
%sqnp=floor(sqrt(np)); xvec=linspace(xmin(1)+eps,xmax(2)-eps,sqnp); 
%yvec=linspace(xmin(2)+eps,xmax(2)-eps,sqnp);  np=sqnp*sqnp;
%for i=1:sqnp
%A((i-1)*sqnp+1:i*sqnp,1)=xvec; A((i-1)*sqnp+1:i*sqnp,2)=yvec(i);
%end

% locations of A and B ICs;
concfactor=prod(xmax-xmin)/np;
A(1:nsource,6+1)=1/nsource; A(1:nsource,6+2)=1/nsource/concfactor; % random x,y
A(1:nsource,1:ndim)=0.51*(xmax-xmin); A(1:nsource,7:8)=1/nsource/concfactor;   % Fixed x,y
%%%%%%%%%%%%%%%%%%%%%%%%

%Start a time loop  %%%%%%%%%%%%%%

tic
prefactor=ds/((4*pi*DMT*dt/beta)^(ndim/2));
for n=1:ntsteps
    %ghostmass=A(:,7:end);  
    sigma = sqrt(4*DMT*dt);  % Char. diffusion dist.
    % Do the mass-transfer
    % First chop up space: ideally each chunk is 3*sigma
    nxmt=ones(1,3);  if(gpu>0); nxmt=gpuArray(nxmt); end;
  for i=1:ndim
    pad(i)=4*sigma; if(gpu>0); pad=gpuArray(pad); end;
    %nxmt(i)=floor((xmax(i)-xmin(i))/dxmt(i)); 
    % manual override
    nxmt(i)=constnx;  % subdomains in each direction
    dxmt(i)=(xmax(i)-xmin(i))/nxmt(i);  % for constant dx
    for m=1:nxmt(i)    % This can allow variable dx's later (i.e. kdtree)
        xlo(i,m)=xmin(i)+(m-1)*dxmt(i);
        xhi(i,m)=xmin(i)+m*dxmt(i);
    end
  end
  pad(3)=1;
  disp(['number of subdomains is ',num2str(prod(nxmt)),...
        '; tstep ',num2str(n),' of ',num2str(ntsteps), ...
        '; mass = ',num2str(concfactor*sum(A(:,7)))]);
  
  % Figure out the i,j,k (cell) of each particle
  for i=1:ndim
      A(1:np,3+i)=ceil((A(:,i)-xmin(i))./dxmt(i));
  end
  
  % Build matrix of separations including ghosts from adjoining cells
  for i=1:ndim
      sz(i)=nxmt(i);
  end
  massnew=A(:,7:end);   % Pull masses from A; push to massnew!
  masstemp=massnew;     % This holds temporary masses (incl. ghosts)
  lidx=1:prod(nxmt);   
  [is js ks]=ind2sub(sz,lidx);

  for m=1:prod(nxmt)            % go to every cell  (cell loop)
      ijk(1)=is(m); ijk(2)=js(m); ijk(3)=ks(m);   
      lims=zeros(3,2);  if(gpu>0); lims=gpuArray(lims); end;
      for i=1:3
         lims(i,1)=xlo(1,ijk(i))-pad(i); lims(i,2)=xhi(2,ijk(i))+pad(i);  %xlo and xhi
      end
%      xsmall=xlo(1,inow);  xlarge=xhi(1,inow);
%      ysmall=xlo(2,jnow);  ylarge=xhi(2,jnow);
%      zsmall=xlo(3,know);  zlarge=xhi(3,know);
      idxlocal=find(A(:,4)==ijk(1) & A(:,5)==ijk(2) & A(:,6)==ijk(3));

      idxghost=find(A(:,1)>lims(1,1) & A(:,1)<=lims(1,2) & ...
               A(:,2)>lims(2,1) & A(:,2)<=lims(2,2) & ... 
               A(:,3)>lims(3,1) & A(:,3)<=lims(3,2));
if n==0
figure(44)
plot(A(idxghost,1),A(idxghost,2),'+')
hold on; axis([0 10 0 10]); 
xticks(0:2:10); 
yticks([0:2:10]); 
axis square; grid on;
plot(A(idxlocal,1),A(idxlocal,2),'.')
drawnow; hold off; pause;
end

  if(length(idxlocal)>0)
      X=A(idxghost,:); % npad=sum(idxghost);
      [idx,r]=rangesearch(A(idxghost,1:ndim),A(idxlocal,1:ndim),4*sigma);
      for i=1:length(idxlocal)
          partnow=idxlocal(i);
          idxnow=idx{i};
          numnear=length(idxnow);
          if numnear>0
%              size(idxnow)
%              size(r{i})
            Prow=prefactor*exp( (r{i}.*r{i}) /(-4*DMT*dt/beta) );
            rowsum=sum(Prow);
            Prow(1)=Prow(1)-diag(rowsum-1);
            dM=beta*(Prow*A(idxghost(idxnow),7:end)-A(partnow,7:end));
            masstemp(partnow,:)=massnew(partnow,:)+dM;  % keep new mass on all ghosts
            massnew(partnow,:)=masstemp(partnow,:);       % push only local masses up

          end
      end
  end

%      P=squareform(pdist(X(:,1:ndim)));
%      P=prefactor*exp((P.*P)/(-4*DMT*dt/beta));
      %colsum=sum(P); 
%      rowsum=sum(P,2);
%      min(rowsum)
%      max(rowsum)
%      P=P-diag(rowsum-1);
%      sfactor(sfactor<1e-20)=0;  sfactor(sfactor==0)=1;
%      P=P./((colsum+rowsum)/2);
%      rowsum=sum(P,2);
%      P=beta*(diag(rowsum)-P);      
%      P=beta*(diag(P*ones(npad,1))-P);
%      dM=beta*(P*X(:,7:end)-X(:,7:end));
%      masstemp(idxghost,:)=massnew(idxghost,:)+dM;  % keep new mass on all ghosts
%      massnew(idxlocal,:)=masstemp(idxlocal,:);       % push only local masses up
      %A(idxlocal,7:end)=masstemp(idxlocal,:);
%  end  % more than 1 particle if statement

  end  % cell loop
  A(:,7:end)=massnew;
  if(kappa<1)
      if gpu>0
         rvec=randn(np,ndim,'gpuArray');
      else
         rvec=randn(np,ndim);
      end
      A(:,1:ndim)=A(:,1:ndim)+sqrt(2*DRW*dt)*rvec;
      A(A(:,1)<xmin(1),1)=-A(A(:,1)<xmin(1),1);
      A(A(:,2)<xmin(2),2)=-A(A(:,2)<xmin(2),2);
      A(A(:,1)>xmax(1),1)=2*xmax(1)-A(A(:,1)>xmax(1),1);
      A(A(:,2)>xmax(2),2)=2*xmax(2)-A(A(:,2)>xmax(2),2);

  end
end  % Timesteps loop
elapsedtime=toc
speed(nn,1)=constnx;  speed(nn,2)=elapsedtime;
save('speedout_noGPU_dt10_1e5.dat','speed','-ascii');


  %  Compare to analytic soln at particle points
sep=pdist2(A(1:nsource,1:ndim),A(:,1:ndim));
analytic=(4*pi*Dtot*ttot)^(-ndim/2)*exp(sep.*sep/(-4*Dtot*ttot));
analytic=sum(analytic,1);
xplot=linspace(xmin(1),xmax(1),100);
analplot=(4*pi*Dtot*ttot)^(-ndim/2)*exp((xplot-A(1,1)).^2/(-4*Dtot*ttot));

xline=0.5*(xmax(1)-xmin(1)); tol=0.01*(xmax(1)-xmin(1)); 
idxplot = abs(A(:,1)-xline)<tol ;
figure(1)
plot(A(idxplot,2),A(idxplot,7),'o')
hold on
plot(A(idxplot,2),analytic(idxplot),'+')
plot(xplot,analplot)
hold off
figure(2)
semilogy(A(idxplot,2),A(idxplot,7),'o')
hold on
axis([0 10 1e-10 1])
plot(A(idxplot,2),analytic(idxplot),'+')
plot(xplot,analplot)
hold off
RMSE=sqrt(mean( ( A(:,7)-analytic' ).^2))



end   % End the various discretizations

