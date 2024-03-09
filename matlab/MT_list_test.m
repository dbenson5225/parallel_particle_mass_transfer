%  starting from scratch to verify some things ...
% This code is more general, but is set up to recreate the results in the
% Appendix, where we twst the regularly-spaced list-based algorithm for
% accuracy.  Vary number of particles np
clear all;

%gpu=gpuDeviceCount;      %  See if there's a gpu
gpu=0;                   %  Turn off gpu
%define numerical and physical parameters
np=101*101; ens=1;  % # A points, #particles per point, # runs, point radius
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
Sx=[1]
%Sx=[ 37 31 25 19 13 9 5 3 2 ]  %  # subdomains in one direction
%Sx=31
for nn=1:length(Sx);       % Go through all of the discretizations
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
%locations of A and B ICs;
x_0=0.5*(xmax-xmin);
concfactor=prod(xmax-xmin)/np;
volume=prod(xmax(1:ndim)-xmin(1:ndim));

%Uniform random locations %%%%%%%%%%%%%%%%%%%%%%%
% for i=1:ndim
% A(1:np,i)=(xmax(i)-xmin(i))*rand(np,1)+xmin(i);
% end
% 
% A(1:nsource,6+1)=1/nsource; A(1:nsource,6+2)=1/nsource/concfactor; % random x,y
% A(1:nsource,1:ndim)=x_0; A(1:nsource,7:8)=1/nsource/concfactor;   % Fixed x,y

dx=ones(1,3); dx(1)=(volume/np)^(1/ndim); dx(2)=(volume/np)^(1/ndim);

% % % Regularly spaced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 dxtrial=(volume/np)^(1/ndim);
 nx=ones(1,3)
 for i=1:ndim
   nx(i) = floor((xmax(i)-xmin(i))/dxtrial);
   dx(i) = (xmax(i)-xmin(i))/nx(i);
 end
 np=prod(nx)
 count=1;
 for j=1:nx(2)
     for i=1:nx(1)
         A(count,1)=dx(1)*(i-0.5);
         A(count,2)=dx(2)*(j-0.5);
         count=count+1;
     end
 end
 A(count:end,:)=[];
 [val,idx1]=min(abs(A(:,1)-x_0(1)) + abs(A(:,2)-x_0(2))); 
 A(idx1,7:8)=1/concfactor;
 x_0(1)=A(idx1,1); x_0(2)=A(idx1,2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start a time loop  %%%%%%%%%%%%%%

tic
prefactor=ds/((4*pi*DMT*dt/beta)^(ndim/2));
for n=1:ntsteps

    sigma = sqrt(4*DMT*dt);  % Char. diffusion dist.
    % Do the mass-transfer
    % First chop up space: ideally each chunk is 3*sigma
    nxmt=ones(1,3);  if(gpu>0); nxmt=gpuArray(nxmt); end;
  for i=1:ndim
    pad(i)=4*sigma; if(gpu>0); pad=gpuArray(pad); end;

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
%  masstemp=massnew;     % This holds temporary masses (incl. ghosts)
  lidx=1:prod(nxmt);   
  [is js ks]=ind2sub(sz,lidx);

  for m=1:prod(nxmt)            % go to every cell  (cell loop)
      ijk(1)=is(m); ijk(2)=js(m); ijk(3)=ks(m);   
      lims=zeros(3,2);  if(gpu>0); lims=gpuArray(lims); end;
      for i=1:3
         lims(i,1)=xlo(1,ijk(i))-pad(i); lims(i,2)=xhi(2,ijk(i))+pad(i);  %xlo and xhi
      end
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
      [idx,r]=rangesearch(A(idxghost,1:ndim),A(idxlocal,1:ndim),6*sigma);
      for i=1:length(idxlocal)
          partnow=idxlocal(i);
          idxnow=idx{i};
          numnear=length(idxnow);
          if numnear>0
            Prow=prefactor*exp( (r{i}.*r{i}) /(-4*DMT*dt/beta) );
%            if(abs(Prow(1)-prefactor)>1e-10)
%                Prow(1)
%            end
            rowsum=sum(Prow);
            Prow(1)=Prow(1)-diag(rowsum-1);
            dM=beta*(Prow*A(idxghost(idxnow),7:end)-A(partnow,7:end));
            massnew(partnow,:)=massnew(partnow,:)+dM;  % keep new mass on all ghosts

          end
      end
  end

  end  % cell loop
  A(:,7:end)=massnew;

   %  Random walks next ...
  if(kappa<1)
      if gpu>0
         rvec=randn(np,ndim,'gpuArray');
      else
         rvec=randn(np,ndim);
      end
      A(:,1:ndim)=A(:,1:ndim)+sqrt(2*DRW*dt)*rvec;
      for i=1:ndim
        A(A(:,i)<xmin(i),i)=-A(A(:,i)<xmin(i),i);
        A(A(:,i)>xmax(i),i)=2*xmax(i)-A(A(:,i)>xmax(i),i);
      end
  end

end  % Timesteps loop
elapsedtime=toc
speed(nn,1)=constnx;  speed(nn,2)=elapsedtime;
%save('speedout_noGPU_dt10_1e5.dat','speed','-ascii');


  %  Compare to analytic (Gaussian) soln at particle points
IC=x_0;
sep2=(A(:,1:ndim) - IC);
sep2=sum(sep2.*sep2,2);
analytic=(4*pi*Dtot*ttot)^(-ndim/2)*exp(sep2/(-4*Dtot*ttot));
xplot=linspace(xmin(1)+dx(1)/2,xmax(1)-dx(1)/2,100);
analplot=(4*pi*Dtot*ttot)^(-ndim/2)*exp((xplot-x_0(1)).^2/(-4*Dtot*ttot));

xline=IC(1);
tol=dx(2)/2; 
idxplot = abs(A(:,2)-xline)<tol ;
figure(1)
plot(A(idxplot,1),A(idxplot,7),'o')
hold on
plot(A(idxplot,1),analytic(idxplot),'+')
plot(xplot,analplot)
legend('Particles','Analytic at points','Analytic through x_0')
hold off
figure(2)
semilogy(A(idxplot,1),A(idxplot,7),'o')
hold on
axis([0 10 1e-10 1])
plot(A(idxplot,1),analytic(idxplot),'+')
plot(xplot,analplot)
legend('Particles','Analytic at points','Analytic through x_0')
hold off

RMSE=sqrt(mean( ( A(:,7)-analytic ).^2))


end   % End the various discretizations

