%poisson2D
function lhs=poissonNH3D(rhs,epsln,h,opts)
%set up L
n=size(rhs);

Lx = genLapl1D(n(1),h(1));
Ly = genLapl1D(n(2),h(2));
Lz = genLapl1D(n(3),h(3));
L = kron(speye(n(2),n(2)),Lx) + kron(Ly,speye(n(1),n(1)));
L = kron(Lz,speye(n(1)*n(2),n(1)*n(2))) + kron(speye(n(3),n(3)),L);

L=epsln(:).*L;

Dx = genDiffMat(n(1),h(1));
Dy = genDiffMat(n(2),h(2));
Dz = genDiffMat(n(3),h(3));

Dx = kron(speye(n(3)*n(2)),Dx);
Dy = kron(speye(n(3)),kron(Dy,speye(n(1))));
Dz = kron(Dz,speye(n(2)*n(1)));

edgeInd=[];

tic
for kk=1:n(3) %z ind
    for jj=1:n(2) %y ind
        for ii=1:n(1) % x ind
            ind=(kk-1)*n(1)*n(2)+(jj-1)*n(1)+ii;
            %apply zero on the faces, disregard edges, apply dirichlet on
            %faces
            if (ii==1 || jj==1 || kk==1 || ii==n(1) || jj==n(2)|| kk==n(3))
                L(ind,:)=0;
                Dx(ind,:)=0;
                Dy(ind,:)=0;
                Dz(ind,:)=0;
                edgex=(jj==1 || jj==n(2)) && (kk==1 || kk==n(3));
                edgey=(kk==1 || kk==n(3)) && (ii==1 || ii==n(1));
                edgez=(ii==1 || ii==n(1)) && (jj==1 || jj==n(2));
                if (edgex || edgey || edgez)
                    edgeInd=[edgeInd; ind];
                    continue
                else
                    L(ind,ind)=1;
                end
            end
        end
    end
end
toc

tic
%add non-homogenous part
L=L+spdiags(Dx*epsln(:),0,numel(epsln),numel(epsln))*Dx+...
    spdiags(Dy*epsln(:),0,numel(epsln),numel(epsln))*Dy+...
    spdiags(Dz*epsln(:),0,numel(epsln),numel(epsln))*Dz;

%get rid of corner values in the laplacian and reshape rhs
L(edgeInd,:)=[];
L(:,edgeInd)=[];
rhs=rhs(:);
rhs(edgeInd)=[];
toc

%solve
tic
switch opts.method
    case '\'
        lhs=L\rhs;
    case 'cong'
      afun = @(t) L*t;
      lhs = bicgstab(afun,rhs,1E-6,500);
end
toc
%put the corner values back into the lhs and rhs and reshape to mesh
for ii=edgeInd'
    lhs=[lhs(1:ii-1);inf;lhs(ii:end)];
    rhs=[rhs(1:ii-1);inf;rhs(ii:end)];
end
lhs=reshape(lhs,[n(1) n(2) n(3)]);
rhs=reshape(rhs,[n(1) n(2) n(3)]);
end

function Lx = genLapl1D(nx,dx)
Lx = repmat([1 -2 1]/dx^2,nx,1);
Lx = spdiags(Lx,-1:1,nx,nx);
end

function Dx = genDiffMat(nx,dx)
Dx = repmat([-1 0 1]/(2*dx),nx,1);
Dx = spdiags(Dx,-1:1,nx,nx);
end