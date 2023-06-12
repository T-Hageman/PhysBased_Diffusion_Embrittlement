classdef T6
    %Q9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ipcount1D
        ipcount
        
        rectangular
        
        Nbase   % ip, shapefunc
        Gbase   % ip, shapefunc, dx/dy
		G2base  % ip, shapefunc, dxx/dyy/dxy
        
        wbase   % ip
        xbase   % ip
        ybase   % ip
    end
    
    methods
        function obj = T6(ipcount1D, rect, zeroWeight)

            [x, w] = obj.getIpscheme(ipcount1D, zeroWeight);
			obj.rectangular = false;
            
            if (zeroWeight)
                ipcount1D = ipcount1D+2;
            end
            obj.ipcount1D = ipcount1D;

			xbase = x(:,1);
			ybase = x(:,2);
			wbase = w;

			obj.ipcount = length(w);
			
			if true
            	syms s t 
            	for i=1:length(xbase)
                	x = xbase(i); y = ybase(i);
                	x1 = 0; x2 = 1/2; x3 = 1;
	
					%Ntriangle6nodes=[2*(1-s-t)*((1/2)-s-t) 2*s*(s-(1/2)) 2*t*(t-(1/2)) 4*s*(1-s-t) 4*s*t 4*t*(1-s-t)];
					Ntriangle6nodes=0.5*[2*(1-s-t)*((1)-s-t) 2*s*(s-(0)) 2*t*(t-(0)) 4*s*(1-s-t) 4*s*t 4*t*(1-s-t)];
	
					Nbase(i,:) = double(subs(subs(Ntriangle6nodes, x), y));
					Gbase(i,:,1) = double(subs(subs(diff(Ntriangle6nodes,s), x), y));
                	Gbase(i,:,2) = double(subs(subs(diff(Ntriangle6nodes,t), x), y));
	
					G2base(i,:,1) = double(subs(subs(diff(diff(Ntriangle6nodes,s),s), x), y));
                	G2base(i,:,2) = double(subs(subs(diff(diff(Ntriangle6nodes,t),t), x), y));		
					G2base(i,:,3) = double(subs(subs(diff(diff(Ntriangle6nodes,s),t), x), y));
				end
				save('./Shapes/@T6/NG.mat','Nbase','Gbase','G2base')
			else
				load('./Shapes/@T6/NG.mat','Nbase','Gbase','G2base')
			end
            
            obj.xbase = xbase;
            obj.ybase = ybase;
            obj.wbase = wbase;
            obj.Nbase = Nbase;
            obj.Gbase = Gbase;
			obj.G2base= G2base;
            
        end
        
        function [N, G, w] = getVals(obj, X, Y)
            N = obj.Nbase;
            
            dXdXi = obj.Gbase(:,:,1)*X; dXdEta = obj.Gbase(:,:,2)*X;
            dYdXi = obj.Gbase(:,:,1)*Y; dYdEta = obj.Gbase(:,:,2)*Y;

            J(:,1,1) = dXdXi; J(:,2,1) = dXdEta;
            J(:,1,2) = dYdXi; J(:,2,2) = dYdEta;

            G = 0.0*obj.Gbase;
            for i=1:size(obj.Nbase, 1)
                Jinv = inv(squeeze(J(i,:,:)));
                for j=1:size(obj.Nbase, 2)
                    G(i,j,:) = Jinv*squeeze(obj.Gbase(i,j,:));
                end

                w(i) = obj.wbase(i)*abs(det(squeeze(J(i,:,:))));
            end
            

		end

		function G2 = getG2(obj, X, Y)
            G2 = 0*obj.G2base;
            
            dXdXi = obj.Gbase(:,:,1)*X; dXdEta = obj.Gbase(:,:,2)*X;
            dYdXi = obj.Gbase(:,:,1)*Y; dYdEta = obj.Gbase(:,:,2)*Y;

			J(:,1,1) = dXdXi; J(:,2,1) = dXdEta;
            J(:,1,2) = dYdXi; J(:,2,2) = dYdEta;

			dXdXiXi   = obj.G2base(:,:,1)*X;
			dXdXiEta  = obj.G2base(:,:,3)*X;
			dXdEtaEta = obj.G2base(:,:,2)*X;
			dYdXiXi   = obj.G2base(:,:,1)*Y;
			dYdXiEta  = obj.G2base(:,:,3)*Y;
			dYdEtaEta = obj.G2base(:,:,2)*Y;

            for i=1:size(obj.Nbase, 1)
				J2(1,1) = dXdXi(i)^2;
				J2(2,1) = dXdEta(i)^2;
				J2(3,1) = dXdXi(i)*dXdEta(i);
				J2(1,2) = dYdXi(i)^2;
				J2(2,2) = dYdEta(i)^2;
				J2(3,2) = dYdXi(i)*dYdEta(i);
				J2(1,3) = 2*dXdXi(i)*dYdXi(i);
				J2(2,3) = 2*dXdEta(i)*dYdEta(i);
				J2(3,3) = dXdEta(i)*dYdXi(i)+dYdEta(i)*dXdXi(i);

				C(1,1) = dXdXiXi(i); C(1,2) = dYdXiXi(i);
				C(2,1) = dXdEtaEta(i); C(2,2) = dYdEtaEta(i);
				C(3,1) = dXdXiEta(i); C(3,2) = dYdXiEta(i);
                for j=1:size(obj.Nbase, 2)
					B(1,1) = obj.G2base(i,j,1);
					B(2,1) = obj.G2base(i,j,2);
					B(3,1) = obj.G2base(i,j,3);

					dXY = inv(squeeze(J(i,:,:)))*squeeze(obj.Gbase(i,j,:));

					B = B - C*dXY;

                    G2(i,j,:) = inv(J2)*B;
				end
            end
            

        end
        
        function xy = getIPGlobal(obj, X,Y)
            
			xy(1,:) = obj.Nbase*X;
			xy(2,:) = obj.Nbase*Y;

%             if obj.rectangular 
%                 xy(1,:) = X(1) + (0.5*(obj.xbase+1))*(X(3)-X(1));
%                 xy(2,:) = Y(1) + (0.5*(obj.ybase+1))*(Y(7)-Y(1));
%             else
%                 %% not implemented
%             end

        end
 
    end
    
    
    methods (Access = private)
        function [x, w] = getIpscheme(obj, ipcount1D, zeroWeight)
            if (ipcount1D == 1)
                x = [1/3, 1/3];
                w = 0.5;
            elseif (ipcount1D == 2)
                x = [2/3, 1/6; 1/6, 1/6; 1/6, 2/3];
				w = [1/3; 1/3; 1/3]*0.5;
            elseif (ipcount1D == 3)
                x = [1/3, 1/3; 3/5, 1/5; 1/5, 1/5; 1/5, 3/5];
                w = [-27/48; 25/48; 25/48; 25/48]*0.5;                     
            else
                error("Higer order ip schemes not implemented in Shapes.Q9");
            end
% 			N = ipcount1D;
% 			v = [0 0; 0 1; 1 0];
% 
% 			n=1:N;  nnk=2*n+1; A=[1/3 repmat(1,1,N)./(nnk.*(nnk+2))];
% 			n=2:N; nnk=nnk(n); B1=2/9; nk=n+1; nnk2=nnk.*nnk;
% 			B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2); ab=[A' [2; B1; B']]; s=sqrt(ab(2:N,2));
% 			[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
% 			[X,I]=sort(diag(X)); x=(X+1)/2; wx=ab(1,2)*V(1,I)'.^2/4;
% 			N=N-1; N1=N+1; N2=N+2;  y=cos((2*(N:-1:0)'+1)*pi/(2*N+2));
% 			L=zeros(N1,N2);  y0=2;  iter=0;
% 			while max(abs(y-y0))>eps    
%     			L(:,1)=1;    L(:,2)=y;   
%     			for k=2:N1
%         			L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
%     			end
%     			Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
%     			y0=y;    y=y0-L(:,N2)./Lp;  iter=iter+1;
% 			end
% 			cd=[ 1, 0, 0; -1, 0, 1; 0, 1,-1]*v; 
% 			t=(1+y)/2;  Wx=abs(det(cd(2:3,:)))*wx;  Wy=1./((1-y.^2).*Lp.^2)*(N2/N1)^2;
% 			[tt,xx]=meshgrid(t,x); yy=tt.*xx;
% 			X=cd(1,1)+cd(2,1)*xx+cd(3,1)*yy;    Y=cd(1,2)+cd(2,2)*xx+cd(3,2)*yy;
% 			W=Wx*Wy';
% 
% 			clear x;
% 			x(:,1)=X(:);
% 			x(:,2)=Y(:);
% 			w = W(:);
            
            if (zeroWeight)
               x = [x; 0, 0; 0, 1; 1, 0];
               w = [w; 0; 0; 0];
            end
            
        end
    end
end

