classdef Q9
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
        function obj = Q9(ipcount1D, rect, zeroWeight)

            [x1D, w1D] = obj.getIpscheme(ipcount1D, zeroWeight);
            
            if (zeroWeight)
                ipcount1D = ipcount1D+2;
            end
            obj.ipcount1D = ipcount1D;
            obj.ipcount = obj.ipcount1D*obj.ipcount1D;
            obj.rectangular = rect;
            
            k = 0;
            for j=1:obj.ipcount1D
                for i=1:obj.ipcount1D
                    k=k+1;
                    xbase(k) = x1D(i);
                    ybase(k) = x1D(j);
                    wbase(k) = w1D(i)*w1D(j);
                end
            end
            
			syms x y 
            for i=1:length(xbase)
                xb = xbase(i); yb = ybase(i);
				
                x1 = -1; x2 = 0; x3 = 1;
                Nx = [(x-x2)*(x-x3)/((x1-x2)*(x1-x3)); 
                      (x-x1)*(x-x3)/((x2-x1)*(x2-x3)); 
                      (x-x1)*(x-x2)/((x3-x1)*(x3-x2))];
                Ny = [(y-x2)*(y-x3)/((x1-x2)*(x1-x3)); 
                      (y-x1)*(y-x3)/((x2-x1)*(x2-x3)); 
                      (y-x1)*(y-x2)/((x3-x1)*(x3-x2))];
				N = kron(Ny, Nx);
                
                Nbase(i,:) = double(subs(subs(N, xb), yb));
				Gbase(i,:,1) = double(subs(subs(diff(N,x), xb), yb));
                Gbase(i,:,2) = double(subs(subs(diff(N,y), xb), yb));
	
				G2base(i,:,1) = double(subs(subs(diff(diff(N,x),x), xb), yb));
                G2base(i,:,2) = double(subs(subs(diff(diff(N,y),y), xb), yb));		
				G2base(i,:,3) = double(subs(subs(diff(diff(N,x),y), xb), yb));
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
            
            if obj.rectangular 
            	G(:,:,1) = obj.Gbase(:,:,1)*2/(X(3)-X(1));
            	G(:,:,2) = obj.Gbase(:,:,2)*2/(Y(7)-Y(1));
                w = obj.wbase * (X(3)-X(1))/2 * (Y(7)-Y(1))/2;
            else %% probably requires correcting
                dXdXi = obj.Gbase(:,:,1)*X; dXdEta = obj.Gbase(:,:,2)*X;
                dYdXi = obj.Gbase(:,:,1)*Y; dYdEta = obj.Gbase(:,:,2)*Y;

                J(:,1,1) = dXdXi; J(:,1,2) = dXdEta;
                J(:,2,1) = dYdXi; J(:,2,2) = dYdEta;

                G = 0.0*obj.Gbase;
                for i=1:size(obj.Nbase, 1)
                    Jinv = inv(squeeze(J(i,:,:)));
                    for j=1:size(obj.Nbase, 2)
                        G(i,j,:) = Jinv*squeeze(obj.Gbase(i,j,:));
                    end

                    w = obj.wbase(i)*det(squeeze(J));
                end
            
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
            
            if obj.rectangular 
                xy(1,:) = X(1) + (0.5*(obj.xbase+1))*(X(3)-X(1));
                xy(2,:) = Y(1) + (0.5*(obj.ybase+1))*(Y(7)-Y(1));
            else
                %% not implemented
            end

        end
 
    end
    
    
    methods (Access = private)
        function [x1D, w1D] = getIpscheme(obj, ipcount1D, zeroWeight)
            if (ipcount1D == 1)
                x1D = 0;
                w1D = 2;
            elseif (ipcount1D == 2)
                x1D = [-1/sqrt(3); 1/sqrt(3)];
                w1D = [1; 1];                
            elseif (ipcount1D == 3)
                x1D = [-sqrt(3/5); 0; sqrt(3/5)];
                w1D = [5/9; 8/9; 5/9];                     
            elseif (ipcount1D ==4)
                x1D = [-sqrt(3/7+2/7*sqrt(6/5)); -sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7+2/7*sqrt(6/5))];
                w1D = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];                     
            elseif (ipcount1D == 5)
                x1D = [-1/3*sqrt(5+2*sqrt(10/7)); -1/3*sqrt(5-2*sqrt(10/7)); 0; 1/3*sqrt(5-2*sqrt(10/7)); 1/3*sqrt(5+2*sqrt(10/7))];
                w1D = [(322-13*sqrt(70))/900; (322+13*sqrt(70))/900; 128/225; (322+13*sqrt(70))/900; (322-13*sqrt(70))/900];                     
            else
                error("Higer order ip schemes not implemented in Shapes.Q9");
            end
            
            if (zeroWeight)
               x1D = [-1; x1D; 1];
               w1D = [0; w1D; 0];
            end
            
        end
    end
end

