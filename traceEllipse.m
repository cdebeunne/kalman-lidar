function hE = traceEllipse(X0,Cx,couleur, hE)
% function traceEllipse(X0,Cx,couleur)
% Display an ellipse with 1 sigma probability centered in X0 with covariance Cx
% Input:
% - X0: 2x1 : center of the ellipse
% - Cx: 2x2 : covariance matrix
% - couleur : {'b','r'...} or 1x3 vector RGB
% example : traceEllipse([1;2],[4,2;2,5],'r');

if (Cx(1,1) ~= 0 && Cx(2,2) ~= 0)
    M=10;
    % valeurs propres donnent le rayon au carré
    [V,D] = eig(Cx);
    
    L2 = abs(D(1,1));
    l2 = D(2,2);
    L = sqrt(L2);
    
    if L~=0
        k=1;
        for x=-L:L/M:L
            y = real(sqrt(l2*(1-x^2/L2)));
            Y1(:,k) = [x;y];
            Y2(:,k) = [x;-y];
            k=k+1;
        end;
        
        
        Y1 = V*Y1 + X0*ones(1,2*M+1);
        hold on;
        Y2 = V*Y2 + X0*ones(1,2*M+1);
        
        if isempty(hE)
            hE.e1 = plot(Y1(1,:),Y1(2,:),couleur);            
            hE.e2 = plot(Y2(1,:),Y2(2,:),couleur);           
        else
            set(hE.e1, 'xdata',Y1(1,:) ,'ydata', Y1(2,:))            
            set(hE.e2, 'xdata',Y2(1,:) ,'ydata', Y2(2,:))
        end
        
    end
    
end

