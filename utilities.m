classdef utilities
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
   
    methods (Static)
        function animator(ax, draw, seq, savename)
            %UNTITLED Summary of this function goes here
            %   Detailed explanation goes here
            
            [~,N] = size(seq);
            
            if nargin == 4
                w = VideoWriter(savename);
                open(w);
            end
            
            for n = 1:N
                draw(ax,seq(:,n));
                if nargin == 4
                    writeVideo(w, getframe(ax));
                else
                    pause(0.1);
                end
            end
            if nargin == 4
                close(w);
            end
        end
        
        function diagnostics(model, t, q, f, name)
            % Calculate the distances
            d = zeros(size(f,1),size(q,2));
            for k = 1:size(q,2)
               [n,a] = model.contactNormal(q(:,k));
               d(:,k) = n*q(:,k) - a;
            end
            % Make figures and plot the distances and normal forces together
            for k = 1:size(f,1)
               figure('Name',name);
               yyaxis left;
               plot(t,f(k,:),'-');
               ylabel('Normal Impulse (Ns)');
               ylim([-1,1]*max(abs(f(k,:))));
               yyaxis right;
               plot(t,d(k,:),'-');
               ylim([-1,1]*max(abs(d(k,:))));
               ylabel('Distance');
               xlabel('Time (s)');                
            end
            
        end
        
        function errorPlots(t,f,r,name)
           figure('Name',name); 
           mf = sqrt(sum(f.^2,1));
           yyaxis left;
           plot(t,mf,'-');
           ylabel('Force Magnitude');
           ylim([-1,1]*max(mf));
           yyaxis right;
           plot(t,r,'-');
           ylim([-1,1]*max(abs(r)));
           ylabel('Complementarity Residual');
           xlabel('Time (s)');
           
            
        end
    end
end

