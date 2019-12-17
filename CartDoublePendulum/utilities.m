classdef utilities
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
   
    methods (Static)
        function trajectoryAnimator(model, seq, ax, savename)
            
            % Create a handle to the axis, if there isn't one already
            if nargin < 3 || isempty(ax)
               figure();
               ax = gca;
            end
            % Create a video
            if nargin >= 4
                w = VideoWriter(savename);
                open(w);
            end
            
            % Draw all the components of the trajectory first
            T = size(seq, 2);
            data = cell(1,T);
            xlims = zeros(1,2);
            ylims = zeros(1,2);
            for n = 1:T
               [x,y] = model.draw(seq(:,n)); 
                data{n} = [x;y];
                xlims(1) = min([xlims(1), x]);
                xlims(2) = max([xlims(2),x]);
                ylims(1) = min([ylims(1),y]);
                ylims(2) = max([ylims(2),y]);
            end
            % Extend the limits by 5% 
            xlims(1) = xlims(1) * (1 - sign(xlims(1)) * 0.05);
            xlims(2) = xlims(2) * (1 + sign(xlims(2)) * 0.05);
            ylims(1) = ylims(1) * (1 - sign(ylims(1)) * 0.05);
            ylims(2) = ylims(2) * (1 + sign(ylims(2)) * 0.05);
            % Make the limits square
            xr = range(xlims);
            yr = range(ylims);
            if xr > yr
               ylims = ylims * (xr / yr); 
            else
               xlims = xlims * (yr/xr);
            end
            
          
            % Now draw the terrain
            [xterrain, yterrain] = model.terrain.draw(xlims, T);
            plot(ax, xterrain, yterrain, 'k-', 'LineWidth',2);
            hold on;
            % Draw the last point in the trajectory
            plot(ax,data{end}(1,:),data{end}(2,:),'r--','LineWidth',2);
            
            % Draw the first point on the trajectory
            sketch = plot(ax, data{1}(1,:), data{1}(2,:),'b-','LineWidth',2);
            % Set the axis limits
            set(ax,'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
            
            % Now draw all the frames in sequence and create a video
            for n = 1:T
               set(sketch, 'XData',data{n}(1,:),'YData',data{n}(2,:));
               drawnow;
               if nargin == 4
                  writeVideo(w, getframe(ax));
               else
                   pause(0.1);
               end
            end
            % Close the video to save it.
            if nargin == 4
               close(w); 
            end
            
            
        end
                
        function animator(ax, draw, seq, savename)
            %UNTITLED Summary of this function goes here
            %   Detailed explanation goes here
            
            [~,N] = size(seq);
            
            if nargin >= 4
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

