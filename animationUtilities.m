classdef animationUtilities
    
    methods (Static)
        function multiTrajectoryAnimator(model, seq, target, savename, labels, format)
            
            % Strip the extension from the savename
            C = strsplit(savename,'.');
            savename=C{1};
            if nargin < 6
                format = 'Motion JPEG AVI';
                ext = '.avi';
            else
                switch format
                    case 'MPEG-4'
                        ext = '.mp4';
                    case 'Motion JPEG 2000'
                        ext = '.mj2';
                    otherwise
                        ext = '.avi';
                end
            end
            % Create a handle to the axis
            figure();
            ax = gca;
            % Pause to allow the figure to appear
            pause(1);
            % Check for multiple sequences
            if ~iscell(seq)
                seq = {seq};
            end
            numTraj = length(seq);
            if ~iscell(labels)
                labels = {labels};
            end
            % Check for labels
            if nargin >=5 && ~isempty(labels)
                if length(labels) ~= numTraj
                    error('There must be an equal number of trajectories and labels');
                end
            else
                labels = cell(1,numTraj);
                for k = 1:numTraj
                    labels{k} = sprintf('Trajectory %d',k);
                end
            end
            
            N = size(seq{1},2);
            for k = 2:numTraj
                if N ~= size(seq{k},2)
                    error('All sequences must have the same number of columns');
                end
            end
            % Create a video
            if nargin >= 4 && ~isempty(savename)
                w = VideoWriter([savename, ext], format);
                open(w);
            end
            xdata = cell(numTraj,N);
            ydata = cell(numTraj,N);
            xlims = zeros(1,2);
            ylims = zeros(1,2);
            
            % Draw the target position, if given
            if nargin > 2 && ~isempty(target)
                [x, y] = model.draw(target);
                % Update the target position
                xlims(1) = min([xlims(1), x]);
                xlims(2) = max([xlims(2), x]);
                ylims(1) = min([ylims(1), y]);
                ylims(2) = max([ylims(2), y]);
                % Plot the target
                plot(ax,x,y,'r--','LineWidth',2,'DisplayName','target');
                hold on;
            end
            
            % Draw all the frames first
            for k = 1:numTraj
                % Draw all the components of all trajectories first
                for n = 1:N
                    [x,y] = model.draw(seq{k}(:,n));
                    xdata{k,n} = x;
                    ydata{k,n} = y;
                    %data{k, n} = [x;y];
                    xlims(1) = min([xlims(1), x]);
                    xlims(2) = max([xlims(2), x]);
                    ylims(1) = min([ylims(1), y]);
                    ylims(2) = max([ylims(2), y]);
                end
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
            [xterrain, yterrain] = model.terrain.draw(xlims, N);
            plot(ax, xterrain, yterrain, 'k-', 'LineWidth',2,'DisplayName','terrain');
            hold on;
            % Draw the first point on each of the trajectories
            colors = lines(numTraj);
            sketches = cell(1,numTraj);
            for k = 1:numTraj
                sketches{k} = plot(ax, xdata{k,1}, ydata{k,1},'LineStyle','-','Color',colors(k,:),'LineWidth',2,'DisplayName',labels{k});
            end
            hold off;
            % Set the axis limits
            set(ax,'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
            legend show;
            legend boxoff;
            % Now draw all the frames in sequence and create a video
            for n = 1:N
                set([sketches{:}], {'XData'},xdata(:,n), {'YData'},ydata(:,n));
                set(ax,'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
                drawnow;
                if nargin >= 4 && ~isempty(savename)
                    writeVideo(w, getframe(ax));
                else
                    pause(0.1);
                end
            end
            % Close the video to save it.
            if nargin >= 4 && ~isempty(savename)
                close(w);
            end
        end  
        function plotTrajectoryFrames(model, seq, nframes, target, labels)
            
            
            % Create a figure
            figure();
            ax = gca;
            
            % Check the number of trajectories
            if ~iscell(seq)
                seq = {seq};
            end
            % Total number of trajectories
            numTraj = length(seq);
            % Check for labels
            if nargin >=5 && ~isempty(labels)
                if length(labels) ~= numTraj
                    error('There must be an equal number of trajectories and labels');
                end
            else
                labels = cell(1,numTraj);
                for k = 1:numTraj
                    labels{k} = sprintf('Trajectory %d',k);
                end
            end
            % Check the size of each trajectory
            N = size(seq{1},2);
            for k = 2:numTraj
                if N ~= size(seq{k},2)
                    error('All sequences must have the same number of columns');
                end
            end
            % Subsample the number of frames to plot
            idx = round(linspace(1,N,nframes));
            
            % Get the data to plot
            xdata = cell(numTraj,nframes);
            ydata = cell(numTraj,nframes);
            xlims = zeros(1,2);
            ylims = zeros(1,2);
            for k = 1:numTraj
                % Draw all the components of all trajectories first
                for n = 1:nframes
                    [x,y] = model.draw(seq{k}(:,idx(n)));
                    xdata{k,n} = [x(:);NaN];
                    ydata{k,n} = [y(:);NaN];
                    xlims(1) = min([xlims(1), x]);
                    xlims(2) = max([xlims(2), x]);
                    ylims(1) = min([ylims(1), y]);
                    ylims(2) = max([ylims(2), y]);
                end
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
            [xterrain, yterrain] = model.terrain.draw(xlims, N);
            plot(ax, xterrain, yterrain, 'k-', 'LineWidth',2,'DisplayName','terrain');
            hold on;
            % Draw the target position, if given
            if nargin > 2 && ~isempty(target)
                [xT, yT] = model.draw(target);
                plot(ax,xT, yT,'r--','LineWidth',2,'DisplayName','target');
            end
            % Set the axis limits
            set(ax,'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
            legend show;
            legend boxoff;
            % Now draw the selected keyframes
            for k = 1:numTraj
                plot(ax, cat(1,xdata{k,:}), cat(1,ydata{k,:}), 'LineWidth',2,'DisplayName',labels{k});
            end
            
        end
        function visualizeTrajectory(plant,q,savename)
            %% VisualizeTrajectory: Visualizes a single trajectory of a dynamical system
            %
            pause(1);
            %% Get the data
            % Gather the frames together and determine the framesize
            N = size(q,2);
            xdata = cell(1,N);
            ydata = cell(1,N);
            xlims = [inf, -inf];
            ylims = [inf, -inf];
            for n = 1:N
                [xdata{n}, ydata{n}] = plant.visualize(q(:,n));
                xlims(1) = min([xlims(1); xdata{n}(:)]);
                xlims(2) = max([xlims(2); xdata{n}(:)]);
                ylims(1) = min([ylims(1); ydata{n}(:)]);
                ylims(2) = max([ylims(2); ydata{n}(:)]);
            end
            % Extend the limits by 5%
            xlims = xlims .* (1 + [-1,1].*sign(xlims).*0.05);
            ylims = ylims .* (1 + [-1, 1].*sign(ylims).*0.05);
            % Make the limits square
            xr = range(xlims);
            yr = range(ylims);
            if xr > yr
                ylims = ylims * (xr/yr);
            else
                xlims = xlims * (yr/xr);
            end
            % Get the terrain
            [xterrain, yterrain] = plant.terrain.draw(xlims,N);
            % Extend the terrain so we can draw it as a patch
            xterrain = [xterrain(:); xterrain(end); xterrain(1); xterrain(1)];
            yterrain = [yterrain(:); ylims(1); ylims(1); yterrain(1)];
            
            %% Create the animation
            % Get a figure
            figure('visible','on');
            ax = gca;
            if nargin > 2 && ~isempty(savename)
                writer = VideoWriter(savename);
                writer.open();
            end
            
            % Draw the terrain and the first frame
            patch(ax, xterrain, yterrain, [0.4,0.4,0.4]);
            hold on;
            NPatch = size(xdata{1}, 2);
            patches = cell(1,NPatch);
            colors = lines(NPatch);
            for k = 1:NPatch
                patches{k} = patch(ax, xdata{1}(:,k), ydata{1}(:,k), colors(k,:));
            end
            set(ax, 'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
            % Get the frame
            drawnow limitrate;
            hold on;
            if nargin > 2 && ~isempty(savename)
                pause(0.1);
                writer.writeVideo(getframe(ax));
            end
            % Draw all the frames and create a video
            for n = 2:N
                for k = 1:NPatch
                    set(patches{k}, 'xdata', xdata{n}(:,k), 'ydata', ydata{n}(:,k));
                end
                set(ax, 'XLimMode','manual','YLimMode','manual');
                drawnow limitrate;
                if nargin > 2 && ~isempty(savename)
                    pause(0.1);
                    writer.writeVideo(getframe(ax));
                end
            end
            
            % Close the video to save it
            if nargin > 2 && ~isempty(savename)
                writer.close();
            end
            
        end
        function visualizeMultiTrajectory(model, seq, savename, labels)
            %% VisualizeMultiTrajectory

            % Create a handle to the axis
            figure();
            ax = gca;
            % Pause to allow the figure to appear
            pause(1);
            % Check for multiple sequences
            if ~iscell(seq)
                seq = {seq};
            end
            numTraj = length(seq);
            if ~iscell(labels)
                labels = {labels};
            end
            % Check for labels
            if nargin >=4  && ~isempty(labels)
                if length(labels) ~= numTraj
                    error('There must be an equal number of trajectories and labels');
                end
            else
                labels = cell(1,numTraj);
                for k = 1:numTraj
                    labels{k} = sprintf('Trajectory %d',k);
                end
            end
            
            N = size(seq{1},2);
            for k = 2:numTraj
                if N ~= size(seq{k},2)
                    error('All sequences must have the same number of columns');
                end
            end
            % Create a video
            if nargin >3 && ~isempty(savename)
                writer = VideoWriter(savename);
                writer.open();
            end
            xdata = cell(numTraj,N);
            ydata = cell(numTraj,N);
            xlims = zeros(1,2);
            ylims = zeros(1,2);
            
            
            % Draw all the frames first
            for k = 1:numTraj
                % Draw all the components of all trajectories first
                for n = 1:N
                    [xdata{k,n},ydata{k,n}] = model.visualize(seq{k}(:,n));
                    xlims(1) = min([xlims(1); xdata{k,n}(:)]);
                    xlims(2) = max([xlims(2); xdata{k,n}(:)]);
                    ylims(1) = min([ylims(1); ydata{k,n}(:)]);
                    ylims(2) = max([ylims(2); ydata{k,n}(:)]);
                end
            end
            
            % Extend the limits by 5%
            xlims = xlims .* (1 + [-1,1].*sign(xlims).*0.05);
            ylims = ylims .* (1 + [-1, 1].*sign(ylims).*0.05);
            % Make the limits square
            xr = range(xlims);
            yr = range(ylims);
            if xr > yr
                ylims = ylims * (xr / yr);
            else
                xlims = xlims * (yr/xr);
            end
            
            % Now draw the terrain
            [xterrain, yterrain] = model.terrain.draw(xlims, N);
            % Extend the terrain so we can draw it as a patch
            xterrain = [xterrain(:); xterrain(end); xterrain(1); xterrain(1)];
            yterrain = [yterrain(:); ylims(1); ylims(1); yterrain(1)];
            % Draw the terrain and the first frame
            patch(ax, xterrain, yterrain, [0.4,0.4,0.4],'DisplayName','terrain');
            hold on;
            % Draw the first point on each of the trajectories
            sketches = cell(1,numTraj);
            colors = lines(numTraj);
            for k = 1:numTraj
                sketches{k} = patch(ax, xdata{k,1},ydata{k,1},colors(k,:),'LineWidth',1.5,'DisplayName',labels{k});
            end
            % Set the axis limits
            set(ax,'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
            legend show;
            legend boxoff;
            % Now draw all the frames in sequence and create a video
            for n = 1:N
                set([sketches{:}], {'XData'},xdata(:,n), {'YData'},ydata(:,n));
                set(ax,'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
                drawnow;
                if nargin >= 3 && ~isempty(savename)
                    writer.writeVideo(getframe(ax));
                else
                    pause(0.01);
                end
            end
            % Close the video to save it.
            if nargin >= 3 && ~isempty(savename)
               writer.close();
            end
            
        end
        function visualizeFrames(model, seq, nframes, labels)
             % Create a figure
            figure();
            ax = gca;
            
            % Check the number of trajectories
            if ~iscell(seq)
                seq = {seq};
            end
            % Total number of trajectories
            numTraj = length(seq);
            % Check for labels
            if nargin >=4 && ~isempty(labels)
                if length(labels) ~= numTraj
                    error('There must be an equal number of trajectories and labels');
                end
            else
                labels = cell(1,numTraj);
                for k = 1:numTraj
                    labels{k} = sprintf('Trajectory %d',k);
                end
            end
            % Check the size of each trajectory
            N = size(seq{1},2);
            for k = 2:numTraj
                if N ~= size(seq{k},2)
                    error('All sequences must have the same number of columns');
                end
            end
            % Subsample the number of frames to plot
            idx = round(linspace(1,N,nframes));
            
            % Get the data to plot
            xdata = cell(numTraj,nframes);
            ydata = cell(numTraj,nframes);
            xlims = zeros(1,2);
            ylims = zeros(1,2);
            for k = 1:numTraj
                % Draw all the components of all trajectories first
                for n = 1:nframes
                    [xdata{k,n}, ydata{k,n}] = model.visualize(seq{k}(:,idx(n)));
                    xlims(1) = min([xlims(1); xdata{k,n}(:)]);
                    xlims(2) = max([xlims(2); xdata{k,n}(:)]);
                    ylims(1) = min([ylims(1); ydata{k,n}(:)]);
                    ylims(2) = max([ylims(2); ydata{k,n}(:)]);
                end
            end
            
            % Extend the limits by 5%
            xlims = xlims .* (1 + [-1,1].*sign(xlims).*0.05);
            ylims = ylims .* (1 + [-1, 1].*sign(ylims).*0.05);
            % Make the limits square
            xr = range(xlims);
            yr = range(ylims);
            if xr > yr
                ylims = ylims * (xr / yr);
            else
                xlims = xlims * (yr/xr);
            end
            
              % Now draw the terrain
            [xterrain, yterrain] = model.terrain.draw(xlims, N);
            % Extend the terrain so we can draw it as a patch
            xterrain = [xterrain(:); xterrain(end); xterrain(1); xterrain(1)];
            yterrain = [yterrain(:); ylims(1); ylims(1); yterrain(1)];
            % Draw the terrain and the first frame
            patch(ax, xterrain, yterrain, [0.4,0.4,0.4],'DisplayName','terrain');
            hold on;
            % Set the axis limits
            set(ax,'XLimMode','manual','YLimMode','manual','xlim',xlims,'ylim',ylims);
            legend show;
            legend boxoff;
            % Now draw the selected keyframes
            colors = lines(numTraj);
            for k = 1:numTraj
                patch(ax, cat(2,xdata{k,:}), cat(2,ydata{k,:}), colors(k,:), 'LineWidth',1.5,'DisplayName',labels{k});
            end
        end
    end
end