classdef PlaneObstacle < Obstacle
    
    methods (Access = public)
        
        function this = PlaneObstacle(varargin)
            
            % parse arguments
            parser = inputParser;

            parser.KeepUnmatched = false;
            parser.PartialMatching = false;
            parser.CaseSensitive = false;

            parser.addParameter('point', [0; 0]); % a point on the plane
            parser.addParameter('normal_angle', 0.0); % angle of the plane's normal
            parser.addOptional('d0', 0.02);
            parser.addOptional('color', {});

            parser.parse(varargin{:});
            args = parser.Results;
            
            this@Obstacle(args.d0, args.color);
            
            this.p0 = args.point;
            this.calcPlaneNormal(args.normal_angle);
            
        end

    end
    
    methods (Access = public)
        
        function plot(this, varargin)
        
            % initialize parser with the names and default values of the input arguments
            inPars = inputParser;

            inPars.KeepUnmatched = true;
            inPars.PartialMatching = false;
            inPars.CaseSensitive = false;

            inPars.addRequired('axes');
            inPars.addOptional('Color', this.color);
            inPars.addOptional('LineWidth', 2);
            inPars.addOptional('LineStyle', '--');

            % Parse input arguments
            inPars.parse(varargin{:});
            args = inPars.Results;
            
            hold(args.axes, 'on');

            ax = args.axes;
            points = this.generatePoints(this.n, this.p0, ax.XLim, ax.YLim);
            points2 = this.generatePoints(this.n, this.p0+this.n*this.d0, ax.XLim, ax.YLim);

            pl_args = {'Parent',ax, 'LineStyle',args.LineStyle, 'HandleVisibility','off'};
            pl = plot(points(1,:), points(2,:), pl_args{:}, 'LineWidth',args.LineWidth);
            if ~isempty(args.Color)
                pl.Color = args.Color;
            end
            
            p_m = [mean(points(1,:)); mean(points(2,:))];
            len = 0.05;
            quiver(p_m(1), p_m(2), len*this.n(1), len*this.n(2), 'LineWidth',args.LineWidth, ...
                'Parent',ax, 'HandleVisibility','off', 'AutoScale','off', 'Color',pl.Color);
            
            plot(points2(1,:), points2(2,:), pl_args{:}, 'LineWidth',0.8*args.LineWidth, 'color',[pl.Color, 0.4]);

        end

    end
    
    methods (Access = protected)
        
        function psi = surface_fun(this, p)
           
            psi = dot(this.n, p-this.p0);
            
        end
        
        function psi_grad = surface_fun_grad(this, p)
            
            psi_grad = this.n;
            
        end
        
        function calcPlaneNormal(this, angle)

            this.n = [cosd(angle); sind(angle)];
            
%             if (angle == 0)
%                 this.n = [1; 0];
%             elseif (isinf(angle))
%                 this.n = [0; 1];
%             else
%                 this.n = [cosd(angle); sind(angle)];
%             end
            
%             if (dot(this.p0, this.n) > 0)
%                 this.n = -this.n;
%             end

        end
        
        function points = generatePoints(this, n, p0, x_lim, y_lim)

            x1 = x_lim(1);
            x2 = x_lim(2);
            y1 = y_lim(1);
            y2 = y_lim(2);

            if (n(2) == 0)
                x1 = p0(1);
                x2 = x1;  
            elseif (n(1) == 0)
                y1 = p0(2);
                y2 = y1;
            else 
                c = dot(n, p0);

                y1b = (-n(1)*x1 + c) / n(2);
                y2b = (-n(1)*x2 + c) / n(2);

                y1 = max([y1 y1b y2b]);
                y2 = min([y2 y1b y2b]);

                if (y1 < y_lim(1)), y1 = y_lim(1); end
                if (y1 > y_lim(2)), y1 = y_lim(2); end
                if (y2 < y_lim(1)), y2 = y_lim(1); end
                if (y2 > y_lim(2)), y2 = y_lim(2); end

                x1 = (-n(2)*y1 + c) / n(1);
                x2 = (-n(2)*y2 + c) / n(1);

            end

            points = [x1 x2; y1 y2];

        end

    end
    
    properties (Access = protected)
       
        n % plane normal
        p0 % a point on the plane
        
    end
    
end