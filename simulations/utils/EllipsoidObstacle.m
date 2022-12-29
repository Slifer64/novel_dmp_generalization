classdef EllipsoidObstacle < Obstacle
    
    methods (Access = public)
        
        function this = EllipsoidObstacle(varargin)
            
            % parse arguments
            parser = inputParser;

            parser.KeepUnmatched = false;
            parser.PartialMatching = false;
            parser.CaseSensitive = false;

            parser.addParameter('center', [0; 0]); % center of ellipsoid
            parser.addParameter('angle', 0.0); % rotation of ellipsoid around z
            parser.addParameter('lambda_x', 1.0); % scaling in x-axis
            parser.addParameter('lambda_y', 1.0); % scaling in y-axis
            parser.addOptional('d0', 0.5);
            parser.addOptional('color', {});

            parser.parse(varargin{:});
            args = parser.Results;
            
            this@Obstacle(args.d0, args.color);
            
            this.c = args.center;
            this.calcSigma(args.angle, args.lambda_x, args.lambda_y);
            
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
            inPars.addOptional('CenterSize', 10);

            % Parse input arguments
            inPars.parse(varargin{:});
            args = inPars.Results;
            
            hold(args.axes, 'on');

            points = this.generatePoints(this.c, this.Sigma);
            points2 = this.generatePoints(this.c, this.Sigma*(1+this.d0));

            pl_args = {'Parent',args.axes,  'HandleVisibility','off'};
            pl = plot(points(1,:), points(2,:), pl_args{:}, 'LineWidth',args.LineWidth);
            
            if ~isempty(args.Color)
                pl.Color = args.Color;
            end
            
            plot(this.c(1), this.c(2), 'LineWidth',args.LineWidth, 'Marker','x', 'LineStyle','None', 'Markersize',args.CenterSize, 'Color',pl.Color, 'Parent',args.axes, 'HandleVisibility','off');
            plot(points2(1,:), points2(2,:), pl_args{:}, 'LineWidth',0.8*args.LineWidth, 'color',[pl.Color 0.4]);

        end

    end
    
    methods (Access = protected)
        
        function psi = surface_fun(this, y)
           
            psi = (y-this.c)'*this.inv_Sigma*(y-this.c) - 1;
            
        end
        
        function psi_grad = surface_fun_grad(this, y)
            
            psi_grad = 2*this.inv_Sigma*(y-this.c);
            
        end
        
        function calcSigma(this, angle, lambda_x, lambda_y)

            theta = angle * pi / 180;
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            Lambda_2 = diag([lambda_x, lambda_y]);
            this.Sigma = R * Lambda_2.^2 * R';
            this.inv_Sigma = inv(this.Sigma);

        end
        
        function points = generatePoints(this, c, Sigma, n_points)

            if (nargin < 4), n_points=200; end
            theta = linspace(0, 2*pi, n_points);
            L = chol(Sigma, 'lower'); % sqrtm(Sigma);
            points = c + L*[cos(theta); sin(theta)];

        end


    end
    
    
    
    properties (Access = protected)
       
        c % center
        Sigma % covariance matrix
        inv_Sigma
        
    end
    
end