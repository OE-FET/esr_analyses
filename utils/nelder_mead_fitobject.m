classdef nelder_mead_fitobject
    
    properties
        fitfunc
        coef0
        coef
        sse
        independent_fitdata
        dependent_fitdata
    end
    
    methods
        function ci = confint(obj, accur)
            if nargin < 2; accur = 'accurate'; end
            % degrees of freedom in fitting problem
            dof     = numel(obj.dependent_fitdata) - numel(obj.coef0);
            % standard deviation of residuals
            sdr     = sqrt(obj.sse/dof);
            % jacobian matrix
            rff     = @(coef) obj.fitfunc(coef, obj.independent_fitdata);
            if strcmp(accur, 'quick')
                % use lsqnonlin with single iteration, much quicker
                [~,~,~,~,~,~,J] = lsqnonlin(rff, obj.coef,[],[], ...
                    optimset('MaxFunEvals', 0, 'Display', 'off'));
            elseif strcmp(accur, 'accurate')
                % use jacobianest from spinach toolbar, more accurate
                J = jacobianest(rff, obj.coef);
            end
            % decomposition J = Q*R with upper triangular matrix R and unitary matrix Q
            [~, R] = qr(J, 0);
            % diagnonal of covariance matrix Sigma = sdr^2*inv(J'*J)
            diag_sigma = sdr^2*sum(inv(R).^2, 2);
            % parameter standrad errors
            ci      = sqrt(diag_sigma)';
        end
        function h = plot(obj)
            
            if iscell(obj.independent_fitdata)
                if ~length(obj.independent_fitdata) == 2
                    error('Can only plot 2D or 3D data.')
                end
                
                X       = obj.independent_fitdata{1};
                Y       = obj.independent_fitdata{2};
                zData   = obj.dependent_fitdata;
                
                x = X(1,:);
                y = Y(:,1);
                % get best-fit curve (in a higher resolution version)
                x_interp        = linspace(min(x), max(x), 2^10);
                y_interp        = linspace(min(y), max(y), 2^10);
                [XPlot, YPlot]  = meshgrid(x_interp, y_interp);
                zFit            = obj.fitfunc(obj.coef, {XPlot, YPlot});
                
                % plot best-fit curve and data
                Xt = X'; Yt = Y';

                figure('Name', '3D fit');
                hold on;
                h1 = scatter3(Xt(:), Yt(:), zData(:), '.k');
                h2 = surf(XPlot, YPlot, zFit', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                legend('Data', 'Fit', 'Location', 'northeast')
                axis tight; grid on;
                view(3);

                figure('Name', '3D fit, x-section');
                offset = max(max(zData))*0.8;
                hold on;
                h3 = stackplot(x, zData, 'yoffset', offset, 'style', '.k');
                h4 = stackplot(x, obj.fitfunc(obj.coef, {X, Y}), ...
                    'yoffset', offset, 'style', '-r');
                legend([h3(1), h4(1)], {'Data', 'Fit'});
                if nargout > 0
                    h = {h1 h2 h3 h4};
                end
            else
                x = obj.independent_fitdata;
                zData = obj.dependent_fitdata;
                
                x_interp = linspace(min(x), max(x), 2^10)';
                zFit = obj.fitfunc(obj.coef, x_interp);
                
                figure('Name', 'Least-squares fit');
                hold on;
                h1 = plot(x, zData, '.k');
                h2 = plot(x_interp, zFit, '-r');
                legend([h1, h2], {'Data', 'Fit'});
                axis tight; grid on;
                if nargout > 0
                    h = {h1 h2};
                end
            end
                
        end
    end
end
