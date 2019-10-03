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
        function se = standarderror(obj, accur)
            % STANDARDERROR of of fit paramters
            %
            % We calculate the standard-error (SE) of fit parameters as
            %
            %   SE = sqrt( sigma_y * inv(H) )
            %
            % where 'sqrt(sigma_y)' is the standard deviation of residuals
            % and 'H' is the Hessian. 'sigma_y' can be calculated as
            %
            %   sigma_y = RSS / dof
            %
            % where 'RSS' is the resiual sum-of-squares and 'dof' is the
            % number of degrees-of-freedom in the fitting problem. The
            % Hessian can be estimated from the Jacobian 'J' of the fit
            % function with respect to the fitting parameters as H ~ J'*J.

            import esr_analyses.*
            import esr_analyses.utils.*

            if nargin < 2; accur = 'accurate'; end

            dof     = numel(obj.dependent_fitdata) - numel(obj.coef0);
            sigma_y = obj.sse/dof;

            % get jacobian matrix with respect to fit parameters
            func = @(coef) obj.fitfunc(coef, obj.independent_fitdata);
            if strcmp(accur, 'quick')
                % use lsqnonlin with single iteration to get J (quick)
                [~,~,~,~,~,~,J] = lsqnonlin(func, obj.coef,[],[], optimset('Display', 'off'));
            elseif strcmp(accur, 'accurate')
                % use jacobianest from spinach toolbox (accurate but slow)
                J = jacobianest(func, obj.coef);
            end
            % decomposition J = Q*R with upper triangular matrix R and unitary matrix Q
            [~, R] = qr(J, 0);
            % diagnonal of covariance matrix Sigma = sigma_y*inv(J'*J)
            diag_sigma = sigma_y * sum(inv(R).^2, 2);
            % diag_sigma = sigma_y*inv(J'*J);
            % parameter standrad errors
            se = sqrt(diag_sigma)';
            se = full(se); % convert sparse to full matrix
        end
        function h = plot(obj)

            import esr_analyses.*
            import esr_analyses.utils.*

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
                zFit            = obj.fitfunc(obj.coef, {X, Y});
                zFitMesh        = obj.fitfunc(obj.coef, {XPlot, YPlot});

                % plot best-fit curve and data
                Xt = X'; Yt = Y';

                figure('Name', '3D fit');
                hold on;
                h1 = scatter3(Xt(:), Yt(:), zData(:), '.k');
                h2 = surf(XPlot, YPlot, zFitMesh', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                legend('Data', 'Fit', 'Location', 'northeast')
                axis tight; grid on;
                view(3);

                figure('Name', '3D fit, x-section');
                hold on;
                [h3, yoffsets] = stack_plot(x, zData, 'style', '.k');
                h4 = stack_plot(x, zFit, 'yoffsets', yoffsets, 'style', '-r');
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
