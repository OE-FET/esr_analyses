classdef nelder_mead_fitobject < handle

    properties
        fitfunc
        coef0
        coef
        sse
        J
        xData
        yData
        yFit
    end

    methods
        function obj = nelder_mead_fitobject(fitfunc, xData, yData, coef0, coef, sse)
            
            obj.fitfunc = fitfunc;
            obj.xData = xData;
            obj.yData = yData;
            obj.coef0 = coef0;
            obj.coef = coef;
            obj.sse = sse;
            obj.yFit = eval_at(obj, obj.xData);
        
        end
        
        function J = jacobian(obj, accur)
            
            import esr_analyses.*
            import esr_analyses.utils.*
            
            if isempty(obj.J)

                if nargin < 2; accur = 'accurate'; end

                % get jacobian matrix with respect to fit parameters
                func = @(coef) obj.fitfunc(coef, obj.xData);
                if strcmp(accur, 'quick')
                    % use lsqnonlin with single iteration to get J (quick)
                    [~,~,~,~,~,~,obj.J] = lsqnonlin(func, obj.coef,[],[], optimset('Display', 'off'));
                elseif strcmp(accur, 'accurate')
                    % use jacobianest from spinach toolbox (accurate but slow)
                    obj.J = jacobianest(func, obj.coef);
                end
                J = obj.J;
            end
        end
        
        function ci = confint(obj, varargin)
            % CONFINT of fit paramters
            %
            %   Returns the 95% confidence intervals for fit coefficients.
            %
            
            import esr_analyses.utils.*
            
            jacobian(obj);

            % calculate residuals
            resid = obj.yData - obj.yFit;
            
            alpha = get_kwarg(varargin, 'alpha', 0.05);
            ci = nlparci(obj.coef, resid, 'jacobian', obj.J, 'alpha', alpha);

        end
        
        function se = standarderror(obj, varargin)
            % CONFINT of fit paramters
            %
            %   Returns the 95% confidence intervals for fit coefficients.
            %

            ci = confint(obj);
            se = diff(ci,1,2)/2;

        end
        
        function h = plot(obj)

            import esr_analyses.*
            import esr_analyses.utils.*

            if iscell(obj.xData)
                if ~length(obj.xData) == 2
                    error('Can only plot 2D or 3D data.')
                end

                x1       = obj.xData{1};
                x2       = obj.xData{2};

                xx1 = x1(1,:);
                xx2 = x2(:,1);
                % get best-fit curve (in a higher resolution version)
                x_interp        = linspace(min(xx1), max(xx1), 2^10);
                x2_interp        = linspace(min(xx2), max(xx2), 2^10);
                [X1Plot, X2Plot] = meshgrid(x_interp, x2_interp);
                yFitInterpol     = eval_at(obj, {x1, x2});
                yFitMesh         = eval_at(obj, {X1Plot, X2Plot});

                % plot best-fit curve and data
                X1t = transpose(x1); X2t = transpose(x2);

                figure('Name', '3D fit');
                hold on;
                h1 = scatter3(X1t(:), X2t(:), obj.yData(:), '.k');
                h2 = surf(X1Plot, X2Plot, yFitMesh', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                legend('Data', 'Fit', 'Location', 'northeast')
                set(gca, 'XLimSpec', 'Tight'); grid on;
                view(3);

                figure('Name', '3D fit, x-section');
                hold on;
                [h3, yoffsets] = stackplot(xx1, obj.yData, 'style', 'b');
                h4 = stackplot(xx1, yFitInterpol, 'yoffsets', yoffsets, 'style', 'r');
                legend([h3(1), h4(1)], {'Data', 'Fit'});
                if nargout > 0
                    h = {h1 h2 h3 h4};
                end
            else
                x_interp = transpose(linspace(min(obj.xData), max(obj.xData), 2^10));
                yFit_interp = eval_at(obj, x_interp);

                figure('Name', 'Least-squares fit');
                hold on;
                h1 = plot(obj.xData, obj.yData);
                h2 = plot(x_interp, yFit_interp);
                legend([h1, h2], {'Data', 'Fit'});
                set(gca, 'XLimSpec', 'Tight'); grid on;
                if nargout > 0
                    h = {h1 h2};
                end
            end

        end
        
        function y = eval_at(obj, x)
            y = obj.fitfunc(obj.coef, x);
        end
    end
end
