classdef GammaCounts
    % activity Object to calculate an instantaneous activity from measurement
    %   Detailed explanation goes here
    
    
    
    properties
        %
        %         Nominal values
        
        %         Gamma-ray energy, in keV
        EGamma
        Lambda
        %         Decay half-life, in seconds
        HalfLife
        %         Number of observed counts, from peak fitting
        NumberOfCounts
        Efficiency
        IGamma
        %         Counting live time, in seconds
        LiveTime
        Mu
        %         Areal density, in mg/cm^2
        RhoDr
        EfficiencyCovarianceData
        %         Activity, in Bq
        Activity
        Mass = NaN;
        CountStartTime = 0;
        EoBTime
        PCov
        POpt
        FileName
        %
        %         Boolean for returning either activity or number of counts
        %
        IsCumulative
        %
        %         Uncertainties in nominal values
        %
        UncertaintyHalfLife
        UncertaintyLambda
        UncertaintyNumberOfCounts
        UncertaintyEfficiency
        UncertaintyIGamma
        UncertaintyLiveTime
        UncertaintyMu
        UncertaintyRhoDr
        %         UncertaintyActivity;
    end
    
    properties (Dependent)
        UncertaintyNumberofDecays
        NumberofDecays
        TimeSinceEoB
    end
    
    properties (SetAccess = private)
        UncertaintyActivity = 5;
    end
    
    methods
        function obj = GammaCounts( E_gamma, t_half, number_of_counts, unc_number_of_counts, I_gamma, unc_I_gamma, live_time, mu, rhodr, file_name, mass, count_start_time, eob_time, unc_rhodr, covariance_data)
            %Construct an instance of this class
            if nargin == 2
                % Make 2D oject array
                obj(E_gamma, t_half) = obj;
            elseif nargin > 2
                % Construct objects
                obj.FileName = file_name;
                obj.EGamma = E_gamma;
                obj.HalfLife = t_half;
                obj.Lambda = (log(2)) ./  obj.HalfLife;
                obj.NumberOfCounts = number_of_counts;
                %                 obj.Efficiency = efficiency;
                obj.IGamma = I_gamma.*.01;
                obj.LiveTime = live_time;
                obj.Mu = mu;
                obj.RhoDr = rhodr;
                obj.Mass = mass;
                obj.CountStartTime = count_start_time;
                obj.EoBTime = eob_time;
                obj.EfficiencyCovarianceData = covariance_data;
                obj.UncertaintyNumberOfCounts = unc_number_of_counts .* obj.NumberOfCounts ./ 100;
                obj.UncertaintyIGamma = unc_I_gamma.*.01;
                obj.UncertaintyRhoDr = unc_rhodr;
                
                % Take 5% uncertainty since none listed in XCOM
                obj.UncertaintyMu = 0.05 .*mu;
                
                % Assume very tiny (0.1%) uncertainty in lifetime
                obj.UncertaintyHalfLife = obj.HalfLife .* 0.001;
                obj.UncertaintyLambda = obj.Lambda .* (obj.UncertaintyHalfLife ./ obj.HalfLife);
                
                % Assume 2s uncertainty in live time
                obj.UncertaintyLiveTime = 2;
                
                % Calculate efficiency
                [eff, unc_eff, pcov, popt] = efficiency_calibration(obj.EGamma, obj.EfficiencyCovarianceData);
                obj.Efficiency = eff;
                obj.UncertaintyEfficiency = unc_eff;
                obj.PCov = pcov;
                obj.POpt = popt;
                
                % Calculate activity automatically
                %                 fprintf("Activity:",obj.calcActivity)
                obj.Activity = obj.calcActivity;
                
                %                 obj.UncertaintyActivity = 5;
                
            else
                fprintf('Warning: object not initialized!\n\n')
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        function obj = set.Activity(obj,val)
            % Check if Activity has not been set
            if isempty(obj.Activity)
                obj.Activity = obj.calcActivity;
                %                 obj.updateUncertaintyActivity;
            else
                obj.Activity = val;
                %                 obj.updateUncertaintyActivity;
                %                 obj.calcUncertaintyActivity
            end
            %             updateUncertaintyActivity(obj);
        end
        %         function obj = calcUncertaintyActivity(obj)
        %             obj.UncertaintyActivity = 2;
        % %             fprintf('test')
        %         end
        
        function value = calcActivity(obj,varargin)
            %
            %
            %   A(t)  =                         lambda * number_of_counts
            %            ----------------------------------------------------------------------
            %            efficiency * I_gamma * (1-exp(-lambda*live_time)) * (exp(-mu*rhodr/2))
            %
            %
            % my_activity = (lambda .* number_of_counts) ./ (efficiency .* I_gamma .* (1-exp(-lambda.*live_time)) .* (exp(-mu .* rhodr .* 0.5e-3)));
            if nargin > 1
                %obj.calcActivity(lambda,numCounts,eff,i_gamma,lt,mu,rhodr)
                value = (varargin{1} .* varargin{2}) ./ (varargin{3} .* varargin{4} .* (1-exp(-varargin{1}.*varargin{5})) .* (exp(-varargin{6} .* varargin{7} .* 0.5e-3)));
            else
                value = (obj.Lambda .* obj.NumberOfCounts) ./ (obj.Efficiency .* obj.IGamma .* (1-exp(-obj.Lambda.*obj.LiveTime)) .* (exp(-obj.Mu .* obj.RhoDr .* 0.5e-3)));
                %                 obj.Activity = value;
            end
            
            %             obj.Activity = value;
            
        end
        function value = get.UncertaintyActivity(obj)
            % Numerical derivative approximation for uncorrelated
            % uncertainty propagation in instantaneous activity
            
            %             obj.Activity = obj.calcActivity;
            
            % Set machine tolerance for numerical derivatives
            % This is a relative scaling parameter!
            epsilon = 1E-2;
            
            % Get nominal values of parameters for activity
            lambda = obj.Lambda;
            numCounts = obj.NumberOfCounts;
            eff = obj.Efficiency;
            i_gamma = obj.IGamma;
            lt = obj.LiveTime;
            mu = obj.Mu;
            rhodr = obj.RhoDr;
            
            delta_x = [lambda,numCounts,eff,i_gamma,lt,mu,rhodr].*epsilon;
            
            lambda = ones(2,7).*lambda;
            numCounts = ones(2,7).*numCounts;
            eff = ones(2,7).*eff;
            i_gamma = ones(2,7).*i_gamma;
            lt = ones(2,7).*lt;
            mu = ones(2,7).*mu;
            rhodr = ones(2,7).*rhodr;
            
            
            % Perturb parameters
            lambda(:,1) = [obj.Lambda + delta_x(1)./2 ;obj.Lambda - delta_x(1)./2];
            numCounts(:,2) = [obj.NumberOfCounts + delta_x(2)./2 ;obj.NumberOfCounts - delta_x(2)./2];
            eff(:,3) = [obj.Efficiency + delta_x(3)./2 ;obj.Efficiency - delta_x(3)./2];
            i_gamma(:,4) = [obj.IGamma + delta_x(4)./2 ;obj.IGamma - delta_x(4)./2];
            lt(:,5) = [obj.LiveTime + delta_x(5)./2 ;obj.LiveTime - delta_x(5)./2];
            mu(:,6) = [obj.Mu + delta_x(6)./2 ;obj.Mu - delta_x(6)./2];
            rhodr(:,7) = [obj.RhoDr + delta_x(7)./2 ;obj.RhoDr - delta_x(7)./2];
            
            
            perturbed_activities = obj.calcActivity(lambda,numCounts,eff,i_gamma,lt,mu,rhodr);
            partial_derivaties = (perturbed_activities(1,:) - perturbed_activities(2,:)) ./delta_x;
            %             (lambda .* numCounts) ./ (eff .* i_gamma .* (1-exp(-lambda.*lt)) .* (exp(-mu .* rhodr .* 0.5e-3)))
            
            value = sqrt(sum(partial_derivaties.^2  .*  [obj.UncertaintyLambda,obj.UncertaintyNumberOfCounts,obj.UncertaintyEfficiency,obj.UncertaintyIGamma,obj.UncertaintyLiveTime,obj.UncertaintyMu,obj.UncertaintyRhoDr].^2  ));
            
            
            %             value = obj.LiveTime + obj.IGamma;
        end
        function value = get.NumberofDecays(obj)
            value = obj.NumberOfCounts / (obj.IGamma * obj.Efficiency);
        end
        function value = get.UncertaintyNumberofDecays(obj)
            %             NumberofDecays = obj.NumberOfCounts / (obj.IGamma * obj.Efficiency);
            
            value = sqrt((obj.UncertaintyNumberOfCounts ./ (obj.IGamma .* obj.Efficiency)).^2 ...
                + (obj.NumberOfCounts ./ ((obj.IGamma).^2 .* obj.Efficiency)).^2 .* (obj.UncertaintyIGamma).^2 ...
                + (obj.NumberOfCounts ./ ((obj.Efficiency).^2 .* obj.IGamma)).^2 .* (obj.UncertaintyEfficiency).^2);
        end
        function value = get.TimeSinceEoB(obj)
          value = (obj.CountStartTime - juliandate(datetime(obj.EoBTime)) ) .* 24;
        end
        %     end
        %
        %     methods(Access = private)
        %         function   updateUncertaintyActivity(obj)
        %             obj.UncertaintyActivity = 2;
        % %             fprintf('test')
        %         end
    end
end

