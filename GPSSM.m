classdef GPSSM
% class for concatenating 2 GP predictions
% Last edited: Armin Lederer, 04/2020
    properties
        gp1;
        gp2;
    end
    methods
        function obj = GPSSM(gp1,gp2)
        % initialize object
        % In:
        %   gp1        obj     first GP moddel
        %   gp2        obj     second GP model
        % Out:
        %   obj        obj     two GPs model
            obj.gp1=gp1;
            obj.gp2=gp2;
        end
        function [mu,sig]=predict(obj,x)
        % predict GP models
        % In:
        %   obj        1  x 1  class object
        %   x          D  x 1  prediction input
        % Out:
        %   mu         2  x 1  posterior mean of both GPs
        %   sig        2  x 1  posterior variance of both GPs
            mu=zeros(2,1);
            sig=zeros(2,1);
            [mu(1),sig(1)]=obj.gp1.predict(x);
            [mu(2),sig(2)]=obj.gp2.predict(x);
        end
    end

end