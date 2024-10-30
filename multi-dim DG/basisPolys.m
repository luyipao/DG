classdef basisPolys
    properties
        X (:,1) double
        Xc (:,1) double
        h
        coeffs (:,:) double
        degree (1,1)
        basisFuncs
        basisInterval
		basisBoundaryValues
    end
    methods
        function obj = basisPolys(mesh, coeffs, degree, basisFuncs, interval)
            % intput:
            %   mesh: 1-D mesh
            %   coeffs: ith colum represent the coeffs of basisFuncs;
            %   degree: max degree of basisFuncs;
            %   basisFuncs: 1-D cell, basis is in the standard
            %   interval;baisFuncs{i}(x) should cal correctly, where x is
            %   vector;
            %   interval: basisFuncs based on this interval: you should cal
            %   based on mesh.;
            % output:
            %   obj
            if size(coeffs) ~= [degree+1,length(mesh)-1]
                error('input error: size ');
            else
                obj.X = mesh;
                obj.Xc = 0.5*(mesh(1:end-1) + mesh(2:end));
                obj.degree = degree;
                obj.basisFuncs = basisFuncs;
                obj.basisInterval = interval;
				obj.basisBoundaryValues = zeros(obj.degree+1, 2);

                assert(isnumeric(obj.degree) && obj.degree >= 0, 'obj.degree 应为非负整数');
                assert(iscell(obj.basisFuncs), 'obj.basisFuncs 应该是一个 cell 数组');
                for idx = 1:numel(obj.basisFuncs)
                    assert(isa(obj.basisFuncs{idx}, 'function_handle'), 'obj.basisFuncs 中的每个元素应为函数句柄');
                end
                assert(isnumeric(obj.basisInterval) && numel(obj.basisInterval) == 2, 'obj.basisInterval 应为包含两个数值元素的数组');
				
                
                for j = 1:obj.degree+1
                    obj.basisBoundaryValues(j,:) = obj.basisFuncs{j}(obj.basisInterval);
                end

                obj.h = abs(mesh(2) - mesh(1));
                obj.coeffs = coeffs;
            end
        end
        function y = getNodeValues(obj)
            % output: y(2, length(mesh) - 1) or y(2, cellNum)
            %   y(1,i) : the value of f, on the left endpoint of I_i;
            %   y(2,i): the value of f, one the right endpoint of I_i;
			A = repmat(obj.coeffs,2,1).*obj.basisBoundaryValues(:);
			y(1,:) = sum(A(1:obj.degree+1,:),1);
            y(2,:) = sum(A(obj.degree+2:end,:),1);
			
		end
        function y = solve(obj, x)
            % input: x can be 2-D matrix at most; output: y = f(x);
            [m,n] = size(x);
            x = x(:);
			% scale
			[x, index] = transform(obj, x);
			
            y = obj.coeffs(:,index) .* cell2mat(cellfun(@(f) f(x), obj.basisFuncs, 'UniformOutput', false))';
            
            y = reshape(sum(y,1), m, n);
        end

        function [x, index] = transform(obj, x)
        % transform function: Projects x onto the equivalent points on the support interval of the basis functions
        % Input: x is the vector of points to be transformed
        % Output: x is the vector of transformed points
            index = discretize(x,obj.X);
            x = (x-obj.Xc(index))/obj.h;
        end
    end
    
end