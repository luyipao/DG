
function outputDoping = initialFunction(inputPosition)
%% INITIALFUNCTION 计算400nm沟道器件的维度掺杂函数N_D
% 这是 $\Phi$ 的初始值，即：
% $\Phi(0,x,w,\mu) = s(w) N_D(x) e^{-w} \mathcal{M}$
% 其中 $\mathcal{M}$ 被选择使得密度的初始值等于掺杂 $N_D(x)$ 
%
% 输入:
%   inputPosition - 可以是2D矩阵， $inputPosition \in [0, 1]$
%
% 输出: 
%   outputDoping - 是inputPosition的初始函数结果，即 $N_D(x)$
% 
%  示例
% inputPosition = linspace(0, 1, 1000);
% outputDoping = initialFunction(inputPosition);
% hold on
% plot(inputPosition, outputDoping);
% legend("initial function");
% hold off

% 定义常量
transitionWidth = 0.03;
highDopingLevel = 5e17;  % 5 * 10^17 cm^-3
lowDopingLevel = 1e15;   % 1 * 10^15 cm^-3
dopingDifference = highDopingLevel - lowDopingLevel;

% 预分配输出数组
outputDoping = zeros(size(inputPosition));

% 使用逻辑索引定义区域
lowDopingRegion = inputPosition <= 0.2 - transitionWidth | inputPosition >= 0.4 + transitionWidth;
midDopingRegion = inputPosition >= 0.2 + transitionWidth & inputPosition <= 0.4 - transitionWidth;
leftTransitionRegion = abs(inputPosition - 0.2) < transitionWidth;
rightTransitionRegion = abs(inputPosition - 0.4) < transitionWidth;

% 应用条件
outputDoping(lowDopingRegion) = highDopingLevel;
outputDoping(midDopingRegion) = lowDopingLevel;

% 处理过渡区域
leftTransitionFactor = (inputPosition(leftTransitionRegion) - 0.2 + transitionWidth) / (2 * transitionWidth);
rightTransitionFactor = (-inputPosition(rightTransitionRegion) + 0.4 + transitionWidth) / (2 * transitionWidth);

outputDoping(leftTransitionRegion) = dopingDifference * (1 - leftTransitionFactor.^3).^3 + lowDopingLevel;
outputDoping(rightTransitionRegion) = dopingDifference * (1 - rightTransitionFactor.^3).^3 + lowDopingLevel;

end
