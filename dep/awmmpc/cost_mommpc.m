function F = cost_mommpc(Yi,T,U,models,mpc)
%COST_MOMMPC_WLSE Weighted least-squares objective for multi-model control.
%   The model and target trajectories are normalized to ensure they are
%   similar in terms of gain and offset. The objective is defined as the
%   sum of squared differences between trajectories and targets and the sum
%   of squared input quantities.
%   
%   INPUTS:
%   Yi: [cell] contains trajectories for the models, has the form:
%       { Y_model1(outputs,t) , ... , Y_modelN(outputs,t) }.
%   T: [vector] time points cooresponding to the trajectories.
%   U: [vector] control inputs.
%   models: [structure] contains model information (e.g. model identifier,
%       initial conditions, control/output variable indices, etc.).
%   mpc: [structure] contains algorithm properties (e.g. prediction
%       horizon, duration, switching times, etc.).
%   
%   OUTPUTS:
%   F: [vector] objective values.
%   
%   Written by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 5/17/2012


% Extract working variables
Nm = length(models);        % Number of prediction models
s = mpc.alg.s;              % Profile of target outputs
Q = mpc.alg.Q;              % Weighting for model deviation
R = mpc.alg.R;              % Weighting for input energy
convertToLog10 = mpc.cost.convertToLog10;% Log10 conversion
Ny = arrayfun(@(x)size(x.C,1),models);% Number of outputs per model

% Extract trajectories and scale for qualitative comparison
Y = []; S = []; s = feval(s,T(:)'); U = U(:); f = zeros(Nm,1);
for j = 1:Nm                % FOR each prediction model
    % Select trajectories and generate reference profiles
    yi = Yi{j}; si = s(1:Ny(j),:);
    % Offset and gain matching for model/target outputs
    A = models(j).A; B = models(j).B; one = ones(1,size(yi,2));
    Yn = (yi - B(:,1)*one)./(A(:,1)*one); Yn = Yn(:); Y = cat(1,Y,Yn);
    Sn = (si - B(:,2)*one)./(A(:,2)*one); Sn = Sn(:); S = cat(1,S,Sn);
    % Evaluate objective function
    f(j) = (Yn-Sn)'*Q*(Yn-Sn) + U'*R*U;
end                         % FOR each prediction model

% Evaluate objective and return outputs
F = f(:)';                  % Weighted objective values separated by model
if convertToLog10, F = log10(F + 1); end% Convert for smoothness

end