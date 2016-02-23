function matparser(obj, varargin)
% assigns new values to existing matrices
%
% MATPARSER(obj, 'option1', value1 ...)
%
% The function modifies the sw.matrix.mat matrix, it assigns new values
% from a given parmeter vector.
%
% Input:
%
% obj           Input structure, spinw class object.
%
% Options:
%
% param         Input vector P with nPar elements that contains the
%               new values to be assignd to elements of sw.matrix.mat
%               matrix.
% mat           Identifies which matrices to be changed according to their
%               label or index. To select matrices with given labels use a
%               cell of strings with nPar elements, for example
%               M = {'J1','J2'}. This will change the diagonal elements of
%               matrices J1 and J2 to a given value that is provided in the
%               'param' option. Alternatively the index of the matrices can
%               be given in a vector, such as [1 2] (index runs according
%               to the order of the previous creation of the matrices using
%               sw.addmatrix).
%               To assign parameter value only to a selected element of a
%               3x3 matrix, a bracket notation can be used in any string,
%               such as 'D(3,3)', in this case only the (3,3) element of
%               the 3x3 matrix of 'D' will be modified, the other elements
%               will be unchanged. To modify multiple elements of a matrix
%               at once, use the option 'selector'.
% selector      Matrix with dimensions of [3 3 nPar]. Each S(:,:,ii)
%               submatrix can contain +/-1 and 0. Where S(:,:,ii) contains
%               ones, the corresponding matrix elements of
%               sw.matrix.mat(:,:,M(ii)) will be changed to the value
%               P(ii)*S(:,:,ii) where P(ii) is the corresponding parameter
%               value. For example do assign a Dzyaloshinskii-Moriya vector
%               to the 'DM' matrix, the following input would be
%               sufficient:
%               P = [0.2 0.35 3.14]
%               M = {'DM' 'DM' 'DM'}
%               S = cat(3,[0 0 0;0 0 1;0 -1 0],[0 0 -1;0 0 0;1 0 0],[0 1 0;-1 0 0;0 0 0])
%               sw.matparser('param',P,'mat',M,'selector',S)
% init          Initialize the matrices of sw.matrix.mat with zeros for all
%               selected labels before assigning paramter values. Default
%               is false.
%
% Output:
%
% The sw object will contain the modified matrix.mat field.
%
% See also SPINW, SPINW.HORACE, SPINW.ADDMATRIX.
%

inpForm.fname  = {'param' 'mat' 'selector' 'init'};
inpForm.defval = {[]      []     []         false};
inpForm.size   = {[1 -1]  [1 -1] []         [1 1]};
inpForm.soft   = {false   true   true       false};

param = sw_readparam(inpForm, varargin{:});

% number of parameters
nPar = numel(param.param);
% input parameters
P    = param.param;

if isempty(param.mat)
    if size(obj.matrix.mat,3)<nPar
        error('sw:matparser:WrongInput','Too many input parameters!');
    end
    M = 1:nPar;
else
    M = param.mat;
end

% check original matrix labels
br1Idx = strfind(obj.matrix.label,'(');
if any(cellfun(@(C)~isempty(C),br1Idx))
    error('sw:matparser:WrongLabel',['The sw object contains matrix labels'...
        ' with () symbols, can lead to ambiguous assignment of parameters!'...
        ' Please change the matrix labels of your sw object!'])
end

% default selector
S0 = zeros(3,3,nPar);

% convert cell of strings into matrix indices
if iscell(M)
    
    
    Mc = M;
    M = zeros(1,nPar);
    for ii = 1:nPar
        % remove text after brackets
        br1Idx = strfind(Mc{ii},'(');
        if ~isempty(br1Idx)
            Mc0 = Mc{ii}(1:(br1Idx(1)-1));
            % add index in default selector
            br2Idx = strfind(Mc{ii},')');
            idx = sscanf(Mc{ii}((br1Idx(1)+1):(br2Idx(1)-1)),'%d,%d');
            if numel(idx) ~= 2
                error('sw:matparser:WrongInput',['Wrong mat parameter, '...
                    'two integers are expected in the bracket!'])
            end
            S0(idx(1),idx(2),ii) = 1;
        else
            Mc0 = Mc{ii};
            S0(:,:,ii) = eye(3);
        end
        matIdx = strcmp(obj.matrix.label,Mc0);
        if sum(matIdx) == 0
            error('sw:matparser:WrongInput','The given matrix label does not exist!')
        elseif sum(matIdx) > 1
            error('sw:matparser:WrongInput','The given matrix label exist multiple times!')
        else
            M(ii) = find(matIdx);
        end
    end
else
    S0 = repmat(eye(3),[1 1 nPar]);
end

% auto create the selector
if isempty(param.selector)
    S = S0;
else
    S = param.selector;
end

if param.init
    obj.matrix.mat(:,:,unique(M)) = 0;
end
% assign matrices
for ii = 1:nPar
    obj.matrix.mat(:,:,M(ii)) = obj.matrix.mat(:,:,M(ii)).*sign((1-abs(S(:,:,ii)))) + S(:,:,ii).*P(ii);
end

end