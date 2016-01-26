function varargout = disp(obj, fid)
% prints the spinw data structure in readable format onto the Command Window
%
% {swDescr} = DISPLAY(obj, {fid})
%
% Input:
%
% obj       spinw class object.
% fid       File identifier (created by fopen function) where the output
%           will be written. Default is stored in obj, the output will be
%           written onto the Command Line. Optional.
%
% Output:
%
% swdescr   If output variable is given, the description of the obj object
%           will be output into the swdescr variable, instead of being
%           written onto the Command Window/file. Optional.
%
% Example:
%
% crystal = spinw;
% swFields = display(crystal);
%
% See also SPINW.
%

Datastruct = datastruct;

choiceStr = {'off' 'on'};
symbStr = choiceStr{obj.symbolic+1};
symmStr = choiceStr{obj.symmetry+1};


if nargout == 1
    swDescr = sprintf('spinw object (symbolic: %s, symmetry: %s, textoutput: %s)\n',symbStr,symmStr,fopen(obj.fileid));
    
    
    for ii = 1:length(Datastruct.mainfield)
        swDescr = [swDescr sprintf('%s\n', Datastruct.mainfield{ii})]; %#ok<*AGROW>
        index = 1;
        plotSize = true;
        while index <= size(Datastruct.subfield,2) && ~isempty(Datastruct.subfield{ii,index})
            swDescr = [swDescr sprintf('    %10s: ',Datastruct.subfield{ii,index})];
            sizefield = Datastruct.sizefield{ii,index};
            typefield = Datastruct.typefield{ii,index};
            
            swDescr = [swDescr sprintf('[')];
            
            for jj = 1:length(sizefield)-1
                objSize = [size(obj.(Datastruct.mainfield{ii}).(Datastruct.subfield{ii,index})) 1 1 1];
                if ischar(sizefield{jj})
                    if plotSize
                        swDescr = [swDescr sprintf('%sx', sizefield{jj})];
                        sF = sizefield{jj};
                        sS = objSize(jj);
                    else
                        swDescr = [swDescr sprintf('%sx', sizefield{jj})];
                    end
                    
                else
                    swDescr = [swDescr sprintf('%dx', sizefield{jj})];
                end
            end
            
            if ischar(sizefield{jj+1})
                if plotSize
                    swDescr = [swDescr sprintf('%s', sizefield{jj+1})];
                    sF = sizefield{jj+1};
                    sS = objSize(jj+1);
                    
                else
                    swDescr = [swDescr sprintf('%s', sizefield{jj+1})];
                end
                
            else
                swDescr = [swDescr sprintf('%d', sizefield{jj+1})];
            end
            swDescr = [swDescr sprintf(' %s]',typefield)];
            if plotSize && exist('sF','var')
                swDescr = [swDescr sprintf('  %s=%d',sF,sS)];
                clear sF sS
                plotSize = false;
            end
            swDescr = [swDescr sprintf('\n')];
            index = index + 1;
        end
    end
    
    varargout{1} = swDescr;
else
    if nargin < 2
        fid = obj.fid;
    end
    
    chem = obj.formula;
    chem.chemform(chem.chemform == ' ') = [];
    chem.chemform(chem.chemform == '_') = [];
    abc  = obj.abc;
    aa   = symbol('angs');
    [~,~,symStr]=sw_gensym(obj.lattice.sym);
    
    swDescr = ['\b<strong>SpinW</strong> object, <a href="matlab:doc @spinw">spinw</a> class:\n'...
        sprintf('<strong>Chemical formula</strong>: %s\n',chem.chemform) ...
        sprintf('<strong>Space group</strong>:      %s\n',symStr)...
        sprintf(['<strong>Lattice</strong>:\n  a=%7.4f ' aa ', b=%7.4f ' aa ', c=%7.4f ' aa '\n'],abc(1),abc(2),abc(3))...
        sprintf(['  ' char(945) '=%6.2f' char(176) ',   ' char(946) '=%6.2f' char(176) ',   ' char(947) '=%6.2f' char(176) '\n'],abc(4),abc(5),abc(6))...
        sprintf('<strong>Magnetic atoms in the cell</strong>: %d\n',numel(obj.matom.S))...
        sprintf('<strong>Mode</strong>:\n  symbolic: %s, symmetry: %s, textoutput: %s\n',symbStr,symmStr,fopen(obj.fileid))...
        ];
    
    % print the text
    fprintf0(fid,swDescr);
end

end