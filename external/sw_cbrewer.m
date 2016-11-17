function cMap = sw_cbrewer(cName,nCol)
% creates nice colormaps based on the original cbrewer() function
%
% cMap = SW_CBREWER(colorName, nColors)
%
% This function produces nice colormaps, for more information on the colors
% please visit http://colorbrewer2.org/
%
% Input:
%
% colorName Name of the colormap in a string, for the possible names, call
%           sw_cbrewer without input and the possible colormaps will be
%           shown.
% nColors   Number of colors in the colormap. The colors will be
%           interpolated from the few predefined base colors.
%
%  This function includes color specifications and designs developed by
%  Cynthia Brewer (http://colorbrewer.org/).
%

% Original cbrewer author: Charles Robert
% email: tannoudji@hotmail.com
% Date: 06.12.2011
%

% initialize the data with the given values from the original cBrewer
cDat.Accent   = [127 201 127;190 174 212;253 192 134;255 255 153;56 108 176;240 2 127;191 91 23;191 91 23];
cDat.Dark2    = [27 158 119;217 95 2;117 112 179;231 41 138;102 166 30;230 171 2;166 118 29;166 118 29];
cDat.Paired   = [166 206 227;31 120 180;178 223 138;51 160 44;251 154 153;227 26 28;253 191 111;255 127 0;202 178 214;106 61 154;255 255 153;255 255 153];
cDat.Pastel1  = [251 180 174;179 205 227;204 235 197;222 203 228;254 217 166;255 255 204;229 216 189;253 218 236;253 218 236];
cDat.Pastel2  = [179 226 205;253 205 172;203 213 232;244 202 228;230 245 201;255 242 174;241 226 204;241 226 204];
cDat.Set1     = [228 26 28;55 126 184;77 175 74;152 78 163;255 127 0;255 255 51;166 86 40;247 129 191;247 129 191];
cDat.Set2     = [102 194 165;252 141 98;141 160 203;231 138 195;166 216 84;255 217 47;229 196 148;229 196 148];
cDat.Set3     = [141 211 199;255 255 179;190 186 218;251 128 114;128 177 211;253 180 98;179 222 105;252 205 229;217 217 217;188 128 189;204 235 197;204 235 197];
cDat.Blues    = [247 251 255;222 235 247;198 219 239;158 202 225;107 174 214;66 146 198;33 113 181;8 81 156;8 81 156];
cDat.BuGn     = [247 252 253;229 245 249;204 236 230;153 216 201;102 194 164;65 174 118;35 139 69;0 109 44;0 109 44];
cDat.BuPu     = [247 252 253;224 236 244;191 211 230;158 188 218;140 150 198;140 107 177;136 65 157;129 15 124;129 15 124];
cDat.GnBu     = [247 252 240;224 243 219;204 235 197;168 221 181;123 204 196;78 179 211;43 140 190;8 104 172;8 104 172];
cDat.Greens   = [247 252 245;229 245 224;199 233 192;161 217 155;116 196 118;65 171 93;35 139 69;0 109 44;0 109 44];
cDat.Greys    = [255 255 255;240 240 240;217 217 217;189 189 189;150 150 150;115 115 115;82 82 82;37 37 37;37 37 37];
cDat.Oranges  = [255 245 235;254 230 206;253 208 162;253 174 107;253 141 60;241 105 19;217 72 1;166 54 3;166 54 3];
cDat.OrRd     = [255 247 236;254 232 200;253 212 158;253 187 132;252 141 89;239 101 72;215 48 31;179 0 0;179 0 0];
cDat.PuBu     = [255 247 251;236 231 242;208 209 230;166 189 219;116 169 207;54 144 192;5 112 176;4 90 141;4 90 141];
cDat.PuBuGn   = [255 247 251;236 226 240;208 209 230;166 189 219;103 169 207;54 144 192;2 129 138;1 108 89;1 108 89];
cDat.PuRd     = [247 244 249;231 225 239;212 185 218;201 148 199;223 101 176;231 41 138;206 18 86;152 0 67;152 0 67];
cDat.Purples  = [252 251 253;239 237 245;218 218 235;188 189 220;158 154 200;128 125 186;106 81 163;84 39 143;84 39 143];
cDat.RdPu     = [255 247 243;253 224 221;252 197 192;250 159 181;247 104 161;221 52 151;174 1 126;122 1 119;122 1 119];
cDat.Reds     = [255 245 240;254 224 210;252 187 161;252 146 114;251 106 74;239 59 44;203 24 29;165 15 21;165 15 21];
cDat.YlGn     = [255 255 229;247 252 185;217 240 163;173 221 142;120 198 121;65 171 93;35 132 67;0 104 55;0 104 55];
cDat.YlGnBu   = [255 255 217;237 248 177;199 233 180;127 205 187;65 182 196;29 145 192;34 94 168;37 52 148;37 52 148];
cDat.YlOrBr   = [255 255 229;255 247 188;254 227 145;254 196 79;254 153 41;236 112 20;204 76 2;153 52 4;153 52 4];
cDat.YlOrRd   = [255 255 204;255 237 160;254 217 118;254 178 76;253 141 60;252 78 42;227 26 28;189 0 38;189 0 38];
cDat.BrBG     = [84 48 5;140 81 10;191 129 45;223 194 125;246 232 195;245 245 245;199 234 229;128 205 193;53 151 143;1 102 94;1 102 94];
cDat.PiYG     = [142 1 82;197 27 125;222 119 174;241 182 218;253 224 239;247 247 247;230 245 208;184 225 134;127 188 65;77 146 33;77 146 33];
cDat.PRGn     = [64 0 75;118 42 131;153 112 171;194 165 207;231 212 232;247 247 247;217 240 211;166 219 160;90 174 97;27 120 55;27 120 55];
cDat.PuOr     = [127 59 8;179 88 6;224 130 20;253 184 99;254 224 182;247 247 247;216 218 235;178 171 210;128 115 172;84 39 136;84 39 136];
cDat.RdBu     = [103 0 31;178 24 43;214 96 77;244 165 130;253 219 199;247 247 247;209 229 240;146 197 222;67 147 195;33 102 172;33 102 172];
cDat.RdGy     = [103 0 31;178 24 43;214 96 77;244 165 130;253 219 199;255 255 255;224 224 224;186 186 186;135 135 135;77 77 77;77 77 77];
cDat.RdYlBu   = [165 0 38;215 48 39;244 109 67;253 174 97;254 224 144;255 255 191;224 243 248;171 217 233;116 173 209;69 117 180;69 117 180];
cDat.RdYlGn   = [165 0 38;215 48 39;244 109 67;253 174 97;254 224 139;255 255 191;217 239 139;166 217 106;102 189 99;26 152 80;26 152 80];
cDat.Spectral = [158 1 66;213 62 79;244 109 67;253 174 97;254 224 139;255 255 191;230 245 152;171 221 164;102 194 165;50 136 189;50 136 189];

% type of colormap
cTypeName = {'Qualitative' 'Sequential' 'Diverging'};
cType     = [ones(1,8) 2*ones(1,19) 3*ones(1,8)];

if nargin == 0
    help sw_cbrewer
end

cNames = fieldnames(cDat);

if nargin == 0 || strcmp(cName,'plot')
    
    hFig = figure('position', [314 327 807 420]);
    
    for ii = 1:3
        subplot(1,3,ii);
        % select color with the right type
        cNameSel = cNames(cType == ii);
        for jj = 1:numel(cNameSel)
            
            ncol = size(cDat.(cNameSel{jj}),1); 
            fg=1./ncol; % geometrical factor
            
            X=fg.*[0 0 1 1];
            Y=0.1.*[1 0 0 1]+(2*jj-1)*0.1;
            F=sw_cbrewer(cNameSel{jj}, ncol);
            
            for icol = 1:ncol
                X2 = X+fg.*(icol-1);
                fill(X2,Y,F(icol, :), 'linestyle', 'none')
                hold all
            end
            text(-0.1, mean(Y), cNameSel{jj}, 'HorizontalAlignment', 'right', 'FontSize', 14)
            xlim([-0.4, 1])
            
        end
        
        title(cTypeName{ii}, 'FontWeight', 'bold', 'FontSize', 16)
        axis off
        set(hFig,'color',[1 1 1])
        set(hFig, 'MenuBar', 'none')
        set(hFig, 'Name', 'sw_cbrewer() color maps')

    end
    text(-2.8,-0.1,'Press Enter to close figure or Ctrl+C to keep it!','fontSize',14)
    pause
    close(hFig)
    return
end

if ~any(ismember(cNames,cName))
    error('sw_cbrewer:WrongInput','The given color does not exists!')
end

if nargin < 2
    % use Matlab default number
    nCol = 64;
end

% number of colors integer
nCol = round(nCol);

cMap0 = cDat.(cName)/255;

% number of colors in the original colomap
nCol0 = size(cMap0,1);

if nCol0 >= nCol
    cMap = cMap0;
else
    interpMethod = 'pchip';
    
    % create the associated X axis (using round to get rid of decimals)
    X = linspace(1,nCol0,nCol);
    
    cMap = zeros(nCol,3);
    for ii = 1:3
        cMap(:,ii) = interp1(1:nCol0,cMap0(:,ii),X,interpMethod)';
    end
    
end

end