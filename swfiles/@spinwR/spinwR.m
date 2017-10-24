classdef spinwR < handle
    %SPINWR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spinw_obj = [];
        baseURL = '';
        username = swpref.getpref('remoteuser').val
    end
    properties(SetAccess = protected)
        status = 'Waiting'
        statusURL = ''
        token = swpref.getpref('remotetoken').val
        token_expire = datetime('now')
    end
    properties(Hidden=true,SetAccess = protected)
        isCalculating = false;
        version = []
        deployed = false;
    end
    
    methods
        function obj = spinwR(sw)
            %SPINWR Construct an instance of this class
            %   Detailed explanation goes here
            obj.baseURL = swpref.getpref('remoteurl').val;
            if isempty(obj.baseURL)
                error('A valid server needs to  be specified.')
            end
            try
                url = strcat(obj.baseURL,'/spinw/version');
                obj.version = webread(url,weboptions('ContentType','json'));
            catch ME
                if strcmp(ME.identifier,'MATLAB:webservices:ExpectedProtocol')
                    error('A valid server needs to  be specified.')
                else
                    rethrow(ME)
                end
            end
            
            obj.spinw_obj = sw;
            obj.login();
        end
        
        
        function login(obj,varargin)
            % No token, so login,
            if isempty(obj.token)
                % Try make user. Catch existing user and get token.
                [obj.username, password] = obj.GetAuthentication(varargin{:});
                try
                    % Make a new user
                    obj.newUser(password)
                catch ME
                    % The user exists, get a token.
                    if strcmp(ME.identifier,'MATLAB:webservices:HTTP409StatusCodeError')
                        obj.status = 'User Exists';
                        obj.getToken(password)
                    else
                        % Something has gone wrong.
                        obj.token = '';
                        obj.status = ME.message;
                        return
                    end
                end
            else
                % A token exists but is it valid?
                url = strcat(obj.baseURL,'/users/quota');
                try
                    % If it's valid we get a new token.
                    [~] = webread(url,weboptions('Username',obj.token,'ContentType','json'));
                    obj.getToken(obj.token) % Re-issue a token.
                catch ME
                    % Catch the Unauthorized error and login.
                    if strcmp(ME.identifier,'MATLAB:webservices:HTTP401StatusCodeError')
                        obj.getToken()
                    else
                        obj.status = ME.message;
                        obj.token = '';
                        return
                    end
                end
            end
            swpref.setpref('remotetoken',obj.token)
            swpref.setpref('remoteuser',obj.username)
        end
        
        
        function newUser(obj,varargin)
            if isempty(varargin)
                password = [];
            elseif length(varargin) == 1
                password = varargin{1};
            elseif length(varargin) == 2
                obj.username = varargin{1};
                password = varargin{2};
            end
            if isempty(password)
                [obj.username, password] = obj.GetAuthentication();
                swpref.setpref('remoteuser',obj.username)
            end
            url = strcat(obj.baseURL,'/users');
            try
                response = webwrite(url,'username',obj.username,'password',password,weboptions('ContentType','json'));
            catch ME
                if strcmp(ME.identifier,'MATLAB:webservices:HTTP409StatusCodeError')
                    obj.status = 'User Exists';
                    response.username = obj.username;
                else
                    obj.status = ME.message;
                    return
                end
            end
            if strcmp(response.username,obj.username)
                url = strcat(obj.baseURL,'/users/token');
                temp = webread(url,weboptions('Username',obj.username,'Password',password,'ContentType','json'));
                obj.token = temp.token;
                obj.token_expire = datetime('now') + seconds(temp.duration);
            end
        end
        
        function quota = view_quota(obj)
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj = obj.getToken();
            end
            url = strcat(obj.baseURL,'/users/quota');
            quota = webread(url,weboptions('Username',obj.token,'ContentType','json'));
        end
        
        function getToken(obj,varargin)
            if isempty(varargin)
                password = [];
            elseif length(varargin) == 1
                password = varargin{1};
            elseif length(varargin) == 2
                obj.username = varargin{1};
                password = varargin{2};
            end
            if isempty(password)
                [obj.username, password] = obj.GetAuthentication();
            end
            url = strcat(obj.baseURL,'/users/token');
            try
                temp = webread(url,weboptions('Username',obj.username,'Password',password,'ContentType','json'));
                obj.token = temp.token;
                obj.token_expire = datetime('now') + seconds(599);
            catch ME
                % Catch the Unauthorized error and create new user.
                if strcmp(ME.identifier,'MATLAB:webservices:HTTP401StatusCodeError')
                    obj.token = '';
                    obj.newUser()
                else
                    obj.status = ME.message;
                    obj.token = '';
                    return
                end
            end
        end
        
        function upload(obj)
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj.getToken();
            end
            url = strcat(obj.baseURL,'/spinw/upload');
            filename = strcat(tempname,'.mat');
            %             sw_obj = obj.spinw_obj;
            d = char(getByteStreamFromArray(obj.spinw_obj));
            %             save(filename,'sw_obj')
            %             f = fopen(filename);
            %             d = char(fread(f)');
            %             fclose(f);
            
            [~,remoteFName, remoteExt] = fileparts(filename);
            opt = weboptions('Username',obj.token,'Password','x',...
                'characterEncoding','ISO-8859-1',...
                'MediaType','application/octet-stream',...
                'RequestMethod','post',...
                'HeaderFields',string({'Content-Length',string(length(d))}),...
                'ContentType','json');
            try
                upload_data = webwrite(sprintf(strcat(url,'/%s%s'),remoteFName,remoteExt), d, opt);
            catch someException
                throw(addCause(MException('uploadToSpinW:unableToUploadFile','Unable to upload file.'),someException));
            end
            obj.status = 'Uploaded File';
            obj.statusURL = upload_data.status;
        end
        
        function spinwave(obj,hkl,varargin)
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj = obj.getToken();
            end
            
            url = strcat(obj.baseURL,'/spinw/spinwave');
            filename = strcat(tempname,'.mat');
            if (obj.version.Deployed)
                sw_opt = struct();
                sw_opt.fun = 'spinwave';
                sw_opt.argin = {obj.spinw_obj,hkl,varargin{:}};
                sw_opt.nargout = 1;
            else
                sw_opt.Q = hkl;
                if ~isempty(varargin)
                    if ~mod(length(varargin),2)
                        for i = 1:2:length(varargin)
                            sw_opt.(varargin{i}) = varargin{i+1};
                        end
                    else
                        error('Uneven parameter/value pairs')
                    end
                end
            end
            %             save(filename,'sw_opt')
            %             f = fopen(filename);
            %             d = char(fread(f)');
            %             fclose(f);
            d = char(getByteStreamFromArray(sw_opt));
            
            [~,remoteFName, remoteExt] = fileparts(filename);
            opt = weboptions('Username',obj.token,'Password','x',...
                'characterEncoding','ISO-8859-1',...
                'MediaType','application/octet-stream',...
                'RequestMethod','post',...
                'HeaderFields',string({'Content-Length',string(length(d))}),...
                'ContentType','json');
            try
                tempOutput = webwrite(sprintf(strcat(url,'/%s%s'),remoteFName,remoteExt), d, opt);
            catch someException
                throw(addCause(MException('uploadToSpinW:unableToUploadFile','Unable to upload file.'),someException));
            end
            if obj.version.Deployed
                obj.statusURL = tempOutput.status;
            end
            if tempOutput.Calculating && ~tempOutput.Errors
                iscomp = webread(obj.statusURL,weboptions('Username',obj.token,'ContentType','json','timeout',100));
                if strcmp(iscomp.status,'running')
                    obj.status = 'Calculating';
                    obj.isCalculating = true;
                end
            end
        end
        
        function powspec(obj,hkl,varargin)
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj = obj.getToken();
            end
            url = strcat(obj.baseURL,'/spinw/powspec');
            if (obj.version.Deployed)
                sw_opt = struct();
                sw_opt.fun = 'powspec';
                sw_opt.argin = {obj.spinw_obj,hkl,varargin{:}};
                sw_opt.nargout = 1;
            else
                sw_opt.Q = hkl;
                if ~isempty(varargin)
                    if ~mod(length(varargin),2)
                        for i = 1:2:length(varargin)
                            sw_opt.(varargin{i}) = varargin{i+1};
                        end
                    else
                        error('Uneven parameter/value pairs')
                    end
                end
            end
            filename = strcat(tempname,'.mat');
            d = char(getByteStreamFromArray(sw_opt));
            %             save(filename,'sw_opt')
            %             f = fopen(filename);
            %             d = char(fread(f)');
            %             fclose(f);
            
            [~,remoteFName, remoteExt] = fileparts(filename);
            opt = weboptions('Username',obj.token,'Password','x',...
                'characterEncoding','ISO-8859-1',...
                'MediaType','application/octet-stream',...
                'RequestMethod','post',...
                'HeaderFields',string({'Content-Length',string(length(d))}),...
                'ContentType','json');
            try
                tempOutput = webwrite(sprintf(strcat(url,'/%s%s'),remoteFName,remoteExt), d, opt);
            catch someException
                throw(addCause(MException('uploadToSpinW:unableToUploadFile','Unable to upload file.'),someException));
            end
            if tempOutput.Calculating && ~tempOutput.Errors
                iscomp = webread(obj.statusURL,weboptions('Username',obj.token,'ContentType','json','timeout',100));
                if strcmp(iscomp.status,'running')
                    obj.status = 'Calculating';
                    obj.isCalculating = true;
                end
            end
        end
        
        
        function spectra = getResult(obj,varargin)
            spectra = [];
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj = obj.getToken();
            end
            if isempty(varargin)
                timeout = 1;
            else
                timeout = varargin{1};
            end
            iscomp = webread(obj.statusURL,weboptions('Username',obj.token,'ContentType','json','timeout',100));
            if strcmp(iscomp.status,'done')
                obj.isCalculating = false;
                obj.status = 'Calculation complete';
            end
            if obj.isCalculating
                iscomp = webread(obj.statusURL,weboptions('Username',obj.token,'ContentType','json','timeout',100));
                if strcmp(iscomp.status,'running')
                    obj.status = 'Calculating';
                    cont = true;
                    while cont
                        iscomp = webread(obj.statusURL,weboptions('Username',obj.token,'ContentType','json','timeout',100));
                        if strcmp(iscomp.status,'done')
                            cont = false;
                            obj.isCalculating = false;
                            obj.status = 'Calculation complete';
                        else
                            pause(timeout)
                        end
                    end
                end
            end
            filename = tempname;
            file = websave(filename,iscomp.url,weboptions('Username',obj.token));
            if obj.version.Deployed
                temp = load(file);
                if ~isempty(temp.err)
                    rethrow(temp.err)
                end
                spectra = temp.argout{:};
                disp(temp.log)
            else
                load(file,'-mat','spectra')
            end
        end
        
        function evaluate(obj,fn_name,varargin)
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj = obj.getToken();
            end
            url = strcat(obj.baseURL,'/spinw/compute');
            if (obj.version.Deployed)
                sw_opt = struct();
                sw_opt.fun = fn_name;
                sw_opt.argin = {obj.spinw_obj,varargin{:}};
                sw_opt.nargout = 1;
            else
                error('This command only works on a compiled server.')
            end
            filename = strcat(tempname,'.mat');
            d = char(getByteStreamFromArray(sw_opt));
            
            [~,remoteFName, remoteExt] = fileparts(filename);
            opt = weboptions('Username',obj.token,'Password','x',...
                'characterEncoding','ISO-8859-1',...
                'MediaType','application/octet-stream',...
                'RequestMethod','post',...
                'HeaderFields',string({'Content-Length',string(length(d))}),...
                'ContentType','json');
            try
                tempOutput = webwrite(sprintf(strcat(url,'/%s%s'),remoteFName,remoteExt), d, opt);
            catch someException
                throw(addCause(MException('uploadToSpinW:unableToUploadFile','Unable to upload file.'),someException));
            end
            if tempOutput.Calculating && ~tempOutput.Errors
                iscomp = webread(obj.statusURL,weboptions('Username',obj.token,'ContentType','json','timeout',100));
                if strcmp(iscomp.status,'running')
                    obj.status = 'Calculating';
                    obj.isCalculating = true;
                end
            end
        end
        
        
        function [username,password]=GetAuthentication(obj)
            %GetAuthentication prompts a username and password from a user and hides the
            % password input by *****
            %
            %   [user,password] = GetAuthentication;
            %   [user,password] = GetAuthentication(defaultuser);
            %
            % arguments:
            %   defaultuser - string for default name
            %
            % results:
            %   username - string for the username
            %   password - password as a string
            %
            % Created by Felix Ruhnow, MPI-CBG Dresden
            % Version 1.00 - 20th February 2009
            %
            
            defaultuser = obj.username;
            
            hAuth.fig = figure('Menubar','none','Units','normalized','Resize','off','NumberTitle','off', ...
                'Name','Authentication','Position',[0.4 0.4 0.2 0.2],'WindowStyle','normal');
            
            uicontrol('Parent',hAuth.fig,'Style','text','Enable','inactive','Units','normalized','Position',[0 0 1 1], ...
                'FontSize',12);
            
            uicontrol('Parent',hAuth.fig,'Style','text','Enable','inactive','Units','normalized','Position',[0.1 0.8 0.8 0.1], ...
                'FontSize',12,'String','Username:','HorizontalAlignment','left');
            
            
            hAuth.eUsername = uicontrol('Parent',hAuth.fig,'Style','edit','Tag','username','Units','normalized','Position',[0.1 0.675 0.8 0.125], ...
                'FontSize',12,'String',defaultuser,'BackGroundColor','white','HorizontalAlignment','left');
            
            uicontrol('Parent',hAuth.fig,'Style','text','Enable','inactive','Units','normalized','Position',[0.1 0.5 0.8 0.1], ...
                'FontSize',12,'String','Password:','HorizontalAlignment','left');
            
            hAuth.ePassword = uicontrol('Parent',hAuth.fig,'Style','edit','Tag','password','Units','normalized','Position',[0.1 0.375 0.8 0.125], ...
                'FontSize',12,'String','','BackGroundColor','white','HorizontalAlignment','left');
            
            uicontrol('Parent',hAuth.fig,'Style','pushbutton','Tag','OK','Units','normalized','Position',[0.1 0.05 0.35 0.2], ...
                'FontSize',12,'String','OK','Callback','uiresume;');
            
            uicontrol('Parent',hAuth.fig,'Style','pushbutton','Tag','Cancel','Units','normalized','Position',[0.55 0.05 0.35 0.2], ...
                'FontSize',12,'String','Cancel','Callback',@AbortAuthentication);
            
            set(hAuth.fig,'CloseRequestFcn',@AbortAuthentication)
            set(hAuth.ePassword,'KeypressFcn',@PasswordKeyPress)
            
            setappdata(0,'hAuth',hAuth);
            uicontrol(hAuth.eUsername);
            uiwait;
            
            username = get(hAuth.eUsername,'String');
            password = get(hAuth.ePassword,'UserData');
            delete(hAuth.fig);
            
            function PasswordKeyPress(hObject,event)
                hAuth = getappdata(0,'hAuth');
                password = get(hAuth.ePassword,'UserData');
                switch event.Key
                    case 'backspace'
                        password = password(1:end-1);
                    case 'return'
                        uiresume;
                        return;
                    otherwise
                        password = [password event.Character];
                end
                set(hAuth.ePassword,'UserData',password)
                set(hAuth.ePassword,'String',char('*'*sign(password)))
            end
            
            function AbortAuthentication(hObject,event)
                hAuth = getappdata(0,'hAuth');
                set(hAuth.eUsername,'String','');
                set(hAuth.ePassword,'UserData','');
                uiresume;
                
            end
        end
    end
end
