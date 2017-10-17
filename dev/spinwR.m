classdef spinwR < handle
    %SPINWR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        spinw_obj = [];
        baseURL = 'http://127.0.0.1:5000'
        username = '';
    end
    properties(SetAccess = protected)
        status = 'Waiting'
        statusURL = ''
        token = ''
        token_expire = datetime('now')
    end
    properties(Hidden=true)
        isCalculating = false;
    end
    
    methods
        function obj = spinwR(sw,username)
            %SPINWR Construct an instance of this class
            %   Detailed explanation goes here
            obj.spinw_obj = sw;
            obj.username = username;
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
                obj.token_expire = datetime('now') + seconds(599);
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
            temp = webread(url,weboptions('Username',obj.username,'Password',password,'ContentType','json'));
            obj.token = temp.token;
            obj.token_expire = datetime('now') + seconds(599);
        end
        
        function upload(obj)
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj = obj.getToken();
            end
            url = strcat(obj.baseURL,'/spinw/upload');
            filename = strcat(tempname,'.mat');
            sw_obj = obj.spinw_obj;
            save(filename,'sw_obj')
            f = fopen(filename);
            d = char(fread(f)');
            fclose(f);
            
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
            filename = strcat(tempname,'.mat');
            save(filename,'sw_opt')
            f = fopen(filename);
            d = char(fread(f)');
            fclose(f);
            
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
        
        function powspec(obj,hkl,varargin)
            if isempty(obj.token)
                error('You need to login')
            end
            if (obj.token_expire - datetime('now')) < 0
                obj = obj.getToken();
            end
            url = strcat(obj.baseURL,'/spinw/powspec');
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
            filename = strcat(tempname,'.mat');
            save(filename,'sw_opt')
            f = fopen(filename);
            d = char(fread(f)');
            fclose(f);
            
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
            load(file,'-mat','spectra')
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
                
                function AbortAuthentication(hObject,event)
                    hAuth = getappdata(0,'hAuth');
                    set(hAuth.eUsername,'String','');
                    set(hAuth.ePassword,'UserData','');
                    uiresume;
                    
                end
            end
        end
    end
end
