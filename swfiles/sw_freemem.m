function mem = sw_freemem()
% gives the amount of free RAM in bytes
%
% mem = SW_FREEMEM()
%
% If it cannot determine the size of the free memory, it returns zero.
%

mem = 0;

try %#ok<TRYNC>
    if ispc
        memStr = memory;
        mem = memStr.MemAvailableAllArrays;
    elseif ismac
        [~,memStr] = unix('vm_stat | grep free');
        mem = sscanf(memStr(14:end),'%f')*4096;
    elseif isunix
        % get the buffer/cache size
        [~, memStr] = unix('free -b | grep ''-''');
        
        if isempty(memStr)
            % there is no buffer/cache, just get the 'Mem' calues
            [~, memStr] = unix('free -b | grep ''Mem''');
            [~, mem_free] = strtok(memStr);
            mem = sscanf(mem_free,'%f');
            mem = mem(3);
        else
            [~, mem_free] = strtok(memStr(20:end));
            mem = num2str(mem_free);
        end
    end
end

end