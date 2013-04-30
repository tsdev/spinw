function mem = sw_freemem()
% mem = SW_FREEMEM() gives free memory in bytes.
%

if ispc
    memStr = memory;
    mem = memStr.MemAvailableAllArrays;
elseif ismac || isunix
    [~,memStr] = unix('vm_stat | grep free');
    mem = sscanf(memStr(14:end),'%f')*4096;
end

end