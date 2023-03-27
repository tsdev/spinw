out_dir = 'ctf';
VERSION = version('-release');
package_name = ['SpinW_', VERSION];
full_package = ['SpinW_', VERSION, '.ctf'];

opts = compiler.build.ProductionServerArchiveOptions( ...
    ['matlab', filesep, 'call.m'], ...
    'ArchiveName', package_name, ...
    'OutputDir', 'ctf', ...
    'AutoDetectDataFiles', 'on', ...
    'AdditionalFiles', { ...
    ['..', filesep, 'swfiles'], ...
    ['..', filesep, 'external'], ...
    ['..', filesep, 'dat_files']});

compiler.build.productionServerArchive(opts);
movefile([out_dir, filesep, full_package] , ['pyspinw', filesep, 'ctfs']);


version_file = ['pyspinw', filesep, 'ctfs', filesep 'versions'];
data = {};
write_data = true;
if exist(version_file,"file")
    data_str = fileread(version_file);
    if ~isempty(data_str)
        data = jsondecode(data_str);
        for in1 = 1:length(data)
            if strcmp(data(in1).version, VERSION) && strcmp(data(in1).file, full_package)
                write_data = false;
                break
            end
        end
    end
end

if write_data
    fid = fopen(['pyspinw', filesep, 'ctfs', filesep 'versions'], 'w');
    data{end+1} = struct('version', VERSION, 'file', full_package);
    fwrite(fid, jsonencode(data));
    fclose(fid);
end
