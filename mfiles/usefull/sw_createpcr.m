function sw_createpcr(path, pcrFile, perm)
% SW_CREATEPCR(path, pcrFile, perm) creates the structural part of a pcr
% file from a .cif file.
%
% This function will create the atomic positions from a .cif file.
% pcr file is the control file for FullProf Rietveld refinement software.
%
% perm  Permutation of the (x,y,z) coordinates.
%

if nargin == 0
    help sw_createpcr;
    return;
end

if nargin < 3
    perm = 1:3;
end

swdat  = sw(path);
cifdat = cif(path);

mult = cifdat.('atom_site_symmetry_multiplicity');
mult = mult/max(mult);


fid = fopen(pcrFile,'w');
fprintf(fid,'!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes\n');

nAtom = size(swdat.unit_cell.r,2);

uc = swdat.unit_cell;

idx = 1;
for ii = 1:nAtom
    if (ii>1) && strcmp(uc.label{ii},uc.label{ii-1})
        idx = idx + 1;
    else
        idx = 1;
    end
    
    if numel(uc.label{ii})>1
        fprintf(fid,'%s%d',uc.label{ii},idx);
        fprintf(fid,'    %s    ',uc.label{ii});
    else
        fprintf(fid,'%s%d',uc.label{ii},idx);
        fprintf(fid,'     %s     ',uc.label{ii});
        
    end
    
    fprintf(fid,'%9.5f%9.5f%9.5f%9.5f%9.5f%4d%4d%4d%4d\n',uc.r(perm,ii)',0,mult(ii),[0 0 0 0]);
    fprintf(fid,'                  0.00     0.00     0.00     0.00      0.00\n');
end

fclose(fid);

end