%%

swpref.setpref('tid',0,'fid',0)
model = sw_model('triAF',1);

N = 1e4;

tic
spec = model.spinwave(rand(3,N));
toc

tic
spec = model.spinwavefast(rand(3,N));
toc

tic
spmd
    swpref.setpref('tid',0,'fid',0)
    spec = model.spinwavefast(rand(3,N));
end
toc