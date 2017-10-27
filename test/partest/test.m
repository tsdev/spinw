function test(fName,nQ)

if exist(fName,'file')
    delete(fName);
end

if nargin<2
    nQ = 1e5;
end

partest('nQ',nQ,'nWorker',[8 16 24 32],'nThread',4,'nRun',5,'fName',fName,'nSlice',1,'hermit',true);
partest('nQ',nQ,'nWorker',[8 16 24 32],'nThread',4,'nRun',5,'fName',fName,'nSlice',4,'hermit',true);
partest('nQ',nQ,'nWorker',[8 16 24 32],'nThread',8,'nRun',5,'fName',fName,'nSlice',1,'hermit',true);
partest('nQ',nQ,'nWorker',[8 16 24 32],'nThread',8,'nRun',5,'fName',fName,'nSlice',4,'hermit',true);

end