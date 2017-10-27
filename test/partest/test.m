function test(fName)

if exist(fName,'file')
    delete(fName);
end

partest('nQ',1e5,'nWorker',[8 16 24 32],'nThread',4,'nRun',5,'fName',fName,'nSlice',1);
partest('nQ',1e5,'nWorker',[8 16 24 32],'nThread',4,'nRun',5,'fName',fName,'nSlice',4);
partest('nQ',1e5,'nWorker',[8 16 24 32],'nThread',8,'nRun',5,'fName',fName,'nSlice',1);
partest('nQ',1e5,'nWorker',[8 16 24 32],'nThread',8,'nRun',5,'fName',fName,'nSlice',4);

end