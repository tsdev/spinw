% run tests
delete('test.mat');
partest('nQ',1e5,'nWorker', 6,'nThread',4,'nRun',5,'fName','test.mat','nSlice',4);
partest('nQ',1e5,'nWorker',12,'nThread',4,'nRun',5,'fName','test.mat','nSlice',4);
partest('nQ',1e5,'nWorker',24,'nThread',4,'nRun',5,'fName','test.mat','nSlice',4);

partest('nQ',1e5,'nWorker', 6,'nThread',8,'nRun',5,'fName','test.mat','nSlice',4);
partest('nQ',1e5,'nWorker',12,'nThread',8,'nRun',5,'fName','test.mat','nSlice',4);
partest('nQ',1e5,'nWorker',24,'nThread',8,'nRun',5,'fName','test.mat','nSlice',4);

partest('nQ',1e5,'nWorker', 6,'nThread',4,'nRun',5,'fName','test.mat','nSlice',1);
partest('nQ',1e5,'nWorker',12,'nThread',4,'nRun',5,'fName','test.mat','nSlice',1);
partest('nQ',1e5,'nWorker',24,'nThread',4,'nRun',5,'fName','test.mat','nSlice',1);