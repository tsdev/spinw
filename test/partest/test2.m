function test2(fName,nQ)

if exist(fName,'file')
    delete(fName);
end

if nargin<2
    nQ = 1e5;
end

partest2('nQ',nQ,'nWorker',zeros(1,0),'nThread',4,'nRun',5,'fName','test01.mat','nSlice',1,'hermit',true)

end