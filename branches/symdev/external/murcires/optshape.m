function sh=optshape(d,shape)
% Calculates standard deviation for different distribution functions
%
% sh = OPTSHAPE(d,shape)
%
% Input:
%
% shape     Shape (distribution function):
%               'c'     circle, d = diameter,
%               'b'     box, d = width,
%               'd'     discrete, d = number of elements, assuming unit 
%                       distance between elements,
%               'g'     gauss function, d = FWHM,
%               Nx2 vector of (x,y) points
%
% Output:
%
% sh is struct typoe with the following fields:
% d         equal to input d,
% sig       standard deviation,
% FWHM      full width at half maximum,
% ix        inline function of the distribution
%               ix(deltaD)  where deltaD is the step size
%


s2f=(8*log(2))^0.5;
if max(size(shape))<3
    tmpsh=shape;
else
    tmpsh='f';
end;
%if d==0
%    d=max(shape(:,1))-min(shape(:,1));
%    tmpsh='ff';
%end;

sh.d=d;

switch tmpsh
    case 'b'
        sh.d=d;
        sh.sig=d/12^0.5;
        sh.FWHM=sh.sig*s2f;
        sh.ix=(@(res) [(-d/2:res:d/2)' (-d/2:res:d/2)'*0+1]);
    case 'd'
        sh.d=d;
        sh.sig=d/12^0.5;
        sh.FWHM=sh.sig*s2f;
        dh=(d-1)/2;
        tmpiis=[(-dh:1:dh)' ones(d,1)];
        sh.sig=mean(tmpiis(:,1).^2).^0.5;
        sh.ix=(@(res) tmpiis);
        
    case 'c'
        sh.d=d;
        sh.sig=d/4;
        sh.FWHM=sh.sig*s2f;
        sh.ix=(@(res) [(-d/2:res:d/2)' (d^2/4-(-d/2:res:d/2)'.^2).^0.5]);
    case 'g'
        sh.d=d;
        sh.sig=d;
        sh.FWHM=sh.sig*s2f;
        sh.ix=(@(res) [(-d*3:res:d*3)' exp(-1/2/d^2*(-d*3:6*d/100:d*3)'.^2)]);
    case 'gf'
        sh.d=d;
        sh.sig=d/s2f;
        sh.FWHM=d;
        d=sh.sig;
        sh.ix=(@(res) [(-d*3:res:d*3)' exp(-1/2/d^2*(-d*3:6*d/100:d*3)'.^2)]);
    case 'f'
        sh.ix=(@(res) [(min(shape(:,1)):res:max(shape(:,1)))' interp1(shape(:,1),shape(:,2),(min(shape(:,1)):res:max(shape(:,1)))')]);
        sh.sig=sum(shape(:,1)'.^2*shape(:,2)/sum(shape(:,2)))^0.5;
        sh.FWHM=sh.sig*s2f;
        sh.d=sh.FWHM;
    case 'ff'
        sh.ix=(@(res) shape);
        sh.sig=sum(shape(:,1)'.^2*shape(:,2)/sum(shape(:,2)))^0.5;
        sh.FWHM=sh.sig*s2f;
        sh.d=sh.FWHM;
end

end