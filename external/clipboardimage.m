function rgbData = clipboardimage()
% returns the RGB image from clipboard
%
% rgbData = CLIPBOARDIMAGE()
%
% The function at the moment works with print screen images on OXS and
% Windows.
%
% Output:
%
% rgbData   Matrix with dimensions W x H x 3, containing the RGB values
%           from the clipboard. Empty matrix if the clipboard doesn't
%           contain image data.
%
% Example:
%
% figure
% imshow(clipboardimage)
%

% Based on code from Saurabh Kumar, saurabhkumar_@rediffmailcom
%

tKit    = java.awt.Toolkit.getDefaultToolkit();
% get clipboard handle
cbrd    = tKit.getSystemClipboard();
reqObj  = java.lang.Object;
img     = cbrd.getContents(reqObj);
Dflavor = img.getTransferDataFlavors();
imgDfvr = java.awt.datatransfer.DataFlavor.imageFlavor;

% only import if the clipboard stores image data
if Dflavor(1).equals(imgDfvr)
    % got image
    jImage = img.getTransferData(java.awt.datatransfer.DataFlavor.imageFlavor);
    % dimensions of the image
    imSize = [jImage.getHeight jImage.getWidth];
    % convert the javs image to Matlab matrix
    cmykDat = double(permute(reshape(typecast(jImage.getData.getDataStorage, 'uint8'), 4, imSize(2), imSize(1)),[3 2 1]))/255;
    % convert CMYK to RGB the transposed image
    rgbData = bsxfun(@times,cmykDat(:,:,[3 2 1]),cmykDat(:,:,4));
else
    % return empty matrix if there is no image
    rgbData = [];
end