% Encode_Ndataset(filename, TD)
% Encodes the Temporal Difference (TD) events into binary file for the
% neuromorphic datasets. It is the opposite of Read_Ndataset.m


% modified for 31 bit timestamp (up to 2147 seconds or ~30minutes)
function message = Encode_Ndataset_longTs(filename, TD)
%ByteVector = zeros(1,length(TD.ts)*5);
ByteVector = zeros(1,length(TD.ts)*6);
ByteVector(1:6:end) = TD.x; %TD.x-1;
ByteVector(2:6:end) = TD.y; %TD.y-1;
% ByteVector(3:5:end) = bitshift(2-TD.p,7) + double(bitand(bitshift(TD.ts,-16), 127));
%ByteVector(3:5:end) = bitshift(TD.p-1,7) + double(bitand(bitshift(TD.ts,-16), 127));
ByteVector(3:6:end) = bitshift(TD.p-1,7) + double(bitand(bitshift(TD.ts,-24), 127));
ByteVector(4:6:end) = bitand(bitshift(TD.ts,-16), 255);
ByteVector(5:6:end) = bitand(bitshift(TD.ts,-8), 255);
ByteVector(6:6:end) = bitand(TD.ts, 255);
[filename '.bs2']
[fileID, message] = fopen([filename '.bs2'], 'w');
if(fileID ~= -1)    
    fwrite(fileID, ByteVector);
end
fclose(fileID);
return