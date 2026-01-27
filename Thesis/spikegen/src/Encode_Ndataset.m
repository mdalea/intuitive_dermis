% Encode_Ndataset(filename, TD)
% Encodes the Temporal Difference (TD) events into binary file for the
% neuromorphic datasets. It is the opposite of Read_Ndataset.m
function message = Encode_Ndataset(filename, TD)
ByteVector = zeros(1,length(TD.ts)*5);
ByteVector(1:5:end) = TD.x;%-1;
ByteVector(2:5:end) = TD.y;%-1;
% ByteVector(3:5:end) = bitshift(2-TD.p,7) + double(bitand(bitshift(TD.ts,-16), 127));
ByteVector(3:5:end) = bitshift(TD.p-1,7) + double(bitand(bitshift(TD.ts,-16), 127));
ByteVector(4:5:end) = bitand(bitshift(TD.ts,-8), 255);
ByteVector(5:5:end) = bitand(TD.ts, 255);

[fileID, message] = fopen([filename '.bs2'], 'w');
if(fileID ~= -1)    
    fwrite(fileID, ByteVector);
end
fclose(fileID);
return