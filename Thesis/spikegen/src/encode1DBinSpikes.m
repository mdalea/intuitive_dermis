% function to encode spike data into binary event format represented by 40
% bit number
% the first 16 bits (bits 39-24) represent the neuronID
% bit 23 represents the sign of spike event: 0=>OFF event, 1=>ON event
% the last 23 bits (bits 22-0) represent the spike event timestamp in
% microseconds
% 
% Usage:
% encodeSpikes(<name fo spike data file>, <eventStamp = [neuronID, timeID(ms), polarity]>)
function encode1DBinSpikes(filename, eventStamp)
neuronID = uint16(eventStamp(:, 1)) - 1;
polarity = uint32(eventStamp(:, 3))*uint32(hex2dec('80'));
%timeID   = uint32(eventStamp(:, 2)*1000);
timeID   = uint32(eventStamp(:, 2));

eventStream = [bitshift(bitand(neuronID(:)', hex2dec('FF00')), -8)
                        bitand(neuronID(:)', hex2dec('00FF'))
               bitshift(bitand(timeID(:)'  , hex2dec('7F0000')), -16) + polarity(:)'
               bitshift(bitand(timeID(:)'  , hex2dec('FF00')), -8)
                        bitand(timeID(:)'  , hex2dec('FF'))
                        ];

eventStream = uint8(eventStream);

fid = fopen(filename, 'wb');
fwrite(fid, eventStream);
fclose(fid);