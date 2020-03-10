function [mdh,mask] = evalMDH( mdh_blob )
% see pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
% and pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h

if ~isa( mdh_blob, 'uint8' )
    error( [mfilename() ':NoInt8'], 'Binary mdh data must be a uint8 array!' )
end

mdh_blob = mdh_blob([1:20 41:end], :);     % remove 20 unnecessary bytes

Nmeas   = size( mdh_blob, 2 );

mdh.ulPackBit   = bitget( mdh_blob(4,:), 2).';
mdh.ulPCI_rx    = bitset(bitset(mdh_blob(4,:), 7, 0), 8, 0).'; % keep 6 relevant bits
mdh_blob(4,:)   = bitget( mdh_blob(4,:),1);  % ubit24: keep only 1 bit from the 4th byte

% unfortunately, typecast works on vectors, only
data_uint32     = typecast( reshape(mdh_blob(1:76,:),  [],1), 'uint32' );
data_uint16     = typecast( reshape(mdh_blob(29:end,:),[],1), 'uint16' );
data_single     = typecast( reshape(mdh_blob(69:end,:),[],1), 'single' );

data_uint32 = reshape( data_uint32, [], Nmeas ).';
data_uint16 = reshape( data_uint16, [], Nmeas ).';
data_single = reshape( data_single, [], Nmeas ).';
                                                        %  byte pos.
%mdh.ulDMALength               = data_uint32(:,1);      %   1 :   4
mdh.lMeasUID                   = data_uint32(:,2);      %   5 :   8
mdh.ulScanCounter              = data_uint32(:,3);      %   9 :  12
mdh.ulTimeStamp                = data_uint32(:,4);      %  13 :  16
mdh.ulPMUTimeStamp             = data_uint32(:,5);      %  17 :  20
mdh.aulEvalInfoMask            = data_uint32(:,6:7);    %  21 :  28
mdh.ushSamplesInScan           = data_uint16(:,1);      %  29 :  30
mdh.ushUsedChannels            = data_uint16(:,2);      %  31 :  32
mdh.sLC                        = data_uint16(:,3:16);   %  33 :  60
mdh.sCutOff                    = data_uint16(:,17:18);  %  61 :  64
mdh.ushKSpaceCentreColumn      = data_uint16(:,19);     %  66 :  66
mdh.ushCoilSelect              = data_uint16(:,20);     %  67 :  68
mdh.fReadOutOffcentre          = data_single(:, 1);     %  69 :  72
mdh.ulTimeSinceLastRF          = data_uint32(:,19);     %  73 :  76
mdh.ushKSpaceCentreLineNo      = data_uint16(:,25);     %  77 :  78
mdh.ushKSpaceCentrePartitionNo = data_uint16(:,26);     %  79 :  80

mdh.SlicePos                    = data_single(:, 4:10); %  81 : 108
mdh.aushIceProgramPara          = data_uint16(:,41:64); % 109 : 156
mdh.aushFreePara                = data_uint16(:,65:68); % 157 : 164

% inlining of evalInfoMask
evalInfoMask1 = mdh.aulEvalInfoMask(:,1);
mask.MDH_ACQEND             = min(bitand(evalInfoMask1, 2^0), 1);
mask.MDH_RTFEEDBACK         = min(bitand(evalInfoMask1, 2^1), 1);
mask.MDH_HPFEEDBACK         = min(bitand(evalInfoMask1, 2^2), 1);
mask.MDH_SYNCDATA           = min(bitand(evalInfoMask1, 2^5), 1);
mask.MDH_RAWDATACORRECTION  = min(bitand(evalInfoMask1, 2^10),1);
mask.MDH_REFPHASESTABSCAN   = min(bitand(evalInfoMask1, 2^14),1);
mask.MDH_PHASESTABSCAN      = min(bitand(evalInfoMask1, 2^15),1);
mask.MDH_SIGNREV            = min(bitand(evalInfoMask1, 2^17),1);
mask.MDH_PHASCOR            = min(bitand(evalInfoMask1, 2^21),1);
mask.MDH_PATREFSCAN         = min(bitand(evalInfoMask1, 2^22),1);
mask.MDH_PATREFANDIMASCAN   = min(bitand(evalInfoMask1, 2^23),1);
mask.MDH_REFLECT            = min(bitand(evalInfoMask1, 2^24),1);
mask.MDH_NOISEADJSCAN       = min(bitand(evalInfoMask1, 2^25),1);
mask.MDH_VOP                = min(bitand(mdh.aulEvalInfoMask(2), 2^(53-32)),1); % was 0 in VD

mask.MDH_IMASCAN            = ones( Nmeas, 1, 'uint32' );

noImaScan = (   mask.MDH_ACQEND             | mask.MDH_RTFEEDBACK   | mask.MDH_HPFEEDBACK       ...
              | mask.MDH_PHASCOR            | mask.MDH_NOISEADJSCAN | mask.MDH_PHASESTABSCAN    ...
              | mask.MDH_REFPHASESTABSCAN   | mask.MDH_SYNCDATA                                 ... 
              | (mask.MDH_PATREFSCAN & ~mask.MDH_PATREFANDIMASCAN) );

mask.MDH_IMASCAN( noImaScan ) = 0;

end % of evalMDH()