classdef Twix < handle
    % Constructs Siemens TWIX object.
    % Properties and methods specific to CISS sequence development for
    % arbitrary sampling patterns and PICS reconstructions.
    
    properties (Access = public)
        measuredLineIndices % index for phase encoding
        measuredPartitionIndices % index for slice encoding
        numberOfAverages
        numberOfReadouts
        numberOfLines
        numberOfPartitions
        numberOfSets
        numberOfAcquisitions % acquisitions means total # of lines measured
        numberOfFakeAcquisitions % includes both measured & synthesized lines
        numberOfChannels
        kspaceData % raw k-space data
        reconstructedKSpaceDataToWrite % reconstructed raw k-space
        channelImageDataEcho1 % original images of each channel (first trufi acquisition)
        reconstructedChannelImageDataEcho1 % recon images of each channel
        sosImageEcho1 % sum of squares image of original data
        reconstructedSosImageEcho1 % reconstructed sum of squares image
        channelImageDataEcho2 % for the second trufi acquisition
        reconstructedChannelImageDataEcho2
        sosImageEcho2
        reconstructedSosImageEcho2
        fakeMeasurementLinesToWrite
        samplingPattern
        mdhScanHeader
        twixObject
        filePath
        numberOfMdh
        dataByteLength
        byteLengthOfSingleMdhLine
        numberOfScans
        byteOffsetToHeader
        byteOffsetToMeasurement
        byteLengthOfHeader
        byteLengthOfMeasurement
        byteDifference
        MID
        FID
        patientName
        protocolName
        totalFileBytes
        acqEndMeasurement
        fakeAcqEndMeasurement
    end
    
    properties (Access = private)
        numberOfRaidFiles = 64
        mdhScanHeaderByteLength = 192
        channelHeaderByteLength = 32
        numberOfNumberSlots = 2
        numberByteLength = 4
    end
    
    methods
        function twix = Twix(path, averageFlag)
            % constructs a twix object from the meas.dat file specified

            twix.filePath = path;
            twix.twixObject = mapVD(twix.filePath);
            twix.numberOfMdh = length(twix.twixObject);
            if twix.numberOfMdh > 1
                twix.twixObject = twix.twixObject{end};
            end
            twix.numberOfAverages = twix.twixObject.image.NAve;
            twix.numberOfReadouts = twix.twixObject.image.NCol;
            twix.numberOfLines = twix.twixObject.image.NLin;
            twix.numberOfPartitions = twix.twixObject.image.NPar;
            twix.numberOfSets = twix.twixObject.image.NSet;
            twix.numberOfAcquisitions = twix.twixObject.image.NAcq;
            twix.numberOfChannels = twix.twixObject.image.NCha;
            twix.measuredLineIndices =...
                twix.twixObject.image.Lin(1:twix.numberOfAcquisitions/twix.numberOfSets);
            twix.measuredPartitionIndices =...
                twix.twixObject.image.Par(1:twix.numberOfAcquisitions/twix.numberOfSets);
            if averageFlag
                twix.twixObject.image.flagDoAverage = 1;
            end
            twix.kspaceData = squeeze(twix.twixObject.image());
            twix.dataByteLength = twix.numberOfNumberSlots*...
                                  twix.numberByteLength*...
                                  twix.numberOfReadouts;
            twix.byteLengthOfSingleMdhLine =...
                twix.mdhScanHeaderByteLength...
               +twix.numberOfChannels*(twix.channelHeaderByteLength...
                                      +twix.dataByteLength);
                                  
            fid = fopen(twix.filePath,'r','l','US-ASCII');
            fseek(fid,0,'eof');
            twix.totalFileBytes = ftell(fid);
            fseek(fid,0,'bof');
            fread(fid,1,'uint32');
            twix.numberOfScans = fread(fid,1,'uint32');
            for i=1:twix.numberOfRaidFiles
                twix.MID(i) = fread(fid,1,'uint32');
                twix.FID(i) = fread(fid,1,'uint32');
                twix.byteOffsetToHeader(i) = fread(fid,1,'uint64');
                twix.byteLengthOfMeasurement(i) = fread(fid,1,'uint64');
                twix.patientName{i} = fread(fid,64,'uint8=>char');
                twix.protocolName{i} = fread(fid,64,'uint8=>char');
            end
            for i=1:twix.numberOfRaidFiles
                fseek(fid, twix.byteOffsetToHeader(i), 'bof');
                twix.byteLengthOfHeader(i) = fread(fid,1,'uint32');
                twix.byteOffsetToMeasurement(i) = twix.byteOffsetToHeader(i)...
                                                 +twix.byteLengthOfHeader(i);
            end
            fclose(fid);
        end
        
        function permuteKSpaceData(twix, permutation)
            % permute kspace data with given permutation
            
            twix.kspaceData = permute(twix.kspaceData, permutation);
        end
        
        function permuteReconstructedKSpaceData(twix, permutation)
            % permute reconstructed kspace data with given permutation
            
            twix.reconstructedKSpaceDataToWrite =...
                permute(twix.reconstructedKSpaceDataToWrite, permutation);
        end
        
        function sz = getKSpaceDataSize(twix)
            % outputs size of kspace data
            
            sz = size(twix.kspaceData);
        end
        
        function sz = getReconstructedKSpaceDataSize(twix)
            % outputs size of reconstructed kspace data
            
            sz = size(twix.reconstructedKSpaceDataToWrite);
        end
        
        function repeatReconstructedData(twix, repeats)
            % outputs size of reconstructed kspace data
            
            twix.reconstructedKSpaceDataToWrite =...
                repmat(twix.reconstructedKSpaceDataToWrite,...
                       [ones(1,ndims(twix.reconstructedKSpaceDataToWrite)) repeats]);
        end
        
        function extractSamplingPattern(twix)
            % computes the sampling mask from the measured line and
            % partition indices
            
            twix.samplingPattern = zeros(twix.numberOfLines,...
                                         twix.numberOfPartitions);
            for i=1:twix.numberOfAcquisitions/twix.numberOfSets
                twix.samplingPattern(twix.measuredLineIndices(i),...
                                     twix.measuredPartitionIndices(i)) = 1;
            end
        end
        
        function extractHeaderInfo(twix)
            fid = fopen(twix.filePath,'r','l','US-ASCII');
%             fseek(fid,0,'eof');
%             twix.totalFileBytes = ftell(fid);
%             fseek(fid,0,'bof');
%             fread(fid,1,'uint32');
%             twix.numberOfScans = fread(fid,1,'uint32');
%             for i=1:twix.numberOfRaidFiles
%                 twix.MID(i) = fread(fid,1,'uint32');
%                 twix.FID(i) = fread(fid,1,'uint32');
%                 twix.byteOffsetToHeader(i) = fread(fid,1,'uint64');
%                 twix.byteLengthOfMeasurement(i) = fread(fid,1,'uint64');
%                 twix.patientName{i} = fread(fid,64,'uint8=>char');
%                 twix.protocolName{i} = fread(fid,64,'uint8=>char');
%             end
%             for i=1:twix.numberOfRaidFiles
%                 fseek(fid, twix.byteOffsetToHeader(i), 'bof');
%                 twix.byteLengthOfHeader(i) = fread(fid,1,'uint32');
%                 twix.byteOffsetToMeasurement(i) = twix.byteOffsetToHeader(i)...
%                                                  +twix.byteLengthOfHeader(i);
%             end
            twix.byteDifference = twix.byteLengthOfMeasurement(twix.numberOfScans)...
                -twix.numberOfAcquisitions*twix.byteLengthOfSingleMdhLine;

            fseek(fid, twix.byteOffsetToMeasurement(twix.numberOfScans), 'bof');
            for i=1:twix.numberOfAcquisitions
                twix.mdhScanHeader(i).ulFlagsAndDMALength = fread(fid, 1, 'uint32');
                twix.mdhScanHeader(i).lMeasUID = fread(fid, 1, 'int32');
                twix.mdhScanHeader(i).ulScanCounter = fread(fid, 1, 'uint32');
                twix.mdhScanHeader(i).ulTimeStamp = fread(fid, 1, 'uint32');
                twix.mdhScanHeader(i).ulPMUTimeStamp = fread(fid, 1, 'uint32');
                twix.mdhScanHeader(i).ushSystemType = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).ulPTABPosDelay = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).lPTABPosX = fread(fid, 1, 'int32');
                twix.mdhScanHeader(i).lPTABPosY = fread(fid, 1, 'int32');
                twix.mdhScanHeader(i).lPTABPosZ = fread(fid, 1, 'int32');
                twix.mdhScanHeader(i).ulReserved1 = fread(fid, 1, 'uint32');
                twix.mdhScanHeader(i).aulEvalInfoMask = fread(fid, 2, 'uint32');
                twix.mdhScanHeader(i).ushSamplesInScan = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).ushUsedChannels = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).sLC = fread(fid, 14, 'uint16');
                twix.mdhScanHeader(i).sCutOff = fread(fid, 2, 'uint16');
                twix.mdhScanHeader(i).ushKSpaceCentreColumn = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).ushCoilSelect = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).fReadOutOffcentre = fread(fid, 1, 'float');
                twix.mdhScanHeader(i).ulTimeSinceLastRF = fread(fid, 1, 'uint32');
                twix.mdhScanHeader(i).ushKSpaceCentreLineNo = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).ushKSpaceCentrePartitionNo = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).sSD = fread(fid, 14, 'uint16');
                twix.mdhScanHeader(i).aushIceProgramPara = fread(fid, 24, 'uint16');
                twix.mdhScanHeader(i).aushReservedPara = fread(fid, 4, 'uint16');
                twix.mdhScanHeader(i).ushApplicationCounter = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).ushApplicationMask = fread(fid, 1, 'uint16');
                twix.mdhScanHeader(i).ulCRC = fread(fid, 1, 'uint32');
                for j=1:twix.numberOfChannels
                    twix.mdhScanHeader(i).channelHeader(j).ulTypeAndChannelLength = fread(fid, 1, 'uint32');
                    twix.mdhScanHeader(i).channelHeader(j).lMeasUID = fread(fid, 1, 'int32');
                    twix.mdhScanHeader(i).channelHeader(j).ulScanCounter = fread(fid, 1, 'uint32');
                    twix.mdhScanHeader(i).channelHeader(j).ulReserved1 = fread(fid, 1, 'uint32');
                    twix.mdhScanHeader(i).channelHeader(j).ulSequenceTime = fread(fid, 1, 'uint32');
                    twix.mdhScanHeader(i).channelHeader(j).ulUnused2 = fread(fid, 1, 'uint32');
                    twix.mdhScanHeader(i).channelHeader(j).ulChannelId = fread(fid, 1, 'uint16');
                    twix.mdhScanHeader(i).channelHeader(j).ulUnused3 = fread(fid, 1, 'uint16');
                    twix.mdhScanHeader(i).channelHeader(j).ulCRC = fread(fid, 1, 'uint32');
                    measurementData = fread(fid, twix.numberOfNumberSlots*twix.numberOfReadouts, 'float');
                    twix.mdhScanHeader(i).channelHeader(j).readoutLine =...
                        complex(measurementData(1:2:end), measurementData(2:2:end));
                end
            end
            fseek(fid,twix.byteOffsetToMeasurement(twix.numberOfScans)+...
                      twix.numberOfAcquisitions*twix.byteLengthOfSingleMdhLine,'bof');
            twix.acqEndMeasurement = fread(fid,twix.byteDifference,'uint8=>uint8');
            fseek(fid,twix.byteOffsetToMeasurement(twix.numberOfScans)+...
                      twix.numberOfAcquisitions*twix.byteLengthOfSingleMdhLine,'bof');
            twix.fakeAcqEndMeasurement.ulFlagsAndDMALength = fread(fid, 1, 'uint32');
            twix.fakeAcqEndMeasurement.lMeasUID = fread(fid, 1, 'int32');
%             twix.fakeAcqEndMeasurement.ulScanCounter = fread(fid, 1, 'uint32');
%             twix.fakeAcqEndMeasurement.ulTimeStamp = fread(fid, 1, 'uint32');
%             twix.fakeAcqEndMeasurement.ulPMUTimeStamp = fread(fid, 1, 'uint32');
            twix.fakeAcqEndMeasurement.ushSystemType = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.ulPTABPosDelay = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.lPTABPosX = fread(fid, 1, 'int32');
            twix.fakeAcqEndMeasurement.lPTABPosY = fread(fid, 1, 'int32');
            twix.fakeAcqEndMeasurement.lPTABPosZ = fread(fid, 1, 'int32');
            twix.fakeAcqEndMeasurement.ulReserved1 = fread(fid, 1, 'uint32');
            twix.fakeAcqEndMeasurement.aulEvalInfoMask = fread(fid, 2, 'uint32');
            twix.fakeAcqEndMeasurement.ushSamplesInScan = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.ushUsedChannels = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.sLC = fread(fid, 14, 'uint16');
            twix.fakeAcqEndMeasurement.sCutOff = fread(fid, 2, 'uint16');
            twix.fakeAcqEndMeasurement.ushKSpaceCentreColumn = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.ushCoilSelect = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.fReadOutOffcentre = fread(fid, 1, 'float');
            twix.fakeAcqEndMeasurement.ulTimeSinceLastRF = fread(fid, 1, 'uint32');
            twix.fakeAcqEndMeasurement.ushKSpaceCentreLineNo = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.ushKSpaceCentrePartitionNo = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.sSD = fread(fid, 14, 'uint16');
            twix.fakeAcqEndMeasurement.aushIceProgramPara = fread(fid, 24, 'uint16');
            twix.fakeAcqEndMeasurement.aushReservedPara = fread(fid, 4, 'uint16');
            twix.fakeAcqEndMeasurement.ushApplicationCounter = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.ushApplicationMask = fread(fid, 1, 'uint16');
            twix.fakeAcqEndMeasurement.ulCRC = fread(fid, 1, 'uint32');
            fclose(fid);
        end
        
        function createArtificalScanLinesToWrite(twix, measUID)
            twix.numberOfFakeAcquisitions = twix.numberOfSets...
                                           *twix.numberOfLines...
                                           *twix.numberOfPartitions...
                                           *twix.numberOfAverages;
            initialTimeStamp = twix.mdhScanHeader(1).ulTimeStamp;
            timeToAdd = twix.mdhScanHeader(2).ulTimeStamp-twix.mdhScanHeader(1).ulTimeStamp;
            initialPmuTimeStamp = twix.mdhScanHeader(1).ulPMUTimeStamp;
            pmuTimeToAdd = twix.mdhScanHeader(2).ulPMUTimeStamp-twix.mdhScanHeader(1).ulPMUTimeStamp;
            initialSequenceTimeStamp = twix.mdhScanHeader(1).channelHeader(1).ulSequenceTime;
            sequenceTimeToAdd = twix.mdhScanHeader(2).channelHeader(1).ulSequenceTime...
                               -twix.mdhScanHeader(1).channelHeader(1).ulSequenceTime;
            timeStampsToWrite = zeros(twix.numberOfFakeAcquisitions+1,1);
            pmuTimeStampsToWrite = zeros(twix.numberOfFakeAcquisitions+1,1);
            sequenceTimeStampsToWrite = zeros(twix.numberOfFakeAcquisitions+1,1);
            channelIds = [twix.mdhScanHeader(1).channelHeader.ulChannelId];
            for i=1:twix.numberOfFakeAcquisitions
                timeStampsToWrite(i) = initialTimeStamp+(i-1)*timeToAdd;
                pmuTimeStampsToWrite(i) = initialPmuTimeStamp+(i-1)*pmuTimeToAdd;
                sequenceTimeStampsToWrite(i) = initialSequenceTimeStamp+(i-1)*sequenceTimeToAdd;
            end
            timeStampsToWrite(end) = initialTimeStamp+i*timeToAdd;
            pmuTimeStampsToWrite(end) = initialPmuTimeStamp+i*pmuTimeToAdd;
            sequenceTimeStampsToWrite(end) = initialSequenceTimeStamp+i*sequenceTimeToAdd;
            acqIncrementer = 0;
            for i=1:twix.numberOfSets
                for j=1:twix.numberOfLines
                    for k=1:twix.numberOfPartitions
                        for l=1:twix.numberOfAverages
                            acqIncrementer = acqIncrementer+1;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulFlagsAndDMALength =...
                                twix.mdhScanHeader(1).ulFlagsAndDMALength;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lMeasUID =...
                                measUID;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulScanCounter =...
                                acqIncrementer;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulTimeStamp =...
                                timeStampsToWrite(acqIncrementer);
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulPMUTimeStamp =...
                                pmuTimeStampsToWrite(acqIncrementer);
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushSystemType =...
                                twix.mdhScanHeader(1).ushSystemType;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulPTABPosDelay =...
                                twix.mdhScanHeader(1).ulPTABPosDelay;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lPTABPosX =...
                                twix.mdhScanHeader(1).lPTABPosX;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lPTABPosY =...
                                twix.mdhScanHeader(1).lPTABPosY;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lPTABPosZ =...
                                twix.mdhScanHeader(1).lPTABPosZ;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulReserved1 =...
                                twix.mdhScanHeader(1).ulReserved1;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.aulEvalInfoMask =...
                                twix.mdhScanHeader(1).aulEvalInfoMask;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushSamplesInScan =...
                                twix.mdhScanHeader(1).ushSamplesInScan;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushUsedChannels =...
                                twix.mdhScanHeader(1).ushUsedChannels;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.sLC =...
                                [(j-1) (l-1) 0 (k-1) 0 0 0 (i-1) 0 0 0 0 0 0]';
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.sCutOff =...
                                twix.mdhScanHeader(1).sCutOff;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushKSpaceCentreColumn =...
                                twix.mdhScanHeader(1).ushKSpaceCentreColumn;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushCoilSelect =...
                                twix.mdhScanHeader(1).ushCoilSelect;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.fReadOutOffcentre =...
                                twix.mdhScanHeader(1).fReadOutOffcentre;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulTimeSinceLastRF =...
                                twix.mdhScanHeader(1).ulTimeSinceLastRF;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushKSpaceCentreLineNo =...
                                twix.mdhScanHeader(1).ushKSpaceCentreLineNo;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushKSpaceCentrePartitionNo =...
                                twix.mdhScanHeader(1).ushKSpaceCentrePartitionNo;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.sSD =...
                                twix.mdhScanHeader(1).sSD;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.aushIceProgramPara =...
                                twix.mdhScanHeader(1).aushIceProgramPara;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.aushReservedPara =...
                                twix.mdhScanHeader(1).aushReservedPara;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushApplicationCounter =...
                                twix.mdhScanHeader(1).ushApplicationCounter;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushApplicationMask =...
                                twix.mdhScanHeader(1).ushApplicationMask;
                            twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulCRC =...
                                twix.mdhScanHeader(1).ulCRC;
                            for m=1:twix.numberOfChannels
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulTypeAndChannelLength =...
                                    twix.mdhScanHeader(1).channelHeader(1).ulTypeAndChannelLength;
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).lMeasUID =...
                                    twix.mdhScanHeader(1).channelHeader(1).lMeasUID;
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulScanCounter =...
                                    acqIncrementer;
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulReserved1 =...
                                    twix.mdhScanHeader(1).channelHeader(1).ulReserved1;
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulSequenceTime =...
                                    sequenceTimeStampsToWrite(acqIncrementer);
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulUnused2 =...
                                    twix.mdhScanHeader(1).channelHeader(1).ulUnused2;
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulChannelId =...
                                    channelIds(m);
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulUnused3 =...
                                    twix.mdhScanHeader(1).channelHeader(1).ulUnused3;
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulCRC =...
                                    twix.mdhScanHeader(1).channelHeader(1).ulCRC;
                                readoutLineToWrite = zeros(twix.numberOfNumberSlots*twix.numberOfReadouts,1);
                                readoutLineToWrite(1:2:end) = real(twix.reconstructedKSpaceDataToWrite(i,j,k,l,m,:));
                                readoutLineToWrite(2:2:end) = imag(twix.reconstructedKSpaceDataToWrite(i,j,k,l,m,:));
                                twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).readoutLine =...
                                    readoutLineToWrite;
                            end
                        end
                    end
                end
            end
%             twix.fakeAcqEndMeasurement.ulFlagsAndDMALength = twix.mdhScanHeader(1).ulFlagsAndDMALength;
%             twix.fakeAcqEndMeasurement.fakeAcqEndMeasurementlMeasUID = measUID;
            twix.fakeAcqEndMeasurement.ulScanCounter = twix.numberOfFakeAcquisitions+1;
            twix.fakeAcqEndMeasurement.ulTimeStamp = timeStampsToWrite(twix.numberOfFakeAcquisitions+1);
            twix.fakeAcqEndMeasurement.ulPMUTimeStamp = pmuTimeStampsToWrite(twix.numberOfFakeAcquisitions+1);
%             twix.fakeAcqEndMeasurement.ushSystemType = twix.mdhScanHeader(1).ushSystemType;
%             twix.fakeAcqEndMeasurement.ulPTABPosDelay = twix.mdhScanHeader(1).ulPTABPosDelay;
%             twix.fakeAcqEndMeasurement.lPTABPosX = twix.mdhScanHeader(1).lPTABPosX;
%             twix.fakeAcqEndMeasurement.lPTABPosY = twix.mdhScanHeader(1).lPTABPosY;
%             twix.fakeAcqEndMeasurement.lPTABPosZ = twix.mdhScanHeader(1).lPTABPosZ;
%             twix.fakeAcqEndMeasurement.ulReserved1 = twix.mdhScanHeader(1).ulReserved1;
%             twix.fakeAcqEndMeasurement.aulEvalInfoMask = twix.mdhScanHeader(1).aulEvalInfoMask;
%             twix.fakeAcqEndMeasurement.ushSamplesInScan = twix.mdhScanHeader(1).ushSamplesInScan;
%             twix.fakeAcqEndMeasurement.ushUsedChannels = twix.mdhScanHeader(1).ushUsedChannels;
%             twix.fakeAcqEndMeasurement.sLC = [(j-1) (l-1) 0 (k-1) 0 0 0 (i-1) 0 0 0 0 0 0]';
%             twix.fakeAcqEndMeasurement.sCutOff = twix.mdhScanHeader(1).sCutOff;
%             twix.fakeAcqEndMeasurement.ushKSpaceCentreColumn = twix.mdhScanHeader(1).ushKSpaceCentreColumn;
%             twix.fakeAcqEndMeasurement.ushCoilSelect = twix.mdhScanHeader(1).ushCoilSelect;
%             twix.fakeAcqEndMeasurement.fReadOutOffcentre = twix.mdhScanHeader(1).fReadOutOffcentre;
%             twix.fakeAcqEndMeasurement.ulTimeSinceLastRF = twix.mdhScanHeader(1).ulTimeSinceLastRF;
%             twix.fakeAcqEndMeasurement.ushKSpaceCentreLineNo = twix.mdhScanHeader(1).ushKSpaceCentreLineNo;
%             twix.fakeAcqEndMeasurement.ushKSpaceCentrePartitionNo = twix.mdhScanHeader(1).ushKSpaceCentrePartitionNo;
%             twix.fakeAcqEndMeasurement.sSD = twix.mdhScanHeader(1).sSD;
%             twix.fakeAcqEndMeasurement.aushIceProgramPara = twix.mdhScanHeader(1).aushIceProgramPara;
%             twix.fakeAcqEndMeasurement.aushReservedPara = twix.mdhScanHeader(1).aushReservedPara;
%             twix.fakeAcqEndMeasurement.ushApplicationCounter = twix.mdhScanHeader(1).ushApplicationCounter;
%             twix.fakeAcqEndMeasurement.ushApplicationMask = twix.mdhScanHeader(1).ushApplicationMask;
%             twix.fakeAcqEndMeasurement.ulCRC = twix.mdhScanHeader(1).ulCRC;
        end
        
        function writeReconstructedKSpaceData(twix, path)
            fid = fopen(path,'r+');
            fseek(fid, twix.byteOffsetToMeasurement(twix.numberOfScans), 'bof');
            for i=1:twix.numberOfSets
                for j=1:twix.numberOfLines
                    for k=1:twix.numberOfPartitions
                        for l=1:twix.numberOfAverages
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulFlagsAndDMALength), 'uint32');
                            fwrite(fid, int32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lMeasUID), 'int32');
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulScanCounter), 'uint32');
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulTimeStamp), 'uint32');
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulPMUTimeStamp), 'uint32');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushSystemType), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulPTABPosDelay), 'uint16');
                            fwrite(fid, int32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lPTABPosX), 'int32');
                            fwrite(fid, int32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lPTABPosY), 'int32');
                            fwrite(fid, int32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.lPTABPosZ), 'int32');
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulReserved1), 'uint32');
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.aulEvalInfoMask), 'uint32');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushSamplesInScan), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushUsedChannels), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.sLC), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.sCutOff), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushKSpaceCentreColumn), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushCoilSelect), 'uint16');
                            fwrite(fid, double(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.fReadOutOffcentre), 'float');
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulTimeSinceLastRF), 'uint32');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushKSpaceCentreLineNo), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushKSpaceCentrePartitionNo), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.sSD), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.aushIceProgramPara), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.aushReservedPara), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushApplicationCounter), 'uint16');
                            fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ushApplicationMask), 'uint16');
                            fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).scanHeader.ulCRC), 'uint32');
                            for m=1:twix.numberOfChannels
                                fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulTypeAndChannelLength), 'uint32');
                                fwrite(fid, int32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).lMeasUID), 'int32');
                                fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulScanCounter), 'uint32');
                                fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulReserved1), 'uint32');
                                fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulSequenceTime), 'uint32');
                                fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulUnused2), 'uint32');
                                fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulChannelId), 'uint16');
                                fwrite(fid, uint16(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulUnused3), 'uint16');
                                fwrite(fid, uint32(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).ulCRC), 'uint32');
                                fwrite(fid, double(twix.fakeMeasurementLinesToWrite.set(i).line(j).partition(k).average(l).channelHeader(m).readoutLine), 'float');
                            end
                        end
                    end
                end
            end
            fwrite(fid, uint8(twix.acqEndMeasurement), 'uint8');
            fseek(fid,twix.byteOffsetToMeasurement(twix.numberOfScans)+...
                      twix.numberOfFakeAcquisitions*twix.byteLengthOfSingleMdhLine,'bof');
            fwrite(fid, uint32(twix.fakeAcqEndMeasurement.ulFlagsAndDMALength), 'uint32');
            fwrite(fid, int32(twix.fakeAcqEndMeasurement.lMeasUID), 'int32');
            fwrite(fid, uint32(twix.fakeAcqEndMeasurement.ulScanCounter), 'uint32');
%             fwrite(fid, uint32(twix.fakeAcqEndMeasurement.ulTimeStamp), 'uint32');
%             fwrite(fid, uint32(twix.fakeAcqEndMeasurement.ulPMUTimeStamp), 'uint32');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushSystemType), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ulPTABPosDelay), 'uint16');
%             fwrite(fid, int32(twix.fakeAcqEndMeasurement.lPTABPosX), 'int32');
%             fwrite(fid, int32(twix.fakeAcqEndMeasurement.lPTABPosY), 'int32');
%             fwrite(fid, int32(twix.fakeAcqEndMeasurement.lPTABPosZ), 'int32');
%             fwrite(fid, uint32(twix.fakeAcqEndMeasurement.ulReserved1), 'uint32');
%             fwrite(fid, uint32(twix.fakeAcqEndMeasurement.aulEvalInfoMask), 'uint32');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushSamplesInScan), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushUsedChannels), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.sLC), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.sCutOff), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushKSpaceCentreColumn), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushCoilSelect), 'uint16');
%             fwrite(fid, double(twix.fakeAcqEndMeasurement.fReadOutOffcentre), 'float');
%             fwrite(fid, uint32(twix.fakeAcqEndMeasurement.ulTimeSinceLastRF), 'uint32');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushKSpaceCentreLineNo), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushKSpaceCentrePartitionNo), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.sSD), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.aushIceProgramPara), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.aushReservedPara), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushApplicationCounter), 'uint16');
%             fwrite(fid, uint16(twix.fakeAcqEndMeasurement.ushApplicationMask), 'uint16');
%             fwrite(fid, uint32(twix.fakeAcqEndMeasurement.ulCRC), 'uint32');

%             fseek(fid,0,'bof');
%             fread(fid,1,'uint32');
%             fread(fid,1,'uint32');
%             for i=1:twix.numberOfRaidFiles
%                 fread(fid,1,'uint32');
%                 fread(fid,1,'uint32');
%                 fread(fid,1,'uint64');
%                 if i==twix.numberOfScans
%                     fwrite(fid, uint64(twix.numberOfFakeAcquisitions*twix.byteLengthOfSingleMdhLine), 'uint64');
%                 else
%                     fread(fid,1,'uint64');
%                 end
%                 fread(fid,64,'uint8=>char');
%                 fread(fid,64,'uint8=>char');
%             end

            fclose(fid);
        end
        
        function overwriteKSpaceDataInFile(twix, path)
            fid = fopen(path,'r+');
            fseek(fid, twix.byteOffsetToMeasurement(twix.numberOfScans), 'bof');
            for i=1:twix.numberOfSets
                for j=1:twix.numberOfLines
                    for k=1:twix.numberOfPartitions
                        for l=1:twix.numberOfAverages
                            fread(fid, twix.mdhScanHeaderByteLength);
                            for m=1:twix.numberOfChannels
                                readoutLineToWrite = zeros(twix.numberOfNumberSlots*twix.numberOfReadouts,1);
                                readoutLineToWrite(1:2:end) = real(twix.reconstructedKSpaceDataToWrite(i,j,k,l,m,:));
                                readoutLineToWrite(2:2:end) = imag(twix.reconstructedKSpaceDataToWrite(i,j,k,l,m,:));
                                fread(fid, twix.channelHeaderByteLength);
                                fwrite(fid, double(readoutLineToWrite), 'float');
                            end
                        end
                    end
                end
            end

            fclose(fid);
        end
        
        function extractChannelImage(twix)
            twix.channelImageDataEcho1 = Fourier.ifft3d(twix.kspaceData(:,:,:,:,1));
            twix.channelImageDataEcho2 = Fourier.ifft3d(twix.kspaceData(:,:,:,:,2));
        end

        function getSosImage(twix)
            twix.channelImageDataEcho1 = Fourier.ifft3d(twix.kspaceData(:,:,:,:,1));
            twix.channelImageDataEcho2 = Fourier.ifft3d(twix.kspaceData(:,:,:,:,2));
            twix.sosImageEcho1 = sqrt(sum(abs(twix.channelImageDataEcho1).^2, 4));
            twix.sosImageEcho2 = sqrt(sum(abs(twix.channelImageDataEcho2).^2, 4));
        end
        
        function extractReconstructedChannelImage(twix)
            twix.reconstructedChannelImageDataEcho1 = Fourier.ifft3d(twix.reconstructedKSpaceDataToWrite(:,:,:,:,1,1));
            twix.reconstructedChannelImageDataEcho2 = Fourier.ifft3d(twix.reconstructedKSpaceDataToWrite(:,:,:,:,2,1));
        end

        function getReconstructedSosImage(twix)
            twix.reconstructedChannelImageDataEcho1 = Fourier.ifft3d(twix.reconstructedKSpaceDataToWrite(:,:,:,:,1,1));
            twix.reconstructedChannelImageDataEcho2 = Fourier.ifft3d(twix.reconstructedKSpaceDataToWrite(:,:,:,:,2,1));
            twix.reconstructedSosImageEcho1 = sqrt(sum(abs(twix.reconstructedChannelImageDataEcho1).^2, 4));
            twix.reconstructedSosImageEcho2 = sqrt(sum(abs(twix.reconstructedChannelImageDataEcho2).^2, 4));
        end
    end
end