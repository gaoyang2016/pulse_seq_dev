function [mdh_blob, filePos, isEOF] = loop_mdh_read( fid )
% Goal of this function is to gather all mdhs in the dat file and store them
% in binary form, first. This enables us to evaluate and parse the stuff in
% a MATLAB-friendly (vectorized) way. We also yield a clear separation between
% a lengthy loop and other expressions that are evaluated very few times.
%
% The main challenge is that we never know a priori, where the next mdh is
% and how many there are. So we have to actually evaluate some mdh fields to
% find the next one.
%
% All slow things of the parsing step are found in the while loop.
% => It is the (only) place where micro-optimizations are worthwhile.
%
% The current state is that we are close to sequential disk I/O times.
% More fancy improvements may be possible by using workers through parfeval()
% or threads using a java class (probably faster + no toolbox):
% http://undocumentedmatlab.com/blog/explicit-multi-threading-in-matlab-part1

    byteMDH = 184;
    szScanHeader    = 192; % [bytes]
    szChannelHeader =  32; % [bytes]

    cPos            = ftell(fid);
    n_acq           = 0;
    allocSize       = 4096;
    ulDMALength     = byteMDH;
    isEOF           = false;
    percentFinished = 0;
    progress_str    = '';
    prevLength      = numel( progress_str );

    mdh_blob = zeros( byteMDH, 0, 'uint8' );
    szBlob   = size( mdh_blob, 2 );
    filePos  = zeros(0, 1, class(cPos));  % avoid bug in Matlab 2013b: https://scivision.co/matlab-fseek-bug-with-uint64-offset/

    % get file size
    fseek(fid,0,'eof');
    fileSize = ftell(fid);
    fseek(fid,cPos,'bof');

    % ======================================
    %   constants and conditional variables
    % ======================================
        bit_0 = uint8(2^0);
        bit_5 = uint8(2^5);
        mdhStart = 1-byteMDH;
        
        u8_000 = zeros( 3, 1, 'uint8'); % for comparison with data_u8(1:3)

        % 20 fill bytes in VD (21:40)
        evIdx   = uint8(    21  + 20); % 1st byte of evalInfoMask
        dmaIdx  = uint8((29:32) + 20); % to correct DMA length using NCol and NCha
        dmaOff  = szScanHeader;
        dmaSkip = szChannelHeader;
    % ======================================

    t0 = tic;
    while true
        % Read mdh as binary (uint8) and evaluate as little as possible to know...
        %   ... where the next mdh is (ulDMALength / ushSamplesInScan & ushUsedChannels)
        %   ... whether it is only for sync (MDH_SYNCDATA)
        %   ... whether it is the last one (MDH_ACQEND)
        % evalMDH() contains the correct and readable code for all mdh entries.
        try
            % read everything and cut out the mdh
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
        catch exc
            warning( [mfilename() ':UnxpctdEOF'],  ...
                      [ '\nAn unexpected read error occurred at this byte offset: %d (%g GiB)\n'...
                        'Will stop reading now.\n'                                             ...
                        '=== MATLABs error message ================\n'                         ...
                        exc.message                                                            ...
                        '\n=== end of error =========================\n'                       ...
                       ], cPos, cPos/1024^3 )
            isEOF = true;
            break
        end

        bitMask = data_u8(evIdx);   % the initial 8 bit from evalInfoMask are enough

        if   isequal( data_u8(1:3), u8_000 )    ... % probably ulDMALength == 0
          || bitand(bitMask, bit_0)                % MDH_ACQEND

            % ok, look closer if really all *4* bytes are 0:
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );

            if ulDMALength == 0 || bitand(bitMask, bit_0)
                cPos = cPos + ulDMALength;
                % jump to next full 512 bytes
                if mod(cPos,512)
                    cPos = cPos + 512 - mod(cPos,512);
                end
                break;
            end
        end
        if bitand(bitMask, bit_5)  % MDH_SYNCDATA
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            cPos = cPos + ulDMALength;
            continue
        end

        % pehses: the pack bit indicates that multiple ADC are packed into one
        % DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
        % since this code assumes one adc (x NCha) per DMA, we have to correct
        % the "DMA length"
        %     if mdh.ulPackBit
        % it seems that the packbit is not always set correctly
        NCol_NCha = double( typecast( data_u8(dmaIdx), 'uint16' ) );  % [ushSamplesInScan  ushUsedChannels]
        ulDMALength = dmaOff + (8*NCol_NCha(1) + dmaSkip) * NCol_NCha(2);

        n_acq = n_acq + 1;

        % grow arrays in batches
        if n_acq > szBlob
            mdh_blob( :, end + allocSize ) = 0;
            filePos( end + allocSize ) = 0;
            szBlob = size( mdh_blob, 2 );
        end
        mdh_blob(:,n_acq) = data_u8;
        filePos( n_acq )  = cPos;

        if (100*cPos)/fileSize > percentFinished + 1
            percentFinished = floor((100*cPos)/fileSize);
            elapsed_time    = toc(t0);
            time_left       = (fileSize/cPos-1) * elapsed_time;
            prevLength      = numel(progress_str);
            progress_str    = sprintf('    %3.0f %% read in %4.0f s; estimated time left: %4.0f s \n',...
                                      percentFinished,elapsed_time, time_left);
            fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
        end

        cPos = cPos + ulDMALength;
    end % while true

    if isEOF
        n_acq = n_acq-1;    % ignore the last attempt
    end

    filePos( n_acq+1 ) = cPos;  % save pointer to the next scan

    % discard overallocation:
    mdh_blob = mdh_blob(:,1:n_acq);
    filePos  = reshape( filePos(1:n_acq+1), 1, [] ); % row vector

    prevLength   = numel(progress_str) * (~isEOF);
    elapsed_time = toc(t0);
    progress_str = sprintf('    100 %% read in %4.0f s; estimated time left:    0 s \n', elapsed_time);
    fprintf([repmat('\b',1,prevLength) '%s'],progress_str);

end % of loop_mdh_read()