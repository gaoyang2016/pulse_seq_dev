classdef Fourier
    
    % Performs discrete FFT and IFFT of 3D data. Assumes inputs have
    % the first three dimensions of [ RO , PE , SE ], where RO = readout,
    % PE = phase encoding, SE = slice encoding
    
    properties
    end
    
    methods (Static)
        
        function fourier_3D = fft3d( spatial_3D )

            dims = size( spatial_3D ) ;
            fourier_3D = ifftshift( fft( fftshift( spatial_3D  ) ,...
                [] , 1 ) ) ;
            fourier_3D = ifftshift( fft( fftshift( fourier_3D ) ,...
                [] , 2 ) ) ;
            fourier_3D = 1 / sqrt( prod( dims( 1 : 3 ) ) ) *...
                ifftshift( fft( fftshift( fourier_3D ) , [] , 3 ) ) ;

        end
        
        function spatial_3D = ifft3d( fourier_3D )

            dims = size( fourier_3D ) ;
            spatial_3D = fftshift( ifft( ifftshift( fourier_3D  ) ,...
                [] , 1 ) ) ;
            spatial_3D = fftshift( ifft( ifftshift( spatial_3D ) ,...
                [] , 2 ) ) ;
            spatial_3D = sqrt( prod( dims( 1 : 3 ) ) ) *...
                fftshift( ifft( ifftshift( spatial_3D ) , [] , 3 ) ) ;

        end
        
    end
    
end

