classdef Bart
    
    % Performs BART reconstruction.
    
    properties
    end
    
    methods (Static)
        
        function sensitivity_maps = maps( undersampled_kspace )
            
            sensitivity_maps = bart( 'ecalib -k 6 -m 2 -S' ,...
                undersampled_kspace ) ;
            
        end
        
        function recon = pics_recon( undersampled_kspace ,...
                sensitivity_maps , lambda )
        
            command = ['pics -R W:7:0:' num2str(lambda) ' -S -d 5 -i 50'];
            
            recon = squeeze( bart( command ,...
                undersampled_kspace , sensitivity_maps ) ) ;
            
        end
        
        function recon = pi_recon( undersampled_kspace ,...
                sensitivity_maps )

              recon = squeeze( bart(...
                'pics -S -d 5 -i 50' ,...
                undersampled_kspace , sensitivity_maps ) ) ;

        end
        
        function recon = cs_recon( undersampled_kspace )
            
            for chan = 1 : size( undersampled_kspace , 4 )
                
                command = ['pics -R W:7:0:' num2str(lambda) ' -S -d 5 -i 50'];
                
                recon( : , : , : , chan ) = squeeze( bart( command ,...
                    undersampled_kspace( : , : , : , chan ) , ones(...
                    size( undersampled_kspace( : , : , : , 1 ) ) ) ) ) ;
                
            end
            
        end
        
        function kspace_per_channel = channel_data( recon ,...
                sensitivity_maps )
            
            kspace_per_channel = zeros( size( recon ) ) ;
            
            for channel = 1 : size( sensitivity_maps , 4 )
    
                kspace_per_channel( : , : , : , channel ) =...
                    recon( : , : , : , 1 )...
                    .* sensitivity_maps( : , : , : , channel ) ;
     
            end
            
        end
        
    end
    
end

