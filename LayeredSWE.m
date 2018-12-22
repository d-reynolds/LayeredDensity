
% RELEASE NOTES
%   Written by Dylan Reynolds (reyno18@uw.edu), Feb 2018)
%
% SYNTAX
% bulk_density = LayeredSWE(Depth,Ta,class,timeSeries,assimilationDays)%
%
% INPUTS
%   
% Depth = time series of snowpack depth in m
% Ta = air temperature in C
% class = Sturm snowpack classificaiton (see Sturm et al. 2010)
% timeSeries = matlab datenums of observations
% assimilationDays = rough estimate for number of days until new layer of
%       snow assumes bulk density of snowpack. Default value is 14 days

function bulk_density = LayeredSWE(Depth,Ta,class,timeSeries,assimilationDays)


% This function seeks to reduce spikes in SWE caused by using empirical
% density methods. These density methods often increase bulk snow density with
% snow depth. Immedietly after a new snowfall, however, this assumption may
% be voided due to thed deposition of a relatively low density layer on top
% of the existing snow. This function treats newly depositied as a separate
% layer that increases in density expotentially to assume the bulk snowpack
% density

%Arbitrary depth limit to determine if a new layer was formed in the depth
%record
DEPTH_LIM = 0.03;
T_LIM = 3;

%Timestep of depth record, in hours
TIME_STEP = 1;

if nargin < 5
    assimilationDays = 14; 
end

bulk_density = zeros(length(Depth),length(assimilationDays));

for k = 1:length(assimilationDays)
    
    %Initialize arrays
    layers = [];
    layerNewDen = [];
    layerAges = [];
    layerDepths = [];
    layerDepthInitial = [];
        
    last_Depth = 0;
    last_Density = 0;

    %If time series is in the time_builder format, extract just the datenum
    %portion
    if size(timeSeries,2) > 1
        timeSeries = timeSeries(:,7);
    end
    
    for i = 1:length(bulk_density)
        %Bulk snow density based on snow-class method (Sturm et al. 2010)
        sturm_bulk = CalcSturmDensity(Depth(i).*100,timeSeries(i),class);
        
        %% Do layer evolution
        layer_n = length(layers);
        if layer_n > 0
            for j = 1:layer_n
                %If layer age > assimilationDays then layer can be assimilated into snowpack
                if layerAges(j) >= assimilationDays(k)
                    layers = layers(1:j-1);
                    layerNewDen = layerNewDen(1:j-1);
                    layerAges = layerAges(1:j-1);
                    layerDepths = layerDepths(1:j-1);
                    layerDepthInitial = layerDepthInitial(1:j-1);
                else
                    %Age layer accordingly and calculate "assimilation factor"
                    %which determines how much the new layer looks like old
                    %snowpack
                    layerAges(j) = layerAges(j) + ((1/24)*TIME_STEP);
                    assimilation_factor = log(1 + ((exp(1)-1)*layerAges(j)/assimilationDays(k)));
                    
                    layers(j) = layerNewDen(j) + ...
                        (sturm_bulk-layerNewDen(j))*assimilation_factor;
                    
                    %final depth is the depth of the layer when it has reached
                    %density of old snowpack. Assuming SWE is conserved
                    layerDepthFinal = layerDepthInitial(j)*layerNewDen(j)/sturm_bulk;
                    layerDepths(j) = layerDepthInitial(j) - ...
                        (layerDepthInitial(j)-layerDepthFinal)*assimilation_factor;
                end
            end
            
        end
        
        %% New Snow Density modeling (From Hedstrom and Pomeroy 1998)
        %See if we had accumulation and add a new layer if so
        if (Depth(i) - last_Depth) > DEPTH_LIM && Ta(i) < T_LIM 
            newSnowDen = 0.0679 + 0.0513*exp(Ta(i)/2.6);
            layerNewDen = [newSnowDen; layerNewDen];
            layers = [layerNewDen(1); layers];
            layerAges = [0; layerAges];
            layerDepthInitial = [Depth(i) - last_Depth; layerDepthInitial];
            layerDepths = [layerDepthInitial(1); layerDepths];
        end
        
        %%
        %calculate bulk density as average of all layer density weighted by depth
        if Depth(i) > 0
            bulk_density(i,k) = (sum(layerDepths.*layers) + ...
                max((Depth(i)-sum(layerDepths)),0)*sturm_bulk)/Depth(i);
            if ~isnan(Depth(i))
                last_Depth = Depth(i);
                last_Density = bulk_density(i,k);
            end
        else
            bulk_density(i,k) = last_Density;
        end
    end
end

end

