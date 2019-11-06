function [ wavelength ] = dispersion_relation(frequency, waterdepth)
%%
%DISPERSION_RELATION computes the wavelength belonging to a wave of a
%certain frequency at a certain waterdepth
    
    if waterdepth >= 25
        wavelength = fzero(@disprel,(9.81*2*pi)/frequency^2);
    else
        wavelength = fzero(@disprel,(2*pi)/sqrt(frequency/(9.81*waterdepth)));
    end
        
    
    function [ res ] = disprel(wavelength)
        res = 9.81*((2*pi)/wavelength)*tanh((2*pi)/wavelength*waterdepth) - frequency^2;
    end
       
end

