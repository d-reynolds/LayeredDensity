% RELEASE NOTES
%   Written by Dylan Reynolds (reyno18@uw.edu), Feb 2018)
%
% SYNTAX
% SturmDOY = GetSturmDOY(DateNumVector)

% INPUTS
% datenums = matlab datenums of observations

function SturmDOY = GetSturmDOY(datenums)

% Calculates day of year (DOY) according to the scheme in Sturm et
% al. 2010
LEAP = [2000, 2004, 2008, 2012, 2016, 2020];
LONG_MONTHS = [1, 3, 5, 7, 8, 10, 12];
SHORT_MONTHS = [4, 6, 9, 11];
SturmDOY = zeros(length(datenums),1);

for i = 1:length(SturmDOY)
    [Y, M, D] = datevec(datenums(i));
    DOY = 0;
    %generate actual day of year
    for j = 1:M
        if j == M
            DOY = DOY + D;
        elseif find(j == LONG_MONTHS)
            DOY = DOY + 31;
        elseif find(j == SHORT_MONTHS)
            DOY = DOY + 30;
        elseif j == 2
            if find(Y == LEAP)
                DOY = DOY + 29;
            else
                DOY = DOY + 28;
            end
        end
    end
    
    
    %post-analysis to make sure we aren't in 'negative' months
    if M >= 10
        if find(Y == LEAP)
            DOY = DOY - 367;
        else
            DOY = DOY - 366;
        end
    end
    SturmDOY(i,1) = DOY;
end

end

