% This function courtesy of:
%   Kanghyun Chu (2019). IMAGECF Complex Field Visualization (amplitude, phase)
%   (https://www.mathworks.com/matlabcentral/fileexchange/69930-imagecf-complex-field-visualization-amplitude-phase),
%   MATLAB Central File Exchange.

% ORIGINAL LICENSE
% Copyright (c) 2019, Kanghyun Chu
% Copyright (c) 2018, Kanghyun
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


% rgb = RTh2rgb(Zn,th,bw);
function [rgb] = RTh2rgb(R_norm,Th,zc)
    sz = size(R_norm);

    R = R_norm(:);
    c = cos(Th(:));
    s = sin(Th(:));

    switch zc
        case 'w' % white
            r = 1 - (1/2 + sqrt(6)/4 * (- 2*c/sqrt(6)    )).* R ;
            g = 1 - (1/2 + sqrt(6)/4 * (+ c/sqrt(6) - s/sqrt(2) )).* R ;
            b = 1 - (1/2 + sqrt(6)/4 * (+ c/sqrt(6) + s/sqrt(2) )).* R ;

        otherwise
            r = (1/2 + sqrt(6)/4 * ( 2*c/sqrt(6)    )).*R ;
            g = (1/2 + sqrt(6)/4 * (- c/sqrt(6) + s/sqrt(2) )).*R ;
            b = (1/2 + sqrt(6)/4 * (- c/sqrt(6) - s/sqrt(2) )).*R ;
    end

    rgb = reshape( [r;g;b], [sz,3] );
end