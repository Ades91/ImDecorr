% [ind,A] = getDcorrMax(d)
% ---------------------------------------
%
% Return the local maxima of the decorrelation function d
%
% Inputs:
%  d        	Decorrelation function
%
% Outputs:
%  ind        	Position of the local maxima
%  A			Amplitude of the local maxima
%
% ---------------------------------------
%
% A detailled description of the method can be found in : 
% "Descloux, A., K. S. Grussmayer, and A. Radenovic. "Parameter-free image 
% resolution estimation based on decorrelation analysis."
% Nature methods (2019): 1-7."
% 
%   Copyright ?? 2018 Adrien Descloux - adrien.descloux@epfl.ch, 
%   Ecole Polytechnique Federale de Lausanne, LBEN,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.
%
%  	This program is free software: you can redistribute it and/or modify
%  	it under the terms of the GNU General Public License as published by
% 	the Free Software Foundation, either version 3 of the License, or
%  	(at your option) any later version.
%
%  	This program is distributed in the hope that it will be useful,
%  	but WITHOUT ANY WARRANTY; without even the implied warranty of
%  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  	GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%  	along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [ind,A] = getDcorrLocalMax(d)

Nr = numel(d);
if Nr < 2
    ind = 1;
    A = d(1);
else
    % find maxima of d
    [A,ind] = max(d);
    while numel(d) > 1
         if ind == numel(d) 
            d(end) = [];
            [A,ind] = max(d);
         elseif ind == 1
             break;
         elseif (A - min(d(ind:end))) >= 0.0005
             break;
         else
             d(end) = [];
            [A,ind] = max(d);
         end     
    end
    if isempty(ind)
        ind = 1;
        A = d(1);
    else
        A = d(ind);
        ind = ind(end);
        
    end
end