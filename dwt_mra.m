%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function res_w = dwt_mra(base_a, base_d, signal, level, maxlevel, support, res)
    n = size(signal, 2);
    translation = 2;
    rows = n/2;
    cols = n; 
    if (n < support)
        colnew = 0;
        while (colnew < support)
            colnew = colnew + (n + n + (n/2));
        end
        cols = colnew*2;
    end
    coef_a = zeros(rows,cols);
    coef_d = zeros(rows,cols);     
    for i = 1 : rows
        for j = 1 : support 
            if (((i-1)*translation)+j > cols)
                despl = (((i-1)*translation)+j)-(cols);
                coef_a(i, despl) = base_a(j);
                coef_d(i, despl) = base_d(j);                    
            else
                despl = ((i-1)*translation)+j;
                coef_a(i,despl) = base_a(j);
                coef_d(i,despl) = base_d(j);
            end
        end            
    end        
    sigstretch = signal;
    if (n < support)
        it = 1;
        sigloop = 1;
        padding = [];
        while (size(padding,2) < (cols-size(signal,2)))         
            padding(it) = signal(sigloop);
            if (sigloop == size(signal,2))
                sigloop = 1;
            else
                sigloop = sigloop + 1;
            end
            it = it + 1;            
        end
        sigstretch = [signal padding];
    end          
    if level == maxlevel     
        tmpapp = coef_a * sigstretch';
        tmpdet = coef_d * sigstretch';
        % Here you can apply operations to last level
        % approximations and details
        res(level+1,1:size(tmpapp,1)) = tmpapp;
        res(level,1:size(tmpdet,1)) = tmpdet;             
        res_w = res;
        return;
    else      
        tmpres = coef_a * sigstretch';
        tmpres_d = coef_d * sigstretch';
        res(level,1:size(tmpres_d, 1)) = tmpres_d';            
        res = dwt_mra(base_a, base_d, tmpres', level+1, maxlevel, support, res);
        % Here you can apply operations to the detail coefficients of all
        % levels, res contains the previous low-level details from the
        % transform, since this is a recursive algorithm
        res_w = res;                        
    end
end