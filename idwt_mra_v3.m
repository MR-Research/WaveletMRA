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

function res_w = idwt_mra_v3(base_a, base_d, signal, level, maxlevel, support)
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
        % Here you can apply filtering/noise removal operations to last level
        % approximations and details        
        if (n < support)
            z = zeros(rows,support);
            for i = 1 : rows
                z(i,:) = ((tmpapp(i) * base_a(1:support,1)) + (tmpdet(i) * base_d(1:support,1)))'; % detail with coif
            end

            zans = zeros(1,n);
            j2 = 0;
            for k = 1 : n
                cont = 0;
                for i = 1 : rows
                    if cont == 0
                        j2 = 1;
                    else
                        j2 = cont + 1;
                    end
                    for j = 1 : support
                        if k == j2
                            zans(k) = zans(k) + z(i,j);
                        end
                        if j2 >= n
                            j2 = 1;
                        else
                            j2 = j2 + 1;
                        end
                    end
                    cont = cont + 2;
                end
            end
        else           
            zans = coef_a' * tmpapp + coef_d' * tmpdet;
        end
        res_w = zans';  
        return;
    else      
        tmpres = coef_a * sigstretch';
        tmpres2 = coef_d * sigstretch';          
        res_t = idwt_mra_v3(base_a, base_d, tmpres', level+1, maxlevel, support);     
        % Here you can apply filtering/noise removal operations to the detail coefficients of all
        % levels, res contains the previous low-level details from the
        % transform, since this is a recursive algorithm        
        if (n < support)
            z = zeros(rows,support);
            for i = 1 : rows
                z(i,:) = ((res_t(i) * base_a(1:support,1)) + (tmpres2(i) * base_d(1:support,1)))'; % detail with coif
            end
         
            zans = zeros(1,n);
            j2 = 0;
            for k = 1 : n
                cont = 0;
                for i = 1 : rows
                    if cont == 0
                        j2 = 1;
                    else
                        j2 = cont + 1;
                    end
                    for j = 1 : support
                        if k == j2
                            zans(k) = zans(k) + z(i,j);
                        end
                        if j2 >= n
                            j2 = 1;
                        else
                            j2 = j2 + 1;
                        end
                    end
                    cont = cont + 2;
                end
            end
        else           
            if (size(res_t,1) == 1)
                res_t = res_t';
            end
            zans = coef_a' * res_t + coef_d' * tmpres2;
        end
        res_w = zans';                        
    end
end