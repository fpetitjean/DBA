%*******************************************************************************
 % Copyright (C) 2013 Francois PETITJEAN, Ioannis PAPARRIZOS
 % This program is free software: you can redistribute it and/or modify
 % it under the terms of the GNU General Public License as published by
 % the Free Software Foundation, version 3 of the License.
 % 
 % This program is distributed in the hope that it will be useful,
 % but WITHOUT ANY WARRANTY; without even the implied warranty of
 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 % GNU General Public License for more details.well
 % 
 % You should have received a copy of the GNU General Public License
 % along with this program.  If not, see <http://www.gnu.org/licenses/>.
 %*****************************************************************************/ 

% Optimization of the DBA


function dba = optimdba (sequences)

%--------------------------------------------------------------------------
% Medoid Index + SOS
%--------------------------------------------------------------------------
    index = -1;
    lowestInertia = Inf;
    storage = zeros(length(sequences));
    for i=1:length(sequences)
        %-------------------Sum Of Squares --------------------------------
        sos = 0.0;
        for j=1:length(sequences)
            if i<j
                dist = dtw(sequences{i},sequences{j});
                storage(i,j) = dist;
                sos = sos + dist * dist;
            elseif i>j
                dist = storage(j,i);
                sos = sos + dist * dist;

            end
        end
        %-------------- End of sum of squares -----------------------------
        tmpInertia = sos;
        if (tmpInertia < lowestInertia)
            index = i;
            lowestInertia = tmpInertia;
        end
    end

%--------------------------------------------------------------------------
% Final output
%--------------------------------------------------------------------------
    dba = repmat(sequences{index},1);
    for i=1:5
        dba=DBA_one_iteration(dba,sequences);
    end
end

%--------------------------------------------------------------------------
% DBA_one_iteration
%--------------------------------------------------------------------------

function average = DBA_one_iteration(averageS,sequences)

    tupleAssociation = cell (1, size(averageS,2));
    for t=1:size(averageS,2)
        tupleAssociation{t}=[];
    end
    % check this sizes - in this case 4*3
    costMatrix = zeros(1000,1000);
    pathMatrix = zeros(1000,1000);

    for k=1:length(sequences)
        sequence = sequences{k};
        costMatrix(1,1) = distanceTo(averageS(1),sequence(1));
        pathMatrix(1,1) = -1;
        for i=2:size(averageS,2)
            costMatrix(i,1) = costMatrix(i-1,1) + distanceTo(averageS(i),sequence(1));
            pathMatrix(i,1) = 2;
        end

        for j=2:size(sequence,2)
            costMatrix(1,j) = costMatrix(1,j-1) + distanceTo(sequence(j),averageS(1));
            pathMatrix(1,j) = 1;
        end

        for i=2:size(averageS,2)
            for j=2:size(sequence,2)
                indiceRes = ArgMin3(costMatrix(i-1,j-1),costMatrix(i,j-1),costMatrix(i-1,j));
                pathMatrix(i,j)=indiceRes;

                if indiceRes==0
                    res = costMatrix(i-1,j-1);
                elseif indiceRes==1
                    res = costMatrix(i,j-1);
                elseif indiceRes==2
                    res = costMatrix(i-1,j);
                end

                costMatrix(i,j) = res + distanceTo(averageS(i),sequence(j));

            end
        end

        i=size(averageS,2);
        j=size(sequence,2);

        while(true)
            tupleAssociation{i}(end+1) = sequence(j);
            if pathMatrix(i,j)==0
                i=i-1;
                j=j-1;
            elseif pathMatrix(i,j)==1
                j=j-1;
            elseif pathMatrix(i,j)==2
                i=i-1;
            else
                break;
            end
        end

    end

    for t=1:size(averageS,2)
        averageS(t) = mean(tupleAssociation{t});
    end

    average = averageS;

end

function value = ArgMin3(a,b,c)

    if (a<b)
        if (a<c)
            value=0;
            return;
        else
            value=2;
            return;
        end
    else
        if (b<c)
            value=1;
            return;
        else
            value=2;
            return;
        end
    end

end


function dist = distanceTo(a,b) %if a = b dist = 0
    dist=(a-b)*(a-b);
end

# I didn't test this example
function ex = test()
    sequences = {};
    sequences{100}=[];
    for i=1:100
        length = randi(100);
        sequences{i}=rand(1,length);
    end
    mean=optimdba(sequences);
end
 

	    
