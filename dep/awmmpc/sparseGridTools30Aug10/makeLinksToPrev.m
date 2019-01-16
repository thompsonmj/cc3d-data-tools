function [linksToPrev, coordTrie] = makeLinksToPrev(dim, termDegs, termCoords, termInd, termStart, termLen, numTerms)

% Make a trie to hold the info for each term.  
coordTrie = makeCoordTrie(dim, termDegs, termCoords, termStart, termLen, numTerms); 

linksToPrev = cell(numTerms, 1);
% Loop through for each set of coordinates
for trieInd = 2:size(coordTrie, 1)
    curDegMat = coordTrie{trieInd, 1};
    curTermIndMat = termInd(coordTrie{trieInd, 2});
    prevEntries = coordTrie{trieInd, 4};
    numCoords = size(curDegMat, 2);
    links = zeros(size(curDegMat));
    
    % Find the sets that point to previous entries
    for j=1:numCoords
        curRows = (curDegMat(:, j) == 1);
        prevDegMat = coordTrie{prevEntries(j), 1};
        prevTermIndMat = termInd(coordTrie{prevEntries(j), 2});
        
        % Match the current rows with the previous entries
        % First do rows that point to other trie entries (deg = 1)
        if (numCoords > 1)
            newDegs = curDegMat(curRows, [1:j-1, j+1:numCoords]);
            % We could use intersect as here:
            % [c, ia, ib] = intersect(prevDegMat, curDegMat(curRows, [1:j-1, j+1:numCoords]), 'rows');
            % but prevDegMat and curDegMat already have sorted rows, so the
            % following excerpt from intersect.m is faster:            
            [ia, ib] = intersectSortedRows(prevDegMat, newDegs);
            curRowInds = find(curRows);
            links(curRowInds(ib), j) = prevTermIndMat(ia);
        else
            links(curRows) = 1;
        end
        % Then others (deg > 1)
        newDegMat = curDegMat(~curRows, :);
        newDegMat(:, j) = newDegMat(:, j) - 1;
        [ia, ib] = intersectSortedRows(curDegMat, newDegMat);
        newRowInds = find(~curRows);
        links(newRowInds(ib), j) = curTermIndMat(ia);       
    end
    
    linksToPrev(curTermIndMat) = mat2cell(links, ones(size(links, 1), 1), size(links, 2));
end
end

% Find the indices for the intersecting rows of two matrices whose rows are
% already sorted.
function [ia, ib] = intersectSortedRows(A, B)
    % d = 1 if differences between rows, 0 if equal.
    rows = size(B, 1);
    d = (B(1:rows-1,:)~=B(2:rows,:));
    d = any(d,2);
    d(rows,1) = true;
    ib = find(d);

    [c, ndx] = sortrows([A; B]);
    rowsC = size(c, 1);
    d = (c(1:rowsC-1,:) == c(2:rowsC,:));
    d = find(all(d,2));

    n = size(A,1);
    ia = ndx(d);
    ib = ib(ndx(d+1)-n);
end

function coordTrie = makeCoordTrie(dim, termDegs, termCoords, termStart, termLen, numTerms) 
% Build a trie to store info about the terms
% coordTrie{j, 1} contains an nxk matrix with rows [deg1, ..., degk] for each
%    terminal branch corresponding to a term that depends only on the
%    coordinates leading to this entry, has the given degrees, and starts at polyDeg(index)
% coordTrie{j, 2} contains an nx1 vector giving the term index for each row
%    in  coordTrie{j, 1}
% coordTrie{j, 3} contains a column vector [ind1;...; indd] of indices into coordTrie, one
%    for each of coordinates 1, ..., d. These point to the entry obtained
%    from the current entry by adding one additional variable.
% coordTrie{j, 4) contains a row vector [ind1,..., indk] of indices into
%    coordTrie, one for each coordinate appearing in this entry.  These
%    point to the entry obtained from the current entry by removing one
%    variable.

numToAdd = max(nchoosek(dim, 2), 20);
coordTrie = cell(numToAdd, 3);
coordTrie{1, 1} = 0;
coordTrie{1, 2} = 1;
nextOpenTrieInd = 2;
maxTrieInd = numToAdd;
trieEntryCoords = cell(numToAdd, 1);

% Sort the terms based on the the coordinates in the term
for len=1:max(termLen)
    termsThisLen = (termLen == len);
    termsThisLen(1) = 0; % exclude the constant term
    termsThisLenInds = find(termsThisLen);
    
    % Get all sets of coordinates of this length
    coordMat = cell2mat(termCoords(termsThisLen));
    coordMat = reshape(coordMat, len, [])';
    
    % Determine unique sets of coords and work on each one
    [uniqueCoords, m, n] = unique(coordMat, 'rows');    
    for row=1:size(uniqueCoords, 1)
        
        curTrieInd = 1;
        % Get the coordinates and length for this term
        curCoords = uniqueCoords(row, :);

        % Build the trie to the end of this term
        for c = curCoords

            % Get the index for this coordinate
            nextTrieInds = coordTrie{curTrieInd, 3};
            % Make a vector of indices if needed
            if (isempty(nextTrieInds))
                nextTrieInds = zeros(dim, 1);
            end
            nextTrieInd = nextTrieInds(c);

            % Make a trie entry if needed
            if (nextTrieInd == 0)
                nextTrieInds(c) = nextOpenTrieInd;
                nextTrieInd = nextOpenTrieInd;
                coordTrie{curTrieInd, 3} = nextTrieInds;

                % Make extra space if needed
                if (nextOpenTrieInd > maxTrieInd)
                    maxTrieInd = maxTrieInd + numToAdd;
                    coordTrie = [coordTrie; cell(numToAdd, 4)];
                    trieEntryCoords = [trieEntryCoords; cell(numToAdd, 1)];
                end
                nextOpenTrieInd = nextOpenTrieInd + 1;
            end
            curTrieInd = nextTrieInd;
        end

        % Save the degrees and indices in the correct entry        
        curTerms = termsThisLenInds(n == row);
        curDegMat = cell2mat(termDegs(curTerms));
        curDegMat = reshape(curDegMat, len, [])';        
        curIndMat = termStart(curTerms);
        [sortedDegMat, sortInd] = sortrows(curDegMat);

        coordTrie{curTrieInd, 1} = sortedDegMat;
        coordTrie{curTrieInd, 2} = curIndMat(sortInd);
        coordTrie{curTrieInd, 4} = zeros(1, length(curCoords));
        trieEntryCoords{curTrieInd} = curCoords;        
    end
    
end

% Make links from each set of coordinates to each set of coordinates 
% obtained by removing one of the coordinates.
numEntries = nextOpenTrieInd-1;
coordTrie = coordTrie(1:numEntries, :);
trieEntryCoords = trieEntryCoords(1:numEntries);
numCoordsAll = cellfun(@length, trieEntryCoords);

% Loop through each trie entry
for j=2:numEntries
    % Get the coordinates for this entry
    curCoords = trieEntryCoords{j};
    n = length(curCoords);
    % If only one, then point to the constant
    if (n == 1)
        coordTrie{j, 4}(1) = 1;
        continue;
    end
    % Find entries with one fewer coordinate and get those coordinates
    prevLevel = (numCoordsAll == n-1);
    prevLevelInds = find(prevLevel);
    prevLevelCoords = cell2mat(trieEntryCoords(prevLevel));

    % Remove one coordinate at a time
    for k=1:n
        prevCoords = curCoords([1:k-1, k+1:n]);
        % Find the entry that matches this reduced set of coordinates
        loc = (prevLevelCoords(:, 1) == prevCoords(1));
        for m=2:n-1
            loc = loc & (prevLevelCoords(:, m) == prevCoords(m));
        end
        % Link to that entry
        prevInd = prevLevelInds(loc);
        coordTrie{j, 4}(k) = prevInd;
    end
end

end
