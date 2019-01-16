function termNum = getTermNumber(coords, degs, coordTrie, termInd)

% Follow the trie to the end of this term
curTrieInd = 1;
len = length(coords);
for j=1:len
    c = coords(j);
    % Get the index for this coordinate
    nextTrieInds = coordTrie{curTrieInd, 3};
    curTrieInd = nextTrieInds(c);
end

% Get the degrees for this set of coordinates and compare to given degrees
curDegMat = coordTrie{curTrieInd, 1};
loc = (curDegMat(:, 1) == degs(1));
for j=2:len
    loc = loc & (curDegMat(:, j) == degs(j));
end

curIndMat = coordTrie{curTrieInd, 2};
curInd = curIndMat(loc);
termNum = termInd(curInd);

end