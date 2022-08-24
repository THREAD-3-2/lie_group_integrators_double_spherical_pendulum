function B = getBlock(Mat, i, j)
%
% :param Mat: a matrix
% :param i: index
% :param j: index
%
% :returns: 3x3 block in position i-j of a matrix Mat
    
    B = Mat(3 * i - 2 : 3 * i, 3 * j - 2 : 3 * j);

end
