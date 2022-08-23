function B = getBlock(Mat,i,j)
    
    %From a matrix Mat it simply gives out the 3x3 block in position i-j

    B = Mat(3*i-2:3*i, 3*j-2:3*j);

end