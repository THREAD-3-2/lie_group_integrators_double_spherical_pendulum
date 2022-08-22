function res = actionSE3(B, input)
% action on the SE(3) group
%
% :param B: element of the Lie group SE(3)
% :param input: element of the Lie algebra se(3)
%
% :returns: element of the Lie algebra se(3)

res = [B(1:3, 1:3) * input(1:3); B(1:3, 1:3) * input(4:6) + cross(B(:, end), (B(1:3, 1:3) * input(1:3)))];

end
