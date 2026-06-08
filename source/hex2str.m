function [s_out] = hex2str(i_in)
if (i_in < 10)
    s_out = num2str(i_in);
elseif (i_in == 10)
    s_out = 'a';
elseif (i_in == 11)
    s_out = 'b';
elseif (i_in == 12)
    s_out = 'c';
elseif (i_in == 13)
    s_out = 'd';
elseif (i_in == 14)
    s_out = 'e';
elseif (i_in == 15)
    s_out = 'f';
else
    s_out = '_';
end
