function [NA_loc_array] = gene_NA_loc(init_type)
switch init_type
    case 0
        NA_loc_array = [];
    case 1
        NA_loc_array = 121:1: 144;
    case 2
        NA_loc_array = 61:1:84;
    case 3
        NA_loc_array = [11:12:143,12:12:144];
    case 4
        NA_loc_array = [6:12:143,7:12:145];
    case 5
        NA_loc_array = [41:12:104,42:12:104,43:12:104,44:12:104];
    case 6
        NA_loc_array = [1:12:27,2:12:27,11:12:35,12:12:36,109:12:144,119:12:144,110:12:144,120:12:144];
    case 7
        NA_loc_array = 133:1:144;
    case 8
        NA_loc_array = 61:1:72;
    case 9
        NA_loc_array = 12:12:144;
    case 10
        NA_loc_array = 6:12:144;
    case 11
        NA_loc_array = [42:12:104,43:12:104];
    case 12
        NA_loc_array = [1, 2, 11, 12, 13, 24, 121, 132, 133, 134, 143, 144];
end