function [end_word] = output_res(time,res,tau,t_list,p_list)
%for running on cluster, output txt file 
output1 = horzcat(time,res);
output1 = horzcat(output1,tau);
output2 = horzcat(t_list,p_list);

disp(size(output2))

%save the first one
fileID_1 = fopen('res.txt','w');
[rows1, columns1] = size(output1);
for row = 1 : rows1
	for col = 1 : columns1
		fprintf(fileID_1, '%f ', output1(row, col));
	end
	fprintf(fileID_1, '\n');
end
fclose(fileID_1);

%save the second one
fileID_2 = fopen('pressure.txt','w');
[rows2, columns2] = size(output2);
for row = 1 : rows2
	for col = 1 : columns2
		fprintf(fileID_2, '%f ', output2(row, col));
	end
	fprintf(fileID_2, '\n');
end
fclose(fileID_2);


end_word = 'All files saved! Done!';

end