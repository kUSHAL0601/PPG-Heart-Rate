function [] = Directory(train_file_dir, test_file_dir,store_path_dir)
time_left=1200;
Files=dir(test_file_dir);
Output_file="/output_team_21_";
Files=Files(3:end);
for k=1:length(Files)
   files_left=length(Files)+1-k;
   alloted_time=time_left/files_left;
   tic;
   FileName=Files(k).name
   fs="/";
   FileNamefs=fs+FileName;
   HeartRate(train_file_dir,test_file_dir+FileNamefs,store_path_dir+Output_file+FileName,alloted_time);
   z=toc;
   time_left=time_left-z;
   files_left=files_left-1;
end
end