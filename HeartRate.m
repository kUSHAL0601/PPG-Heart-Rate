function [] = HeartRate(train_file, test_file,store_path,time_left)
load(test_file);
pred=FindBPM(sig,time_left).';
save(store_path,'pred');
end
