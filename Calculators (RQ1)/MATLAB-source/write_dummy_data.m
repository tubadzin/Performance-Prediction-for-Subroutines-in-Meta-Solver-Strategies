feats_fid = fopen('data/dummy-feats.csv','w');
fprintf(feats_fid, 'name, feat1, feat2, feat3\n')
for i=1:length(feats), fprintf(feats_fid,strcat(['inst', num2str(i), ', ', num2str(feats(i,1)), ', ', num2str(feats(i,2)), ', ', num2str(feats(i,3)), '\n'] )), end
fclose(feats_fid)

runtime_fid = fopen('data/dummy-runtime.csv','w');
for i=1:length(feats), fprintf(runtime_fid,strcat(['inst', num2str(i), ', ', num2str(feats(i,1) + 0.1*feats(i,2)), '\n'] )), end
fclose(runtime_fid)
