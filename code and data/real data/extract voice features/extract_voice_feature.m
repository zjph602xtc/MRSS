%% set-up VoiceAnalysisToolbox
% run 'setup_voice_analysis_toolbox.m' in the 'VoiceAnalysisToolbox' folder

%% extract voice features
res = zeros(64973, 339);
name = ls('..\mpower_data\voicewave');
name = name(3:end,:);
subject=cell(64973,1);
for j=1:500:64500
    j
    parfor i=1:500
        [features, feature_names] = voice_analysis(['..\mpower_data\voicewave\' name(j+i,:)]);
        res(j+i,:) = features;
        subject{j+i} = name(j+i,:);
    end
    writematrix(res,'..\voice.csv\mpower_data\extracted features\voice_tmp.csv');
    writecell(subject,'..\voice.csv\mpower_data\extracted features\subid_tmp.csv');
end

%% go to 'download mPower Data.R' to continue