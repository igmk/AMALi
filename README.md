# AMALi_processing
Content:

Amali_subroutines : folder that contains routine used by `amali_eval_Wolkenur532.m`

Input : input data for the processing

.gitignore : `.gitignore` file

KARL_background.m : calculates the atmospheric background

amali_eval_Wolkenur532.m : does the processing following the Klett algorithm

amalifaul.m : calls `amali_eval_Wolkenur532.m`

amalireadraw_lukas.m : reads the AMALi binary files

readamali_lukas.m : routine to read AMALi binary data and concatenate them into one matlab file. uses `amalireadraw_lukas.m`

## processing steps
Measurements taken by AMALi need to be preprocessed by `readamali_lukas.m`. Resulting files can be used be `amalifaul.m`.
