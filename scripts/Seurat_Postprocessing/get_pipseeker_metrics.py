import pandas as pd

# each sample has a metrics folder 
# each metrics folder has a file for each sensitivity, barcode_stats.csv, matrix_stats.csv

## path2pipseeker results dir
#path2results = "/data/250211_pipseq_organoids/250211_pipseq_results"
## list of sample names 
#samplenames = ['DSO_D0','DSO_D3','DSO_D6','DSO_D9','HSO_D0','HSO_D3','HSO_D6','HSO_D9']

samplenames = []
with open("filelist_path_samplenicknames.txt",'r') as filelist:
    for f in filelist:

        path2results, samplename, samplenickname = f.strip().split()
        
        p = F"{path2results}/metrics"
        samplenames.append(p)

        # combine metrics within one sensitivity level into one table
        for n in [1,2,3,4,5]:
            barcode_stats = pd.read_csv(F"{p}/barcode_stats.csv",header=None)
            matrix_stats = pd.read_csv(F"{p}/matrix_stats.csv",header=None)
            sensitvity_stats = pd.read_csv(F"{p}/sensitivity_{n}/metrics_summary.csv",header=None)

            df = pd.concat([sensitvity_stats,barcode_stats, matrix_stats])
            df.loc[len(df.index)] = ['sensitivity',n]
            df.loc[len(df.index)] = ['sample_name',samplename]
            df.loc[len(df.index)] = ['sample_nickname',samplenickname]
            df.to_csv(F"{p}/sensitivity_{n}/combined_stats.csv",index=False,header=None)
        # combine data from all 5 sensitivity levels 
        df_tmps = []
        for n in [1,2,3,4,5]:
            df = pd.read_csv(F"{p}/sensitivity_{n}/combined_stats.csv",header=None)
            df.index=df[0]
            df_tmps.append(df[1])
        df_new = pd.concat(df_tmps,axis=1).T
        df_new.to_csv(F"{p}/combined_stats.csv",index=False)

    # combine data from all samplenames into one
    df_allSamples = pd.concat([pd.read_csv(F"{s}/combined_stats.csv") for s in samplenames])
    df_allSamples.to_csv(F"TABLE_METRICS_ALL_SAMPLES.csv",index=False)