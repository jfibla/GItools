# Inspector_resources (local data, not in GitHub)

The content of this folder is **NOT included in the GitHub repository**.  
It is a local “resources” directory that each user must populate in order to run some GItools modules (LD, NonSyn, EWAS, etc.).

You can downoad a compresed file containing all **Inspector_resources** folder and files at figshare

GItools will look here by default:
`<GItools_root>/Inspector_resources/`

If you prefer to store resources elsewhere (e.g., external drive), set:
- `GITOOLS_RESOURCES=/path/to/Inspector_resources`


## Recommended structure

```
Inspector_resources/
  POP/
    EUR.txt
    AFR.txt
    AMR.txt
    EAS.txt
    SAS.txt
    MED.txt
    CSA.txt

  LD_resources/
    <BFILE_PREFIX>.bed
    <BFILE_PREFIX>.bim
    <BFILE_PREFIX>.fam
    # optional: .hh, .log, etc.

  software/
    plink19/
      plink              # mac/linux executable (or plink.exe on Windows)

  dbNSFP5_broad_hg38/
    output_snps_biallelic_dbNSFP5_broad_hg38.ALL.out.gz
    output_snps_biallelic_dbNSFP5_broad_hg38.ALL.out.gz.tbi

  EWAS_cancer/
    ewas_detail_tum_genome.rds
    slim_tc_by_chr/      # large per-chromosome files

  EWAS_disease/
    ewas_detail_dis_genome.rds
    slim_disease_by_chr/ # large per-chromosome files
```


## Notes

### POP keep files (`Inspector_resources/POP/*.txt`)
Each file must contain either:
- **1 column**: `IID`  
or
- **2 columns**: `FID IID`

Example (`EUR.txt`):
```
FID1 IID1
FID2 IID2
...
```


## Optional environment variables

You can override defaults with:

- `GITOOLS_RESOURCES=/path/to/Inspector_resources`
- `GITOOLS_POP_DIR=/path/to/Inspector_resources/POP`
- `GITOOLS_LD_BFILE=/path/to/Inspector_resources/LD_resources/<BFILE_PREFIX>`
- `GITOOLS_PLINK19=/path/to/Inspector_resources/software/plink19/plink`
