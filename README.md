# EucsPhylogenomics

The list of reference assemblies is available in `eucs_refseq.txt`

Folder `roadies_default` contains the output of ROADIES with the following command:
```
python run_roadies.py --cores 16 --config config.yaml
```

Folder `orthofinder` contains single-copy orthologoues and species tree from OrthoFinder with the following commands:
```
Rscript EucsPhylogenomics/extract_genes.R
OrthoFinder/orthofinder -f /data/jeremias/eucs/genes/proteome/
```

---
*Last update: 14 July 2024 by Jeremias Ivan*