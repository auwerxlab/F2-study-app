# F2-study-app

Interactive Shiny application to the study **_A novel mouse cross uncovers candidate therapeutic targets for hepatic steatosis, adiposity, and dyslipidemia_**.

Explore the study's data from the F2 mouse cross: metabolic phenotypes, liver RNA-seq and proteomics, and QTL/LOD genetic mapping through interactive plots and tables.

**Live app:** https://lisp-lms.shinyapps.io/F2_study/

## Running the app locally

To run it from this repository, you must download the `F2_processed_data.tar.gz` (~180 MB) from the [latest release](https://github.com/auwerxlab/F2-study-app/releases/latest) and extract it in the repository root:

```bash
curl -L -o F2_processed_data.tar.gz \
  https://github.com/auwerxlab/F2-study-app/releases/download/v1.0/F2_processed_data.tar.gz
tar -xzf F2_processed_data.tar.gz
```

This creates the `Data/` folder the app reads at startup. Then install the pinned R dependencies and launch the app:

```r
renv::restore()   # install the exact package versions from renv.lock (first run only)
shiny::runApp()
```

**Raw / source data (optional, to change parameters).** The processed bundles above are precomputed with the study's default settings. To recompute them with different parameters, download the raw source data from Zenodo:

**DOI: [10.5281/zenodo.21336282](https://doi.org/10.5281/zenodo.21336282)**

Place the files under `Data/` (extract `qtl_scans.zip` there), then re-run the relevant `Scripts/preprocess_*.R` scripts, which write updated bundles back into `Data/`. Each raw file and how it was produced is documented in the Zenodo record's README.

- For questions related to the
  - study, please contact Giorgia Benegiamo: `Giorgia.Benegiamo@epfl.ch`
  - app, please contact Alaa Badreddine: `Alaa.Badreddine@epfl.ch`
