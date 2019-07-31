Project 18: Data pipeline for extracting signatures of evolutionary dynamics from tumours
======
This is a REAMDE.md file decribing a piece of software that is accompanying the submission report of University of Glasgow MSc Data Science programme's MSc Team Project course, authored by Razak Nart, Mingfeng Liu, Lukas Rubikas, Yuhong Lin and Jesus Vera and supervised by Dr Ke Yuan.

Our task, described in the report and implement in this pipeline, is twofold:  

1. Reconstructing the subclonal composition of tumours in terms of clusters of simple somatic mutations (SSMs) and their corresponding variant allele frequencies (VAFs) using the following tools and libraries:
    1. PyClone, a Python library developed by Andrew Roth et al. and adapted to this data pipeline by Mingfeng Liu,  
    2. Ccube, a R library developed by Ke Yuan and adapted to this data pipeline by Yuhong Lin,  
    3. PhyloWGS, an end-to-end open-source software written with Python/C++, developed by Quaid Morris et al. and adapted to this data pipeline by Lukas Rubikas,  
2. Performing Latent Dirichlet allocation (LDA), implemented for this project by Jesus Vera, and non-negative matrix factorization (NMF), implemented by Razak Nart, using these mutation assignments for each and every SSM in each every tumour sample.

Is this therefore why we will refer to action steps described in (1) as the Stage One of the project and action steps described in (2) as Stage Two.

Our data
-----------
The test data that was used while developing this tool were generated and supplied to us by the project supervisor, Dr Ke Yuan, and was structured as following below:

`data/`:

* `VCF/` folder, containing the mutation information for every tumour sample, formatted as `SAMPLE_NAME.no_real_info.vcf`, in a widely recognized VCF format,  
* `Segments/` folder, containing copy number profiles for each tumour sample, formatted as `SAMPLE_NAME.segmets.txt`,  
* `pp_table.txt` text file, containing the purity and ploidy information for each sample.

It is therefore our requirement that in order to run our pipeline, the input data should follow the same pattern and structure.

Moreover, when working on Stage Two of the project, we needed to have data already processed by Stage One with the 500 result files that would include the proportion values. In order to do this, ccube was run (as it was the fastest one to run). The data is located as follows:

`data/Data/Phase_two/`:

* `Subclonal_Structure/` folder, containing the 500 subclonal structure files, formatted as `SAMPLE_NAME_subclonal_structure.txt`.

Our project structure
-----------
In this data pipeline we combine a number of different tools and algorithms, therefore it is necessary to be aware of how our source code is structured and how specific parts of it can be accessed, modified or updated, should there be a need for it. Our source code is therefore structured as the following:

`team-project-cs/`:

* `docker/`, containing necessary files and instructions how to set up a Docker image for our pipeline with the numerous dependencies we require
* `Pipeline/`, all of the source code of Stage One and Stage Two  
  * `data/`, the input folder for Stage One and Two, with the folder structure and file-naming pattern described in the "Our data" section in this README.md,
  * `ccube/`, containing all relevant files and folders used by/generated while using Ccube
  * `pycloneL/`, containing source code of our adaptation of PyClone library and the post-proccessing source code of its interim output,
  * `phylowgs/`, containing files associated with PhyloWGS-related part of our project, most notably:
    * `phylowgs.source.code.project18`, containing cloned GitHub repo of PhyloWGS with our modifications to make it adaptable to our project,
    * `phylowgs.results` (generated), containing the final output of PhyloWGS (but rather interim output for our pipeline) should a user would like to use PhyloWGS's own inference about the tumour subclonal composition capabilities
  * `LDA/` 
    * `Code/` source code of our implementation of the LDA algorithm, part of Stage Two
  * `NMF/`
    * `Code/` source code of our implementation of the NMF algorithm, part of Stage Two
  * `inputTmp/` (generated), temporary folder for storing interim PyClone results (can be disabled)
  * `results/` (generated), containing the results, converted to a common format, of each of the Stage One tools and possible input folders of Stage Two algorithms:
    * `ccube/`
      * `multiplicity/` (unused) sample SSMs multiplicities, as inferred by Ccube
      * `mutation_assignments/` sample SSMs mutation assignments, as inferred by Ccube
      * `subclonal_structure/` sample subclonal structure, as inferred by Ccube
      * `runtimes.tsv`, a merged file containing the Ccube running times (in secods) for each sample
    * `pyclone/`
      * `multiplicity/` (unused) sample SSMs multiplicities, as inferred by PyClone
      * `mutation_assignments/` sample SSMs mutation assignments, as inferred by PyClone
      * `subclonal_structure/` sample subclonal structure, as inferred by PyClone
      * `runtimes.tsv`, a merged file containing the PyClone running times (in secods) for each sample
    * `phylowgs/` *[1]
      * `mutation_assignments/` sample SSMs mutation assignments, as inferred by PhyloWGS
      * `subclonal_structure/` sample subclonal structure, as inferred by PhyloWGS
      * `runtimes.tsv`, a merged file containing the PhyloWGS running times (in secods) for each sample
   * `pipeline.py`, the entry point of using the pipeline, combining all of the programming logic of our tools.

*[1] `multiplicity/` folder is skipped by the common format output generator of the PhyloWGS interim output

Our output
-----------
Since all of the Stage One tools outputs their results in a widely different formats, it was necessary to agree upon a common output, which our resulting data had to be transformed to. The common format was suggested by the project supervisor Dr Ke Yuan and must follow such pattern:

1. In `subclonal_structure/` folder, for each sample in the input data, resulting file must be named as `SAMPLE_NAME_subclonal_structure.txt` and must have the following collumns:  
    | cluster  | n_ssms | proportion | ccf (optionally) |  
    | - | - | - | - |  
  
2. In `mutation_assignments` folder, for each sample in the input data, resulting file must be named as `SAMPLE_NAME_mutation_assignments.txt` and must have the following collumns:  
   | chr  | pos | cluster | proportion | ccf (optionally) *[2] |  
   | - | - | - | - | - |  
  
3. In `multiplicity` *[3] folder (optional), for each sample in the input data, resulting file must be named as `SAMPLE_NAME_multiplicity.txt` and must have the following collumns:  
   | chr  | pos | tumour_copynumber | multiplicity |  
   | - | - | - | - |  

*[2] in Ccube common format results, this collumn is named `ccfmean`, in PyClone - `average_ccf`, in PhyloWGS it is left as (and calculated as) just `ccf`.  
*[3] `multiplicity/` data is unused in our project and was skipped in the common output of PhyloWGS results.  

The output for Stage Two tools is very different from the previous Stage. The resulting output from Stage Two is a series of plots in PNG format showing the how the selection for the best models was accounted and the results; also, a document-topic CSV file will be saved. Normally, this is what the output should be from both, the NMF and LDA processes:

1. NMF:
    1. best_component.png
	2. nmf_error.png
	3. coe_matrix.png
	4. nmf_topic_result.png
	5. nmf_topic.csv (document-topic file)
2. LDA:
    1. mean_subclonal_500_samples.png
	2. optimal_lda_model.png
	3. heatmat_doc_topics.png
	4. random_topics_results.png
	5. document_topic.csv (document-topic file)
	
Both of the CSV files for the topics will have the format with the colums:  

      |   -   | Topic 1 | Topic 2 | ... | Topic K |
	  | Doc 1 | Val 1.1 | Val 1.2 | ... | Val 1.K |
	  | Doc 2 | Val 2.1 | Val 2.2 | ... | Val 2.K |
	  |  ...  |   ...   | ...     | ... |   ...   |
	  | Doc N | Val N.1 | Val N.2 | ... | Val N.K |

Our dependencies
-----------
We supplied a `docker/` folder containing all the installation instructions for the numerous dependencies we use throughout our pipeline in a form of a Docker file, and a test script which runs it and tests our pipeline with small parameter values. Therefore the only true dependecy the user should install before using the pipeline is Docker (https://www.docker.com/) and a way to run our `test.sh` file, should the user's operating system (most notably Windows) not have such a native capability.

#### IMPORTANT NOTE:
On Windows System it can be run `test_windows.sh`. Nonetheless, please first read carefully the instructions inside the README file in the `docker/` folder.

Running the pipeline
-----------
Running `python pipeline.py --test` inside the `Pipeline/` directory provides a useful summary for the pipeline parameters:  

  
    usage: pipeline.py [-h] [--path WORKPLACE] [--random-samples RANDOM_SAMPLES]
                       [--samples-to-run SELECTED_SAMPLES]
                       [--phylowgs-burnin-samples PHWGS_BURNIN_SAMPLES]
                       [--phylowgs-mcmc-samples PHWGS_MCMC_SAMPLES]
                       [--phylowgs-mh-iterations PHWGS_MH_ITERATIONS]
                       [--pyclone-burnin-samples PYCLONE_BURNIN_SAMPLES]
                       [--pyclone-mcmc-iterations PYCLONE_MCMC_ITERATIONS]
                       [--ccube-max-clusters CCUBE_MAXCLUSTER]
                       [--ccube-vbem-max-iters CCUBE_VBEM_MAX_ITERS]
                       [--ccube-repeats CCUBE_REPEAT]
                       [--ccube-random-seed RANDOM_SEED] [--num-cores NUM_CORES]
                       [--run-pyclone] [--run-ccube] [--run-phylowgs]
                       [--pyclone-delete-tmp-folder]
                       [--lda-nmf-input RESULTS_FOLDER] [--run-nmf] [--run-lda]
    
  
    Pipeline
        
    optional arguments:
      -h, --help            show this help message and exit
      --path WORKPLACE      Specifies working directory for analysis. All paths in
                            the rest of the PyClone and Ccube files are relative
                            to this
      --random-samples RANDOM_SAMPLES
                            Number of randomnly selected samples to run
      --samples-to-run SELECTED_SAMPLES
                            A newline-seperated file of explicitly-stated tumour
                            sample names to be tested with Stage One tools
                            (PyClone, Ccube, PhyloWGS
      --phylowgs-burnin-samples PHWGS_BURNIN_SAMPLES
                            Number of burn-in samples for PhyloWGS (default: 1000)
      --phylowgs-mcmc-samples PHWGS_MCMC_SAMPLES
                            Number of true MCMC samples for PhyloWGS (default:
                            2500)
      --phylowgs-mh-iterations PHWGS_MH_ITERATIONS
                            Number of Metropolis-Hastings iterations for PhyloWGS
                            (default: 5000)
      --pyclone-burnin-samples PYCLONE_BURNIN_SAMPLES
                            Number of burn-in samples for PyClone (10 percent of
                            total MCMC suggested in the official documentation, 50
                            used here)
      --pyclone-mcmc-iterations PYCLONE_MCMC_ITERATIONS
                            Number of MCMC iterations for PyClone (default: 500)
      --ccube-max-clusters CCUBE_MAXCLUSTER
                            Maximum number of clusters for Ccube (default: 6)
      --ccube-vbem-max-iters CCUBE_VBEM_MAX_ITERS
                            Number of VBEM iterations for Ccube (default: 1000)
      --ccube-repeats CCUBE_REPEAT
                            Number of repeated Ccube runs for each candidate
                            number of clusters (default: 1)
      --ccube-random-seed RANDOM_SEED
                            Random seed (used by Ccube), required to run the the
                            tool deterministically
      --num-cores NUM_CORES
                            Number of processor cores to be employed in
                            computations concurrently (used by Ccube and PhyloWGS)
                            (default: 1)
      --run-pyclone         Flag for running PyClone (default: False)
      --run-ccube           Flag for running Ccube (default: False)
      --run-phylowgs        Flag for running PhyloWGS (default: False)
      --pyclone-delete-tmp-folder
                            Flag to delete temporary folder ("./inputTmp")
                            containing configuration files generated while
                            performing PyClone runs (default: False)
      --lda-nmf-input RESULTS_FOLDER
                            Input folder for LDA and/or NMF analysis (may be one
                            of the result folders of Stage One tools)
      --run-nmf             Flag for running NMF (default: False)
      --run-lda             Flag for running LDA (default: False)

Therefore it is possible to launch each of the Stage One tools seperately. For example:  
  
      python pipeline.py --run-phylowgs --num-cores 4 --phylowgs-burnin-samples 2 --phylowgs-mcmc-samples 20 --phylowgs-mh-iterations 50  
      python pipeline.py --run-pyclone --pyclone-burnin-samples 2 --pyclone-mcmc-iterations 20  
      python pipeline.py --run-ccube --ccube-max-clusters 6 --ccube-repeats 1 --ccube-vbem-max-iters 1000
          
Or, all of them together:  
  
      python pipeline.py --run-pyclone --run-ccube --run-phylowgs --pyclone-burnin-samples 2 --pyclone-mcmc-iterations 20  --ccube-max-clusters 6 --ccube-repeats 1 --ccube-vbem-max-iters 1000 --phylowgs-burnin-samples 2 --phylowgs-mcmc-samples 20 --phylowgs-mh-iterations 50
      
If a user wants to use the pipeline for a specific subset of the samples, he can create a newline-seperated file, containing the sample names he wishes to use, for example `samples.to.run.tsv` containing:  
      
      Sample_500_22  
      Sample_500_23  
      Sample_500_24  
      Sample_500_25  
      ...
  
And supplying the file as a `--samples-to-run` parameter:  
  
      python pipeline.py --run-pyclone --run-ccube --run-phylowgs --pyclone-burnin-samples 2 --pyclone-mcmc-iterations 20  --ccube-max-clusters 6  
      --ccube-repeats 1 --ccube-vbem-max-iters 1000 --phylowgs-burnin-samples 2 --phylowgs-mcmc-samples 20 --phylowgs-mh-iterations 50 --samples-to-run ./samples.to.run.tsv
    
This perhaps useful for such purposes of resuming an interrupting lenghty run. The sample names supplied with `samples.to.run.tsv` must match those described in `data/pp_table.txt` and have their corresponding input data named in a same pattern and located in the same folder structre as described in the "Our data" section.

If all of Stage One tools are run as a single command, regardless of their order in the command line, PyClone will be run first, followed by Ccube, followed by PhyloWGS. The results of these commands are stored in `Pipeline/results/pyclone`, `Pipeline/results/ccube` and `Pipeline/results/phylowgs` respectively.

Stage Two algorithms are launched in a similar manner. If the user already has the output from one of the Stage One tools, it is possible to launch Stage Two as a seperate command:    

      python pipeline.py --run-lda --run-nmf --lda-nmf-input /results/ccube/
    
Or any other source, perhaps not generated by the tools implemented in Stage One, assuming their output is following the common format:    

      python pipeline.py --run-lda --run-nmf --lda-nmf-input /some/other/input

#### IMPORTANT NOTE:  
We supplied the submission with previously run 500 ccube results and Stage Two can be started by executing the pipeline with the following line:

	  python pipeline.py --run-lda --run-nmf --lda-nmf-input /data/Data/Phase_two/Subclonal_Structure/

Licence
-----------
The usage of the source code and input data is entirely subject to University of Glasgow and Dr Ke Yuan of the School of Science and Engineering
