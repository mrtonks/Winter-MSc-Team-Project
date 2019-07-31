load("inputTmp/myfile")
yamlScript <- cat("num_iters: ",iteration, "\n",
    "base_measure_params:", "\n",
    "  alpha: 1", "\n",
    "  beta: 1", "\n",
    "concentration:", "\n",
    "  value: 1.0", "\n",
    "  prior:", "\n",
    "    shape: 1.0", "\n",
    "    rate: 0.001", "\n",
    "density: pyclone_beta_binomial", "\n",
    "beta_binomial_precision_params:", "\n",
    "  value: 1000", "\n",
    "  prior:", "\n",
    "    shape: 1.0", "\n",
    "    rate: 0.0001", "\n",
    "  proposal:", "\n",
    "    precision: 0.01", "\n",
    paste0('working_dir: ', pycloneFolder), "\n",
    "trace_dir: trace", "\n",
    "samples:", "\n",
    paste0("  ", sampleName, ":"), "\n",
    paste0("    mutations_file: ", pycloneFolder, "/pyclone_mutations.yaml"), "\n",
    "    tumour_content:", "\n",
    paste0("      value: ", cellularity), "\n",
    "    error_rate: 0.001",
    sep = "",
    file = paste0(pycloneFolder, "/pyclone_configure.yaml")
)
traceFolder = paste0(pycloneFolder,"/trace")
dir.create(traceFolder)
save.image("inputTmp/myfile")
