task getScript {
	command {
		wget https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/fitNull/genesis_nullmodel.R
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File null_script = "genesis_nullmodel.R"
	}
}

task fitNull {
	File phenofile
	String outcomename
	String outcometype
	String covariates
	File sample_ids
	String label
	File kinshipmatrix
	String phenoid 

	File script

	command {
		R --vanilla --args ${phenofile} ${outcomename} ${outcometype} ${covariates} ${sample_ids} ${label} ${kinshipmatrix} ${phenoid} < ${script}
	}

	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk 20 SSD"
		memory: "20G"
	}

	output {
		File model = "${label}_null.RDa"
		File plots = "${label}_plots.pdf"
		
	command {
		wget https://raw.githubusercontent.com/jasenjackson/topmed-t2d-glycemia-public/dev/workflows/fitNull/genesis_nullmodel.R
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File null_script = "genesis_nullmodel.R"
	}
}

task fitNull {
	File phenofile
	String outcomename
	String outcometype
	String covariates
	File sample_ids
	String label
	File kinshipmatrix
	String phenoid 

	File script

	command {
		R --vanilla --args ${phenofile} ${outcomename} ${outcometype} ${covariates} ${sample_ids} ${label} ${kinshipmatrix} ${phenoid} < ${script}
	}

	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		disks: "local-disk 20 SSD"
		memory: "20G"
	}

	output {
		File model = "${label}_null.RDa"
		File plots = "${label}_plots.pdf"
		File stats = "${label}_stats.tsv"
	}
}


workflow nullModel {
	File this_phenofile
	String this_outcomename
	String this_outcometype
	String this_covariates
	File this_sample_ids
	String this_label
	File this_kinshipmatrix
	String this_phenoid

	call getScript
	
	call fitNull {
            input: phenofile=this_phenofile, outcomename=this_outcomename, outcometype=this_outcometype, covariates=this_covariates, sample_ids=this_sample_ids, label=this_label, kinshipmatrix=this_kinshipmatrix, phenoid=this_phenoid, script=getScript.null_script
	}
}
	}
}


workflow nullModel {
	File this_phenofile
	String this_outcomename
	String this_outcometype
	String this_covariates
	File this_sample_ids
	String this_label
	File this_kinshipmatrix
	String this_phenoid

	call getScript
	
	call fitNull {
            input: phenofile=this_phenofile, outcomename=this_outcomename, outcometype=this_outcometype, covariates=this_covariates, sample_ids=this_sample_ids, label=this_label, kinshipmatrix=this_kinshipmatrix, phenoid=this_phenoid, script=getScript.null_script
	}
}
