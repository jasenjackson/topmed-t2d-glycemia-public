task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/aggregateAssociation/groupByGene.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File group_script = "groupByGene.R"
	}
}

task get_groups {
	File gds
	File allGenes
	File panGenes
	File anno
	File state
	File chain
	File groupScript

	command {
        	R --vanilla --args ${gds} ${allGenes} ${panGenes} ${anno} ${state} ${chain} < ${groupScript}
    }

    meta {
            author: "TM"
            email: "tmajaria@broadinstitute.org"
    }

    runtime {
    	   docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		   disks: "local-disk 100 SSD"
		   memory: "30G"
    }

    output {
    	Pair[File,File] out_groups = (gds,"groups.RData")
        # File out_groups = "groups.RData"
    }	

}

workflow wf {
	Array[File] these_gds
	String this_label

	File this_allGenes
	File this_panGenes
	File this_anno
	File this_state
	File this_chain
	
	call getScript
	
	scatter(this_gds in these_gds) {
		call get_groups {
				input: gds=this_gds, allGenes=this_allGenes, panGenes=this_panGenes, anno=this_anno, state=this_state, chain=this_chain, groupScript=getScript.group_script
		}
	}
}