task getScript {
	command {
		
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/dev/workflows/aggregateAssociation/regional_plot.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File plot_script = "regional_plot.R"
	}
}

task plot_top {
		File results_file
		File group_file
        File state_file
        File gene_file
        String out_pref
        File script

        command {
                R --vanilla --args ${results_file} ${group_file} ${state_file} ${gene_file} ${out_pref} < ${script}
        }

        meta {
                author: "jasen jackson"
                email: "jasenjackson97@gmail.com"
        }

        runtime {
    	   docker: "robbyjo/r-mkl-bioconductor@sha256:b88d8713824e82ed4ae2a0097778e7750d120f2e696a2ddffe57295d531ec8b2"
		   disks: "local-disk 100 SSD"
		   memory: "3G"
        }

        output {
                File plots = "top_hits_${out_pref}.pdf"
        }
}


workflow w_assocTest {
	File this_results_file
	Array[Pair[File,File]] these_gds_groups
    File this_state_file
    File this_gene_file
    
    Array[String] these_out_pref = range(length(these_gds_groups))
    Array[Pair[Int,Pair[File,File]]] ind_gds_pair = zip(these_out_pref,these_gds_groups)

    
	
	call getScript

	scatter(this_pair in ind_gds_pair) {
		
		Pair[File,File] gds_group = this_pair.right

		call plot_top {
			input: results_file=this_results_file, group_file=gds_group.right, state_file=this_state_file, gene_file=this_gene_file, out_pref=this_pair.left, script=getScript.plot_script
		}

	}

}