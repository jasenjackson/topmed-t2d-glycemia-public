# adding a task that takes the vcf file and creates a bgzipped vcf and a tabix index file of the vcf.
# can be optional output

task getScript {
	command {
		wget "https://raw.githubusercontent.com/manning-lab/topmed-t2d-glycemia-public/master/methods/fileConvert/vcfToGds.R"
	}

	runtime {
		docker: "tmajarian/alpine_wget@sha256:f3402d7cb7c5ea864044b91cfbdea20ebe98fc1536292be657e05056dbe5e3a4"
	}

	output {
		File outscript = "vcfToGds.R"
	}
}

task vcfToGds {
	File in
	File indat
	String in_base = basename(indat,".vcf")
	Int disksize
	Float memory
	
	command {
		R --vanilla --args ${indat} ${in_base} < ${in}
	}
	runtime {
		docker: "robbyjo/r-mkl-bioconductor@sha256:96d6b7b1ab1ec89d6413c09f7fef3fdaac5e9d4293b259492ddda4cf5883d354"
		disks: "local-disk ${disksize} SSD"
		memory: "${memory}G"
	}
	output { 
		#File out = read_string("output.txt")
		String out = "${in_base}" + ".gds"
	}
}

task vcfToBgvcf {
	File vcf_in
	Int diskSize
	Float Memory

	command {
		bgzip -c ${vcf_in} > ${vcf_in}.gz
		tabix -p vcf ${vcf_in}.gz
	}

	runtime {
		docker: "biowardrobe2/samtools@sha256:e4dad5f38c1b782d3f1608410c07e8dc47fb7b92bc427175a160dfa0813c48d8"
		disks: "local-disk ${diskSize} SSD"
		memory: "${Memory}G"
	}

	output {
		String zip = "${vcf_in}" + ".gz"
		String ind = "${vcf_in}" + ".gz.tbi"
		# File zip = "${vcf_in}.gz"
		# File ind = "${vcf_in}.gz.tbi"
	}
}

workflow w {
	Array[File] infiles
	Int thisDisksize
	Float thisMemory
	Boolean gzip
	Boolean gds

	scatter(this_file in infiles) {

		if(gds) {
			call getScript
			call vcfToGds {input: in=getScript.outscript, indat=this_file, disksize=thisDisksize, memory=thisMemory}
		}

		if(gzip) {
			call vcfToBgvcf {input: vcf_in=this_file, diskSize=thisDisksize, Memory=thisMemory}
		}
	}

	output {
		Array[File]? gdsOut = select_all(vcfToGds.out)
		Array[File]? vcfzip = select_all(vcfToBgvcf.zip)
		Array[File]? vcfind = select_all(vcfToBgvcf.ind)
	}
}
