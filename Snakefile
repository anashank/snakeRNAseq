configfile: "config.yaml"

if config["map_method"] == "hisat2" or config["map_method"] == "bowtie2":
  rule all:
    input:
      expand("sorted/{sample}.sorted.bam.bai",sample=config["samples"])

elif config["map_method"] == "kallisto":
  rule all:
    input:
      expand("output_kallisto/{sample}",sample=config["samples"])

if config["end_type"] == "pe":
  rule trim_fastqc_pe:
      input:
          read1="input/{sample}_1.fq.gz",
          read2="input/{sample}_2.fq.gz"
      output:
          read1_val="trimmed/{sample}_1_val_1.fq.gz",
	  read2_val="trimmed/{sample}_2_val_2.fq.gz",
          report1out = "trimmed/{sample}_1.fq.gz_trimming_report.txt",
          report2out = "trimmed/{sample}_2.fq.gz_trimming_report.txt"
      params:
          p1 = directory("trimmed")
      log:
          trimgalore = "log/{sample}.trimgalore"

      shell:
          "trim_galore -q {config[trim_quality]} -j {config[trim_cores]} --fastqc_args \"--outdir ./fastQC_output/\" {config[trim_params]} --paired -o {params.p1} {input.read1} {input.read2} &> {log.trimgalore}"


  if config["map_method"] == "bowtie2":
    rule mapReads_bowtie2_pe:
        input:
            read1_val="trimmed/{sample}_1_val_1.fq.gz",
  	        read2_val="trimmed/{sample}_2_val_2.fq.gz"
        output:
            "mappedReads_bowtie2/{sample}.sam"
        log:
            bowtie = "log/{sample}.bowtie"
        shell:
            "bowtie2 -x {config[bowtie_idx]} -1 {input.read1_val} -2 {input.read2_val} {config[bowtie_params]} -S {output}  &> {log.bowtie}"


    rule sort_bam:
       input:
           "mappedReads_bowtie2/{sample}.sam"
       output:
           "sorted/{sample}.sorted.bam"
       log:
           samtools = "log/{sample}.sort"
       shell:
          "samtools view -Sb -F 4 {input}| samtools sort {config[sort_params]} -o {output} &> {log.samtools}"

    #Index bam files
    rule index_bam:
        input:  
            "sorted/{sample}.sorted.bam"
        output: 
            "sorted/{sample}.sorted.bam.bai"
        log:    
            index = "log/{sample}.index_bam"
        shell:
            "samtools index {input} &> {log.index}"


  elif config["map_method"] == "hisat2":

    rule mapReads_hisat2_pe:
      input:
          read1_val="trimmed/{sample}_1_val_1.fq.gz",
          read2_val="trimmed/{sample}_2_val_2.fq.gz"
      output:
          "mappedReads_hisat2/{sample}.sam"
      log:
          hisat = "log/{sample}.hisat"
      shell:
          "hisat2 -x {config[hisat_idx]} -1 {input.read1_val} -2 {input.read2_val} {config[hisat_params]} -S {output}  &> {log.hisat}"

    rule sort_bam:
       input:
           "mappedReads_hisat2/{sample}.sam"
       output:
           "sorted/{sample}.sorted.bam"
       log:
           samtools = "log/{sample}.sort"
       shell:
          "samtools view -Sb -F 4 {input}| samtools sort {config[sort_params]} -o {output} &> {log.samtools}"

    rule index_bam:
        input:  
            "sorted/{sample}.sorted.bam"
        output: 
            "sorted/{sample}.sorted.bam.bai"
        log:    
            index = "log/{sample}.index_bam"
        shell:
            "samtools index {input} &> {log.index}"



  elif config["map_method"] == "kallisto":
    rule kallisto_pe:
      input:
          read1_val="trimmed/{sample}_1_val_1.fq.gz",
          read2_val="trimmed/{sample}_2_val_2.fq.gz"
      output:
          directory("output_kallisto/{sample}")
      log:
          kallisto = "log/{sample}.kallisto"
      shell:
          "kallisto quant -i {config[kallisto_idx]} {config[kallisto_params]} -o {output} {input.read1_val} {input.read2_val}  &> {log.kallisto}"


elif config["end_type"] == "se":
  rule trim_fastqc_se:
      input:
          read="input/{sample}.fq.gz",
      output:
          read_val="trimmed/{sample}_trimmed.fq.gz",
          report = "trimmed/{sample}.fq.gz_trimming_report.txt",
      params:
          p1 = directory("trimmed")
      log:
          trimgalore = "log/{sample}.trimgalore"

      shell:
          "trim_galore -q {config[trim_quality]} -j {config[trim_cores]} --fastqc_args \"--outdir ./fastQC_output/\" {config[trim_params]} -o {params.p1} {input.read} &> {log.trimgalore}"
  
  if config["map_method"] == "bowtie2":
    rule mapReads_bowtie2_se:
        input:
            read_val="trimmed/{sample}_trimmed.fq.gz"
        output:
            "mappedReads_bowtie2/{sample}.sam"
        log:
            bowtie = "log/{sample}.bowtie"
        shell:
            "bowtie2 -x {config[bowtie_idx]} -U {input.read_val} {config[bowtie_params]} -S {output}  &> {log.bowtie}"

    rule sort_bam:
       input:
           "mappedReads_bowtie2/{sample}.sam"
       output:
           "sorted/{sample}.sorted.bam"
       log:
           samtools = "log/{sample}.sort"
       shell:
          "samtools view -Sb -F 4 {input}| samtools sort {config[sort_params]} -o {output} &> {log.samtools}"

    rule index_bam:
        input:  
            "sorted/{sample}.sorted.bam"
        output: 
            "sorted/{sample}.sorted.bam.bai"
        log:    
            index = "log/{sample}.index_bam"
        shell:
            "samtools index {input} &> {log.index}"

  elif config["map_method"] == "hisat2":

    rule mapReads_hisat2_se:
      input:
          read_val="trimmed/{sample}_trimmed.fq.gz"
      output:
          "mappedReads_hisat2/{sample}.sam"
      log:
          hisat = "log/{sample}.hisat"
      shell:
          "hisat2 -x {config[hisat_idx]} -U {input.read_val} {config[hisat_params]} -S {output}  &> {log.hisat}"

    rule sort_bam:
       input:
           "mappedReads_hisat2/{sample}.sam"
       output:
           "sorted/{sample}.sorted.bam"
       log:
           samtools = "log/{sample}.sort"
       shell:
          "samtools view -Sb -F 4 {input}| samtools sort {config[sort_params]} -o {output} &> {log.samtools}"

    rule index_bam:
        input:  
            "sorted/{sample}.sorted.bam"
        output: 
            "sorted/{sample}.sorted.bam.bai"
        log:    
            index = "log/{sample}.index_bam"
        shell:
            "samtools index {input} &> {log.index}"


  elif config["map_method"] == "kallisto":
    rule kallisto_se:
      input:
          read_val="trimmed/{sample}_trimmed.fq.gz"
      output:
          directory("output_kallisto/{sample}")
      log:
          kallisto = "log/{sample}.kallisto"
      shell:
          "kallisto quant -i {config[kallisto_idx]} {config[kallisto_params]} -o {output} --single -l {config[kallisto_fragment_length]} -s {config[kallisto_sd]} {input.read_val}  &> {log.kallisto}"






