version development

task fastqc {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[File] reads
    File? read1
    File? read2
    String? outdir
    Boolean? casava
    Boolean? nano
    Boolean? nofilter
    Boolean? extract
    String? java
    Boolean? noextract
    Boolean? nogroup
    String? format
    Int? threads
    File? contaminants
    File? adapters
    File? limits
    Int? kmers
    Boolean? quiet
    String? dir
  }
  command <<<
    set -e
    fastqc \
      ~{if defined(select_first([outdir, "."])) then ("--outdir '" + select_first([outdir, "."]) + "'") else ""} \
      ~{if (defined(casava) && select_first([casava])) then "--casava" else ""} \
      ~{if (defined(nano) && select_first([nano])) then "--nano" else ""} \
      ~{if (defined(nofilter) && select_first([nofilter])) then "--nofilter" else ""} \
      ~{if select_first([extract, true]) then "--extract" else ""} \
      ~{if defined(java) then ("--java '" + java + "'") else ""} \
      ~{if (defined(noextract) && select_first([noextract])) then "--noextract" else ""} \
      ~{if (defined(nogroup) && select_first([nogroup])) then "--nogroup" else ""} \
      ~{if defined(format) then ("--format '" + format + "'") else ""} \
      ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("--threads " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
      ~{if defined(contaminants) then ("--contaminants '" + contaminants + "'") else ""} \
      ~{if defined(adapters) then ("--adapters '" + adapters + "'") else ""} \
      ~{if defined(limits) then ("--limits '" + limits + "'") else ""} \
      ~{if defined(kmers) then ("--kmers " + kmers) else ''} \
      ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
      ~{if defined(dir) then ("--dir '" + dir + "'") else ""} \
      ~{if defined(select_first([read1, reads[0]])) then ("'" + select_first([read1, reads[0]]) + "'") else ""} \
      ~{if defined(select_first([read2, reads[1]])) then ("'" + select_first([read2, reads[1]]) + "'") else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/fastqc@sha256:810db4a1676d79883cc5f29ce31b4db0b6ebed36007fae5fdd1e78a5304639d8"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
    zones: "australia-southeast1-b"
  }
  output {
    File out_R1 = (basename(basename(select_first([read1, reads[0]]), ".fastq.gz"), ".fq.gz") + "_fastqc.zip")
    File out_R1_datafile = (basename(basename(select_first([read1, reads[0]]), ".fastq.gz"), ".fq.gz") + "_fastqc/fastqc_data.txt")
    File out_R1_html = (basename(basename(select_first([read1, reads[0]]), ".fastq.gz"), ".fq.gz") + "_fastqc.html")
    Directory out_R1_directory = (basename(basename(select_first([read1, reads[0]]), ".fastq.gz"), ".fq.gz") + "_fastqc")
    File out_R2 = (basename(basename(select_first([read2, reads[1]]), ".fastq.gz"), ".fq.gz") + "_fastqc.zip")
    File out_R2_datafile = (basename(basename(select_first([read2, reads[1]]), ".fastq.gz"), ".fq.gz") + "_fastqc/fastqc_data.txt")
    File out_R2_html = (basename(basename(select_first([read2, reads[1]]), ".fastq.gz"), ".fq.gz") + "_fastqc.html")
    Directory out_R2_directory = (basename(basename(select_first([read2, reads[1]]), ".fastq.gz"), ".fq.gz") + "_fastqc")
  }
}