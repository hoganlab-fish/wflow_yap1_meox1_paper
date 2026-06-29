version development

task Velocyto_Run10x {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Directory sample_folder
    File gtf_file
    Int? samtools_threads
    Int? samtools_memory
    Int? dtype
  }
  command <<<
    set -e
     \
      ulimit -n 51200; \
      export LC_ALL=C.UTF-8; export LANG=C.UTF-8; \
      velocyto run10x \
      ~{sample_folder} \
      ~{gtf_file} \
      --samtools-threads ~{select_first([samtools_threads, select_first([runtime_cpu, 1])])} \
      --samtools-memory ~{select_first([samtools_memory, select_first([runtime_memory, 4])])} \
      --dtype ~{select_first([dtype, "uint32"])} \
      ;cp -r ~{sample_folder} .;
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/velocyto.py@sha256:b1772569196709c44c643558dd432f4829f14795e07899a45c11a82d26f9c8f5"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
    zones: "australia-southeast1-b"
  }
  output {
    Directory outdir = sample_folder
  }
}