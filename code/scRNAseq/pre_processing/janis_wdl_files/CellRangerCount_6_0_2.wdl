version development

task CellRangerCount {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    String? cellranger_version
    String id
    Array[Directory] fastqs
    Directory transcriptome
    Array[String] sample
    Int? localcores
    Int? localmem
  }
  command <<<
    set -e
     \
      . /etc/profile.d/modules.sh; module load cellranger/~{select_first([cellranger_version, "6.0.2"])}; \
      cellranger count \
      --id ~{id} \
      ~{if length(fastqs) > 0 then "--fastqs " + sep(",", fastqs) else ""} \
      --transcriptome ~{transcriptome} \
      ~{if length(sample) > 0 then "--sample " + sep(",", sample) else ""} \
      ~{if defined(select_first([localcores, select_first([runtime_cpu, 1])])) then ("--localcores " + select_first([localcores, select_first([runtime_cpu, 1])])) else ''} \
      ~{if defined(select_first([localmem, select_first([runtime_memory, 4])])) then ("--localmem " + select_first([localmem, select_first([runtime_memory, 4])])) else ''}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
    zones: "australia-southeast1-b"
  }
  output {
    Directory outdir = id
    File out_matrix = (id + "/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
  }
}