version development

task GatherFilesForMultiqc {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Array[File] inp_files
    Array[File] inp_files2
    String? output_dir
  }
  command <<<
    set -e
     \
      mkdir \
      ~{select_first([output_dir, "output_dir"])} \
      ; \
      cp \
      ~{if length(inp_files) > 0 then "'" + sep("' '", inp_files) + "'" else ""} \
      ~{if length(inp_files2) > 0 then "'" + sep("' '", inp_files2) + "'" else ""} \
      ~{if defined(select_first([output_dir, "output_dir"])) then ("'" + select_first([output_dir, "output_dir"]) + "'") else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
    zones: "australia-southeast1-b"
  }
  output {
    Directory out = select_first([output_dir, "output_dir"])
  }
}