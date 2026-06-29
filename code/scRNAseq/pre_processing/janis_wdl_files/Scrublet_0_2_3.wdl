version development

task Scrublet {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    File scrublet_script
    File matrix
    String hist
    String umap
    String score
    String pred
  }
  command <<<
    set -e
    python \
      ~{scrublet_script} \
      --matrix ~{matrix} \
      --hist ~{hist} \
      --umap ~{umap} \
      --score ~{score} \
      --pred ~{pred}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/scrublet@sha256:2a4a79a2c72725db6560f3c9c823fd1128dd81b811816405cb3b752f83fe3421"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
    zones: "australia-southeast1-b"
  }
  output {
    File out_hist = hist
    File out_umap = umap
    File out_score = score
    File out_pred = pred
  }
}