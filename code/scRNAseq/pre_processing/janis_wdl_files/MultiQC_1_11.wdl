version development

task MultiQC {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Directory directory
    Boolean? force
    String? dirs
    Int? dirsDepth
    Boolean? fullnames
    String? title
    String? comment
    String? filename
    String? outdir
    String? template
    String? tag
    Boolean? view_tags
    Boolean? ignore
    Boolean? ignoreSamples
    Boolean? ignoreSymlinks
    File? sampleNames
    Array[String]? exclude
    Array[String]? module
    Boolean? dataDir
    Boolean? noDataDir
    String? dataFormat
    Boolean? export
    Boolean? flat
    Boolean? interactive
    Boolean? lint
    Boolean? pdf
    Boolean? noMegaqcUpload
    File? config
    File? cl_config
    Boolean? verbose
    Boolean? quiet
  }
  command <<<
    set -e
    multiqc \
      ~{if (defined(force) && select_first([force])) then "--force" else ""} \
      ~{if defined(dirs) then ("--dirs '" + dirs + "'") else ""} \
      ~{if defined(dirsDepth) then ("--dirs-depth " + dirsDepth) else ''} \
      ~{if (defined(fullnames) && select_first([fullnames])) then "--fullnames" else ""} \
      ~{if defined(title) then ("--title '" + title + "'") else ""} \
      ~{if defined(comment) then ("--comment '" + comment + "'") else ""} \
      --filename '~{select_first([filename, "generated"])}' \
      ~{if defined(select_first([outdir, "."])) then ("--outdir '" + select_first([outdir, "."]) + "'") else ""} \
      ~{if defined(template) then ("--template '" + template + "'") else ""} \
      ~{if defined(tag) then ("--tag '" + tag + "'") else ""} \
      ~{if (defined(view_tags) && select_first([view_tags])) then "--view_tags" else ""} \
      ~{if (defined(ignore) && select_first([ignore])) then "--ignore" else ""} \
      ~{if (defined(ignoreSamples) && select_first([ignoreSamples])) then "--ignore-samples" else ""} \
      ~{if (defined(ignoreSymlinks) && select_first([ignoreSymlinks])) then "--ignore-symlinks" else ""} \
      ~{if defined(sampleNames) then ("--sample-names '" + sampleNames + "'") else ""} \
      ~{if (defined(exclude) && length(select_first([exclude])) > 0) then "--exclude '" + sep("' --exclude '", select_first([exclude])) + "'" else ""} \
      ~{if (defined(module) && length(select_first([module])) > 0) then "--module '" + sep("' --module '", select_first([module])) + "'" else ""} \
      ~{if (defined(dataDir) && select_first([dataDir])) then "--data-dir" else ""} \
      ~{if (defined(noDataDir) && select_first([noDataDir])) then "--no-data-dir" else ""} \
      ~{if defined(dataFormat) then ("--data-format '" + dataFormat + "'") else ""} \
      ~{if (defined(export) && select_first([export])) then "--export" else ""} \
      ~{if (defined(flat) && select_first([flat])) then "--flat" else ""} \
      ~{if (defined(interactive) && select_first([interactive])) then "--interactive" else ""} \
      ~{if (defined(lint) && select_first([lint])) then "--lint" else ""} \
      ~{if (defined(pdf) && select_first([pdf])) then "--pdf" else ""} \
      ~{if (defined(noMegaqcUpload) && select_first([noMegaqcUpload])) then "--no-megaqc-upload" else ""} \
      ~{if defined(config) then ("--config '" + config + "'") else ""} \
      ~{if defined(cl_config) then ("--cl_config '" + cl_config + "'") else ""} \
      ~{if (defined(verbose) && select_first([verbose])) then "--verbose" else ""} \
      ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
      '~{directory}'
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/multiqc@sha256:88df23fac5b9eecda9943d922f81b68e30188eb4dd7cbfe9554e952ff5a3b0ee"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 4])}G"
    preemptible: 2
    zones: "australia-southeast1-b"
  }
  output {
    File out_html = "~{select_first([outdir, "."])}/~{select_first([filename, "generated"])}.html"
    Directory out_multiqc_data = "~{select_first([outdir, "."])}/~{select_first([filename, "generated"])}_data"
  }
}