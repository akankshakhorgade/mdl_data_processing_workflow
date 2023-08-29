version 1.0

task pbSkera {
    meta {
        description: "Given hifi reads, spilts MAS 8x or 16x array structure using provided adapters"
    }
# ------------------------------------------------
#Inputs required
    input {
        # Required:
        File hifi_bam
        String sample_id
        File mas_adapters_fasta
        Int num_threads
        String gcs_output_dir

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 2.5*(size(hifi_bam, "GiB"))
    Int default_ram = 8
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 25

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    command <<<
        set -euxo pipefail

        gsutil -m cp ~{hifi_bam} .
        touch read_counts.txt
        echo ~{hifi_bam} >> read_counts.txt
        samtools view -c ~{sample_id}.hifi_reads.bam >> read_counts.txt
        echo ~{outdir}skera/~{sample_id}.skera.bam >> read_counts.txt
        gsutil cp gs://mdl_terra_sandbox/tools/skera /usr/local/bin/
        chmod 777 /usr/local/bin/skera
        /usr/local/bin/skera split -j ~{num_threads} ~{hifi_bam} ~{mas_adapters_fasta} ~{sample_id}.skera.bam
        echo "Copying skera out to gcs path provided..."
        gsutil -m cp ~{sample_id}.skera.* ~{outdir}skera/
        samtools view -c ~{sample_id}.skera.bam >> read_counts.txt
        gsutil -m cp read_counts.txt ~{outdir}skera/

    >>>
# ------------------------------------------------
# Outputs:
    output {
        # Default output file name:
        String skera_out        = "~{outdir}skera/~{sample_id}.skera.bam"
    }

# ------------------------------------------------
# Runtime settings:
    runtime {
    docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:tag3"
    memory: machine_mem + " GiB"
    disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
    preemptible: select_first([preemptible_attempts, 0])
    cpu: select_first([cpu, 2])
    }

}

task pbLima {
    meta {
        description: "Single Cell workflow: given s-reads, lima identifies and corrects barcodes and adapter sequences"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File skera_bam
        String sample_id
        File sc_primers
        Int num_threads
        String gcs_output_dir
        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 2.5*(size(skera_bam, "GiB"))
    Int default_ram = 8
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 25
    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    command <<<
        set -euxo pipefail
        gsutil -m cp ~{skera_bam} .
        echo "lima"
        lima --isoseq -j ~{num_threads} \
        --log-file ~{sample_id}.lima.log \
        ~{sample_id}.skera.bam ~{sc_primers} ~{sample_id}.lima.bam

        gsutil -m cp ~{sample_id}.lima.bam ~{outdir}lima/~{sample_id}.lima.bam
    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String lima_out        = "~{outdir}lima/~{sample_id}.lima.bam"
    }
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:sc1"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 2])
    }
}

task pbIsoseq3 {
    meta {
        description: "Single Cell workflow: isoseq tag, refine and correct modules"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File input_bam
        String sample_id
        File isoseq_design
        File sc_primers
        File ref_barcodes
        Int num_threads
        String gcs_output_dir
        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 1.5*(size(input_bam, "GiB"))
    Int default_ram = 16
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 25
    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "gs:/", "gs://")
    command <<<
        set -euxo pipefail
        gsutil -m cp ~{input_bam} .

        echo "tag"
        isoseq3 tag --design ~{isoseq_design} -j  ~{num_threads} \
        ~{input_bam} \
        ~{sample_id}.tagged.bam

        echo "refine"
        isoseq3 refine --require-polya -j ~{num_threads} \
        ~{sample_id}.tagged.bam \
        ~{sc_primers} \
        ~{sample_id}.refined.bam

        echo "correct"
        isoseq3 correct --barcodes ~{ref_barcodes} -j ~{num_threads} \
        ~{sample_id}.refined.bam \
        ~{sample_id}.corrected.bam

        gsutil -m cp ~{sample_id}.tagged.bam ~{outdir}isoseq/~{sample_id}.tagged.bam
        gsutil -m cp ~{sample_id}.refined.bam ~{outdir}isoseq/~{sample_id}.refined.bam
        gsutil -m cp ~{sample_id}.corrected.bam ~{outdir}isoseq/~{sample_id}.corrected.bam

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String lima_out        = "~{outdir}lima/~{sample_id}.lima.bam"
    }
    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us-east4-docker.pkg.dev/methods-dev-lab/masseq-dataproc/masseq_prod:sc1"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 2])
    }
}

