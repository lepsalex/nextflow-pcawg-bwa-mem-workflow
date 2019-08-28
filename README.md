# nextflow-pcawg-bwa-mem-workflow
Nextflow implementation of https://github.com/ICGC-TCGA-PanCancer/Seqware-BWA-Workflow

### Example Command
`nextflow run workflow.nf -params-file params.json`

#### Kubernetes
`nextflow kuberun -params-file params.k8s.json lepsalex/nextflow-pcawg-bwa-mem-workflow -latest -v nextflow-pv-claim:/mnt/volume/nextflow`