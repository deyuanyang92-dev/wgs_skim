# Histone seed / label for GetOrganelle (-F anonym)

- `seed.fasta`, `label.fasta`: 69 H1 / H2A / H2B / H3 / H4 copies of the histone tandem cluster from six annelid
  genome skims (MJ20260820041 CHY002, CHY006-010: `mitom nuc` on SPAdes scaffolds, contigs >= 100x coverage =
  cluster, >= 150 nt; names `<gene>|<sample>|cov<x>|<contig>`) + the user's Hesionidae H3 seed (6 sequences) and
  label (87; `/media/cody/18T_R1/ZS/Hesionidae_H3genes`, NCBI partial H3 CDS, 2024).
- Why all five genes: with an H3-only label GetOrganelle kept only the 2.2-kb H3 piece of the cluster (CHY006);
  with H1-H4 it kept all five genes.
- Another taxon: give your own with `histone_seed` / `histone_label` (config) or `--seed` / `--label`.
