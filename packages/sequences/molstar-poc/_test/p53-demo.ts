/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Seq*/Cases/Mol* › p53: sequence and structure` document (seqspec:libraries#the-mol-adapter), shared by
// the story, its smoke entry and the compose tests: P04637 with its UniProt features and Proteins API variants
// (lib/seqstar/molstar-poc/_test/fixtures/p53, trimmed), a hydropathy profile and the residues in a letter scheme; five structures as
// complexes mapped onto it (by SIFTS, the AlphaFold model by identity); the PDB coverage in two levels (collapsed, its
// depth by `molstar:sifts_depth`; open, a row per structure by `molstar:sifts_coverage`, each split by chain), every
// coverage row naming its structure, its rows on the chains' own units (instance space), so a chain's residues light
// its own row alone; the model's pLDDT. The structure shown in 3D is marked on its entry (`p53ShownPatch`). No I/O: the
// hosts pass their own loaders (`p53Loaders`).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { Loaders, Patch, SeqViewSpecDoc } from "@molstar/seq-core/format";

/** The structures, in the coverage's order: their complex, source, label, what they cover. */
export const P53_STRUCTURES = [
  {
    complex: "1tup",
    source: { service: "pdbe", id: "1tup" },
    entity: "3",
    label: "1TUP",
    covers: "DNA-binding domain on DNA, 94–312",
  },
  {
    complex: "2ocj",
    source: { service: "pdbe", id: "2ocj" },
    entity: "1",
    label: "2OCJ",
    covers: "DNA-binding domain, 94–312",
  },
  {
    complex: "2xwr",
    source: { service: "pdbe", id: "2xwr" },
    entity: "1",
    label: "2XWR",
    covers: "DNA-binding domain, extended N terminus, 89–293",
  },
  {
    complex: "1c26",
    source: { service: "pdbe", id: "1c26" },
    entity: "1",
    label: "1C26",
    assembly: "1",
    covers: "tetramerization domain, 325–356 (its tetramer)",
  },
  {
    complex: "af-p04637-f1",
    source: { service: "afdb", id: "AF-P04637-F1" },
    entity: "1",
    label: "AlphaFold",
    covers: "AlphaFold model, full length, pLDDT",
  },
] as const;
export type P53Structure = (typeof P53_STRUCTURES)[number];

/** The track refs the story and the tests name. */
export const P53_TRACKS = {
  domains: "t-domains",
  regions: "t-regions",
  ptms: "t-ptms",
  sites: "t-sites",
  variants: "t-variants",
  hydropathy: "t-hydropathy",
  plddt: "t-plddt",
  coverage: "t-coverage",
  letters: "t-residues",
} as const;

const UNIPROT = "https://rest.uniprot.org/uniprotkb";
const VARIATION = "https://www.ebi.ac.uk/proteins/api/variation/P04637";

/** Kyte–Doolittle hydropathy over a centered window of 9 (positions without a whole window left out). */
const KD: Record<string, number> = {
  A: 1.8,
  R: -4.5,
  N: -3.5,
  D: -3.5,
  C: 2.5,
  Q: -3.5,
  E: -3.5,
  G: -0.4,
  H: -3.2,
  I: 4.5,
  L: 3.8,
  K: -3.9,
  M: 1.9,
  F: 2.8,
  P: -1.6,
  S: -0.8,
  T: -0.7,
  W: -0.9,
  Y: -1.3,
  V: 4.2,
};
export function hydropathy(sequence: string, window = 9): { index: number[]; value: number[] } {
  const h = window >> 1;
  const index: number[] = [];
  const value: number[] = [];
  for (let i = h; i < sequence.length - h; i++) {
    let s = 0;
    for (let k = i - h; k <= i + h; k++) s += KD[sequence[k]] ?? 0;
    index.push(i);
    value.push(Math.round((s / window) * 100) / 100);
  }
  return { index, value };
}

const byField = (field: string) => ({
  theme: "field_palette",
  params: { field, palette: { kind: "categorical", colors: "Set2" } },
});
const COVERAGE_COLORS = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"];
/** The PTM kinds' colors: distinct hues, gray only for "other". */
export const PTM_COLORS: Readonly<Record<string, string>> = {
  phosphorylation: "#e15759",
  acetylation: "#4e79a7",
  methylation: "#59a14f",
  ubiquitination: "#b07aa1",
  sumoylation: "#f28e2b",
  glycosylation: "#76b7b2",
  other: "#9d9d9d",
};
/** The host file the PTM sites are read from (the loaders serve it). */
export const PTMS_FILE = "P04637-ptms.tsv";

/** The demo's document over P04637's `sequence` (the profile is computed from it here). */
export function p53DemoDoc(sequence: string): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "p53: sequence and structure" });
  const p53 = b.data.sequence({ ref: "p53", source: { service: "uniprot", resource: "sequence", id: "P04637" } });
  const feats = b.data.annotation({
    ref: "features",
    geometry: "ranges",
    target: p53,
    source: { service: "uniprot", resource: "features", id: "P04637" },
  });
  const variants = b.data.annotation({
    ref: "variants",
    geometry: "ranges",
    target: p53,
    source: { service: "uniprot", resource: "variants", id: "P04637" },
  });
  const hydro = b.data.annotation({ ref: "hydropathy", geometry: "values", target: p53, value: hydropathy(sequence) });
  // The PTM sites with their kind (the fixture's conversion of UniProt's descriptions): 1-based positions.
  const ptms = b.data.annotation({
    ref: "ptms",
    geometry: "ranges",
    target: p53,
    label: "PTM sites",
    source: { host_file: PTMS_FILE, format: "tsv" },
    numbering: "ordinal",
    fields: { beg: "pos" },
  } as never);
  for (const s of P53_STRUCTURES) {
    const af = s.source.service === "afdb";
    b.data.complex({
      ref: s.complex,
      provider: "molstar:structure",
      source: { service: s.source.service, resource: "structure", id: s.source.id },
    });
    const to = `${s.complex}/entity:${s.entity}`;
    b.data.mapping(
      (af
        ? { ref: `${s.complex}.sifts`, from: p53, to, method: "identity" }
        : { ref: `${s.complex}.sifts`, from: p53, to, method: "by_numbering", numbering: "uniprot:P04637" }) as never,
    );
  }
  const coverageOptions = p53CoverageOptions;
  // Each chain's rows on its own unit, placed on P04637 through SIFTS (gaps where it has no residues).
  const coverage = b.data.computed({
    ref: "coverage",
    provider: "molstar:sifts_coverage",
    options: { ...coverageOptions, instance: true },
  } as never);
  // The summary's data: how many structures cover each position, in runs, with a tooltip naming them.
  const depth = b.data.computed({
    ref: "coverage-depth",
    provider: "molstar:sifts_depth",
    options: coverageOptions,
  } as never);
  const plddt = b.data.computed({ ref: "plddt", provider: "molstar:plddt", options: { of: "af-p04637-f1" } } as never);

  const view = b.view({ ref: "main", label: "p53" });
  const sec = view.block({ sequence: p53 });
  // No explicit heights: the section's row sizing gives every row-sized track one row (a lanes track grows by its lanes).
  const track = (ref: string, label: string) => sec.track({ ref, label } as never);
  sec.track({ label: "" }).ruler();
  sec.track({ label: "Sequence" }).letters();
  track(P53_TRACKS.letters, "Residues, hydrophobicity").letters({
    color: { theme: "letter", params: { colors: "Hydrophobicity" } },
  } as never);
  track(P53_TRACKS.domains, "Domains").ranges({
    data: feats,
    filter: { field: "category", values: ["DNA binding", "Motif"] },
    // By type: UniProt's DNA-binding feature has no description, so a palette over labels would leave it uncolored.
    color: byField("category"),
    labels: true,
  } as never);
  // Split by region: collapsed, the regions in their lanes; open, a row per region (its label).
  track(P53_TRACKS.regions, "Regions").ranges({
    data: feats,
    filter: { field: "category", values: ["Region", "Compositional bias"] },
    split: "label",
    color: byField("label"),
  } as never);
  track(P53_TRACKS.ptms, "PTMs").points({
    data: ptms,
    color: { theme: "field_palette", params: { field: "kind", palette: { kind: "categorical", colors: PTM_COLORS } } },
  } as never);
  track(P53_TRACKS.sites, "Binding sites").points({
    data: feats,
    filter: { field: "category", values: ["Binding site", "Site"] },
    color: byField("category"),
  } as never);
  track(P53_TRACKS.variants, "Variants by consequence").points({
    data: variants,
    color: byField("label"),
  } as never);
  track(P53_TRACKS.hydropathy, "Hydropathy").values({
    data: hydro,
    color: {
      theme: "field_palette",
      params: { field: "value", palette: { kind: "continuous", colors: "RdBu", value_domain: [-3, 3], reverse: true } },
    },
  } as never);
  // The PDB coverage in two levels: collapsed, its depth (how many structures cover each position); open, a row per
  // structure, each the union of its chains until it is opened too.
  const g = sec.group({ ref: "g-pdb", label: "PDB coverage", collapsed: true } as never);
  g.track({ ref: P53_TRACKS.coverage, label: "All structures", summary: true } as never).ranges({
    data: depth,
    display: "union",
    color: {
      theme: "field_palette",
      params: {
        field: "value",
        palette: { kind: "continuous", colors: "Blues", value_domain: [0, P53_STRUCTURES.length] },
      },
    },
  } as never);
  for (const [k, s] of P53_STRUCTURES.entries()) {
    g.track({ ref: coverageTrack(s), label: s.label } as never).ranges({
      data: coverage,
      filter: { field: "entry", values: [s.label] },
      display: "union",
      split: "chain",
      color: { theme: "uniform", params: { color: COVERAGE_COLORS[k] } },
    } as never);
  }
  track(P53_TRACKS.plddt, "pLDDT (AlphaFold)").values({ data: plddt } as never);
  return b.getDocument();
}

/** The coverage's and its depth's options (the coverage's rows on the chains' units: `instance`). */
export const p53CoverageOptions = {
  target: "p53",
  accession: "P04637",
  structures: P53_STRUCTURES.map((s) => ({
    complex: s.complex,
    source: s.source,
    label: s.label,
    assembly: "assembly" in s ? s.assembly : null,
  })),
  defaults: { representation: "cartoon", base: "#c8c8c8" },
};
/** Track `s`'s coverage entry's ref. */
export const coverageTrack = (s: P53Structure) => `${P53_TRACKS.coverage}-${s.complex}`;
/** The coverage entries' labels with the one shown in 3D (null: none) marked "· in 3D". */
export function p53ShownPatch(shown: P53Structure | null): Patch {
  const label = (s: P53Structure) => (s === shown ? `${s.label} · in 3D` : s.label);
  const ops = P53_STRUCTURES.map((s) => ({ op: "set" as const, ref: coverageTrack(s), params: { label: label(s) } }));
  return { format: "seq_view_spec_patch", ops };
}

/** Loaders over the fixtures: the structures' bytes by URL, P04637's sequence and its trimmed UniProt files. */
export function p53Loaders(files: {
  readonly structures: ReadonlyMap<string, Uint8Array>;
  readonly fasta: string;
  readonly features: string;
  readonly variants: string;
  /** `P04637-ptms.tsv`: the PTM sites with their kind. */
  readonly ptms: string;
}): Loaders {
  const text: Record<string, string> = {
    [`${UNIPROT}/P04637.fasta`]: files.fasta,
    [`${UNIPROT}/P04637.json`]: files.features,
    [VARIATION]: files.variants,
    [PTMS_FILE]: files.ptms,
  };
  return {
    async load(request) {
      const name = request.uri ?? request.hostFile ?? "";
      const bytes = files.structures.get(name);
      if (bytes) return bytes;
      const t = text[name];
      if (t === undefined) throw new Error(`no fixture for ${name}`);
      return new TextEncoder().encode(t);
    },
  };
}
