/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the V2 tests (no network), written here beside the shared ones in ../fixtures:
// - UniProt responses in the shape of rest.uniprot.org's FASTA and JSON, from the viewer spike's fixtures
//   (lib/seqstar/core/_test/fixtures: UniProtKB CC BY 4.0, wwPDB CC0, Pfam CC0). The ABL1 isoforms are synthetic:
//   the 1IEP construct (residues 225-498) after a made-up stretch, and 1b with a made-up longer first exon.
// - `molstar:structure` and `molstar:secondary_structure` fixture entries over hand-written outputs in the shape
//   seqspec:libraries#the-mol-adapter describes (1BRS, 1IEP with its four Vina poses, 1A3N, an AlphaFold-shaped
//   barnase-barstar model). The structure "files" are JSON of those outputs; the provider applies `chains` and `model`.
// - Small hand-written files for the other cases (antibody domains, an SH2 excerpt, a tRNA-like family, a DNA
//   construct, a kinase MSA, a designed dimer); values such as disorder, effects and pLDDT are synthetic.
// - Cases 19 and 20's streamed genome: stream-provider.ts's synthetic GRCh38 (letters, genes with transcript T1,
//   variants, a phyloP-like signal), served by `fixtureLoaders` under the `hg38_*` stream ids; case 22's family, the
//   committed PF00042 seed (../fixtures/pf00042.json).

import { readFileSync } from "node:fs";
import type { MoleculeInput, UnitInput } from "@molstar/seq-core/model";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { proceduralRow, tileServer } from "@molstar/seq-core/stream/tiles/server";
import {
  createRegistries,
  type FormatEntry,
  int,
  KernelRegistry,
  type Loaders,
  list,
  literal,
  MvsRegistry,
  nullable,
  OptionalField,
  RequiredField,
  ref,
  register,
  SimpleParamsSchema,
  type SvsRegistries,
  type SvsRegistryEntry,
  str,
  union,
} from "../../src/format/index.ts";
import { pf00042, SH2_STOCKHOLM } from "../model/alignment-fixtures.ts";
import { HEM_ATOMS } from "../model/fixtures.ts";
import {
  CODING_CDS,
  codingProtein,
  GENE_LEVELS,
  GRCH38,
  hg38Streams,
  itemsManifest,
  type Scripted,
  scripted,
  TRANSCRIPT,
  transcriptItems,
  VOCABULARY,
} from "./stream-provider.ts";

const spike = (name: string) => JSON.parse(readFileSync(new URL(`../fixtures/${name}.json`, import.meta.url), "utf8"));
export const P53 = spike("p53");
export const BRS = spike("1brs");
export const IEP = spike("1iep");

// --- UniProt -------------------------------------------------------------------------------------------------

const filler = (n: number, seed: string) => Array.from({ length: n }, (_, i) => seed[(i * 7) % seed.length]).join("");
/** ABL1 1a: 224 made-up residues, then the 1IEP construct, so 1IEP's auth numbers are 1a's. */
export const ABL1A = filler(224, "MLEICLKLVGCKSKKGLSSSSSCYLEEALQRPVASDFEPQGLSEAARWNSKENLLAG") + IEP.sequence;
/** ABL1 1b: a first exon 19 residues longer; from 1a index 26 on, equal. */
export const ABL1B = filler(45, "MGQQPGKVLGDQRRPSLPALHFIKGAGKKESSRHGGPHCNVFVEHE") + ABL1A.slice(26);
export const ISOFORM_SEGMENTS: [number, number, null][] = [[26, 45, null]];

const fasta = (header: string, seq: string) => `>${header}\n${seq.replace(/(.{60})/g, "$1\n")}\n`;
const features = (list: [string, number, number, string][]) =>
  JSON.stringify({
    primaryAccession: "x",
    features: list.map(([type, s, e, description]) => ({
      type,
      location: { start: { value: s, modifier: "EXACT" }, end: { value: e, modifier: "EXACT" } },
      description,
    })),
  });
const UNIPROT = "https://rest.uniprot.org/uniprotkb";

/**
 * Case 13's p53 features (synthetic, in the shape of rest.uniprot.org's JSON; no network): types and positions after
 * UniProt's P04637 entry, descriptions abridged, a few made up (the sequence conflicts).
 */
const P53_FEATURES: [string, number, number, string][] = [
  ["Chain", 1, 393, "Cellular tumor antigen p53"],
  ["Modified residue", 15, 15, "Phosphoserine; by ATM"],
  ["Modified residue", 20, 20, "Phosphoserine; by CHEK2"],
  ["Modified residue", 33, 33, "Phosphoserine"],
  ["Modified residue", 46, 46, "Phosphoserine; by HIPK2"],
  ["Modified residue", 392, 392, "Phosphoserine"],
  ...[370, 372, 373, 381, 382, 386].map((p): [string, number, number, string] => ["Cross-link", p, p, "Ubiquitin"]),
  ["Compositional bias", 64, 92, "Pro residues"],
  ["Compositional bias", 359, 384, "Basic and acidic residues"],
  ["Sequence conflict", 155, 155, "T → I"],
  ["Sequence conflict", 337, 337, "R → C"],
  ["Domain", 94, 292, "DNA-binding"],
  ["Domain", 325, 356, "Tetramerization"],
  ["Region", 1, 40, "Transactivation 1"],
  ["Region", 41, 61, "Transactivation 2"],
  ["Region", 1, 68, "Disordered"],
  ["Region", 359, 393, "Basic (repression of DNA binding)"],
  ["Motif", 305, 322, "Bipartite nuclear localization signal"],
  ["Motif", 339, 350, "Nuclear export signal"],
  ...[176, 179, 238, 242].map((p): [string, number, number, string] => ["Binding site", p, p, "Zn(2+)"]),
  ...[175, 245, 248, 249, 273, 282].map((p): [string, number, number, string] => ["Site", p, p, "Mutational hotspot"]),
];
/** Case 13's p53 variants (synthetic, in the shape of the Proteins API's variation; no network). */
const P53_VARIANTS: [number, string, string, string][] = [
  [47, "P", "S", "polymorphism"],
  [72, "P", "R", "polymorphism"],
  ...["R175H", "G245S", "R248Q", "R248W", "R249S", "R273H", "R273C", "R282W"].map(
    (v): [number, string, string, string] => [Number(v.slice(1, -1)), v[0], v.slice(-1), "somatic hotspot"],
  ),
];
const variation = (sequence: string, list: [number, string, string, string][]) =>
  JSON.stringify({
    accession: "P04637",
    sequence,
    features: list.map(([p, wildType, alternativeSequence, kind]) => ({
      type: "VARIANT",
      begin: String(p),
      end: String(p),
      wildType,
      alternativeSequence,
      consequenceType: "missense",
      description: `${wildType}${p}${alternativeSequence} (${kind})`,
    })),
  });

// --- Structures, as `molstar:structure` yields them -------------------------------------------------------------

interface Recorded {
  molecules: MoleculeInput[];
  units: UnitInput[];
  categories?: Record<string, Record<string, unknown[]>>;
  secondary?: { chain: string; start: number; stop: number; kind: string }[];
}
const mvs = (label: string, auth = label) => ({
  mvs: { label_asym_id: label, auth_asym_id: auth, instance_id: "ASM-1" },
});
const authUnit = (
  key: string,
  molecule: string,
  first: number,
  length: number,
  observed?: [number, number][],
): UnitInput => ({
  key,
  molecule,
  observed,
  origin: "experimental",
  numberings: [{ name: "auth", method: "segments", segments: [[0, first, length]], isDefault: true }],
  custom: mvs(key),
});
const label = { name: "label", method: "offset" as const, first: 1 };
export const BARNASE_MATURE = BRS.barnase.sequence.slice(47);
export const BARSTAR_1BRS = "KKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDALTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS";

function brs1(): Recorded {
  const secondary = (chain: string, feats: [string, number, number][], offset: number) =>
    feats
      .filter(([t]) => t === "Helix" || t === "Beta strand")
      .map(([t, s, e]) => ({
        chain,
        start: s - offset - 1,
        stop: e - offset,
        kind: t === "Helix" ? "helix" : "strand",
      }))
      .filter((r) => r.start >= 0);
  return {
    molecules: [
      { key: "entity:1", residues: BARNASE_MATURE, label: "Barnase", numberings: [label] },
      { key: "entity:2", residues: BARSTAR_1BRS, label: "Barstar", numberings: [label] },
    ],
    units: [
      ...["A", "B", "C"].map((c) => authUnit(c, "entity:1", 1, 110, [[2, 110]])),
      ...["D", "E", "F"].map((c) =>
        authUnit(c, "entity:2", 1, 89, [
          [0, 63],
          [65, 89],
        ]),
      ),
    ],
    secondary: ["A", "B", "C"]
      .flatMap((c) => secondary(c, BRS.barnase.features, 47))
      .concat(["D", "E", "F"].flatMap((c) => secondary(c, BRS.barstar.features, 1))),
  };
}

const elements = (atoms: string[]) => atoms.map((a) => (a === "FE" ? "Fe" : a[0]));
/** STI's atom names: element and 1-based index (`N13` is atom 12). */
export const STI_NAMES = (IEP.ligand.atoms as string[]).map((e, i) => `${e}${i + 1}`);
const stiMolecule = (key: string): MoleculeInput => ({
  key,
  kind: "ligand",
  elements: IEP.ligand.atoms,
  chemCompIds: ["STI"],
  label: "imatinib",
});
const atomUnit = (key: string, molecule: string, extra: Partial<UnitInput> = {}): UnitInput => ({
  key,
  molecule,
  origin: "experimental",
  numberings: [{ name: "atom", method: "names", names: STI_NAMES, isDefault: true }],
  ...extra,
});

function iep1(): Recorded {
  return {
    molecules: [
      {
        key: "entity:1",
        residues: IEP.sequence,
        label: "Proto-oncogene tyrosine-protein kinase ABL1",
        numberings: [label],
      },
      stiMolecule("entity:2"),
    ],
    units: [authUnit("A", "entity:1", 225, IEP.sequence.length), atomUnit("C", "entity:2", { custom: mvs("C", "A") })],
  };
}

function poses(): Recorded {
  return {
    molecules: [stiMolecule("entity:ligand")],
    units: [1, 2, 3, 4].map((k, i) =>
      atomUnit(`${i}/ligand`, "entity:ligand", {
        origin: "docked",
        score: IEP.poses[k].score,
        label: IEP.poses[k].label,
      }),
    ),
  };
}

const HBA =
  "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
const HBB =
  "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
function a3n1(): Recorded {
  return {
    molecules: [
      { key: "entity:1", residues: HBA, label: "Hemoglobin subunit alpha", numberings: [label] },
      { key: "entity:2", residues: HBB, label: "Hemoglobin subunit beta", numberings: [label] },
      { key: "entity:3", kind: "ligand", elements: elements(HEM_ATOMS), chemCompIds: ["HEM"], label: "HEM" },
    ],
    units: [
      authUnit("A", "entity:1", 1, 141),
      authUnit("C", "entity:1", 1, 141),
      authUnit("B", "entity:2", 1, 146),
      authUnit("D", "entity:2", 1, 146),
      ...["E", "F", "G", "H"].map((c) => ({
        key: c,
        molecule: "entity:3",
        origin: "experimental",
        numberings: [{ name: "atom", method: "names" as const, names: HEM_ATOMS, isDefault: true }],
      })),
    ],
    secondary: [
      { chain: "A", start: 2, stop: 18, kind: "helix" },
      { chain: "B", start: 4, stop: 17, kind: "helix" },
    ],
  };
}

/**
 * A synthetic pLDDT for an AlphaFold-shaped chain (deterministic and smooth; a stand-in for a model's own scores): about
 * 97 inside the helices and strands (`folded`, 0-based start and stop), 85 in the loops between, tapering toward low at
 * the termini, and low (40 to 68) across a signal peptide or propeptide (the first `floppy` residues). A triangular
 * moving average joins the steps; a faint ripple keeps neighbors from being equal.
 */
function afPlddt(length: number, folded: [number, number][], floppy = 0): number[] {
  const raw = Array.from({ length }, (_, i) => {
    const core = folded.some(([a, b]) => i >= a && i < b);
    const base = i < floppy ? 40 + 0.6 * i : core ? 97 : 85;
    return Math.min(base, 62 + 7 * Math.min(i, length - 1 - i));
  });
  return raw.map((_, i) => {
    let sum = 0;
    let weight = 0;
    for (let k = -2; k <= 2; k++) {
      const j = i + k;
      if (j < 0 || j >= length) continue;
      sum += raw[j] * (3 - Math.abs(k));
      weight += 3 - Math.abs(k);
    }
    return Math.round((sum / weight + 1.2 * Math.sin(i * 0.9)) * 100) / 100;
  });
}

/** An AlphaFold-shaped model of full-length barnase and barstar, with a synthetic pLDDT per residue. */
function afModel(): Recorded {
  const chains = [BRS.barnase, BRS.barstar] as { sequence: string; features: [string, number, number][] }[];
  // The model's secondary structure: UniProt's helices and strands, 0-based. Barnase's signal peptide and propeptide
  // (residues 1-47) are a floppy tail with none.
  const floppy = [47, 0];
  const secondary = chains.flatMap((c, k) =>
    c.features
      .filter(([t, s]) => (t === "Helix" || t === "Beta strand") && s > floppy[k])
      .map(([t, s, e]) => ({ chain: "AB"[k], start: s - 1, stop: e, kind: t === "Helix" ? "helix" : "strand" })),
  );
  const plddt: Record<string, unknown[]> = { label_asym_id: [], label_seq_id: [], metric_value: [] };
  chains.forEach((c, k) => {
    const folded = secondary.filter((r) => r.chain === "AB"[k]).map((r): [number, number] => [r.start, r.stop]);
    afPlddt(c.sequence.length, folded, floppy[k]).forEach((v, i) => {
      plddt.label_asym_id.push("AB"[k]);
      plddt.label_seq_id.push(i + 1);
      plddt.metric_value.push(v);
    });
  });
  const seqs = chains.map((c) => c.sequence);
  return {
    molecules: seqs.map((residues, k) => ({ key: `entity:${k + 1}`, residues, numberings: [label] })),
    units: seqs.map((s, k) => ({ ...authUnit("AB"[k], `entity:${k + 1}`, 1, s.length), origin: "predicted" })),
    categories: { ma_qa_metric_local: plddt },
    secondary,
  };
}

const json = (v: unknown) => JSON.stringify(v);

const structureFormat = (name: string): FormatEntry => ({
  name,
  gives: ["object", "table"],
  reads: ["want", "categoryName"],
  read(bytes, request, ctx) {
    const value = JSON.parse(new TextDecoder().decode(bytes)) as Recorded;
    if (request.categoryName == null) return { kind: "object", value };
    const columns = value.categories?.[request.categoryName];
    if (!columns) {
      ctx.report({ code: "svs.source.category", severity: "error", path: "", message: `no ${request.categoryName}` });
      return { kind: "table", columns: {}, length: 0 };
    }
    return { kind: "table", columns, length: Object.values(columns)[0].length };
  },
});

/** AlphaFold DB's PAE JSON as `@molstar/seq-molstar-poc`'s `afdb_pae` reads it (its current form): an N × N matrix. */
const paeFormat: FormatEntry = {
  name: "afdb_pae",
  gives: ["object"],
  reads: [],
  read(bytes) {
    const rows = JSON.parse(new TextDecoder().decode(bytes))[0].predicted_aligned_error as number[][];
    return { kind: "object", value: { kind: "afdb_pae", size: rows.length, values: Float64Array.from(rows.flat()) } };
  },
};

/** `@molstar/seq-molstar-poc`'s entries as the cases use them, over recorded outputs. */
export const MolstarFixture: SvsRegistryEntry = {
  formats: [...["mmcif", "sdf"].map(structureFormat), paeFormat],
  services: [
    {
      name: "pdbe",
      resources: ["structure", "sifts"],
      resolve: (id) => ({ uri: `https://www.ebi.ac.uk/pdbe/entry-files/download/${id}_updated.cif`, format: "mmcif" }),
    },
    {
      name: "afdb",
      resources: ["structure", "pae"],
      resolve: (id, resource) =>
        resource === "pae"
          ? { uri: `https://alphafold.ebi.ac.uk/files/${id}-predicted_aligned_error_v4.json`, format: "afdb_pae" }
          : { uri: `https://alphafold.ebi.ac.uk/files/${id}-model_v4.cif`, format: "mmcif" },
    },
  ],
  complexProviders: [
    {
      name: "molstar:structure",
      options: SimpleParamsSchema({
        model: OptionalField(union(int, literal("all")), 0, "the model, 0-based; all: units per model"),
        assembly: OptionalField(nullable(str), null, "an assembly id"),
        include: OptionalField(
          list(literal("polymer", "ligand", "ion", "water", "branched")),
          ["polymer", "ligand", "ion", "branched"],
          "",
        ),
        chains: OptionalField(nullable(list(str)), null, "keep only the units with these names"),
        unit_name: OptionalField(str, "label_asym_id", ""),
        numberings: OptionalField(list(literal("label", "auth", "atom")), ["label", "auth", "atom"], ""),
        default_numbering: OptionalField(nullable(str), "auth", ""),
      }),
      formats: ["mmcif", "bcif", "pdb", "sdf", "mol2"],
      async build(_ctx, source, options) {
        const rec = (source as { value: Recorded }).value;
        const o = options as { chains: string[] | null };
        const units = rec.units.filter((u) => !o.chains || o.chains.includes(u.key));
        const used = new Set(units.map((u) => u.molecule));
        return { molecules: rec.molecules.filter((m) => used.has(m.key)), units };
      },
    },
  ],
  providers: [
    {
      // `molstar:pae`: the matrix over the complex's units in order, a grid row per pair (side a scored, b aligned on).
      name: "molstar:pae",
      options: SimpleParamsSchema({
        of: RequiredField(ref("complex"), "the complex"),
        source: RequiredField(ref("source"), "its PAE"),
      }),
      geometry: "grid",
      async compute(ctx, options) {
        const o = options as { of: string; source: string };
        const m = (ctx.source(o.source) as { value: { size: number; values: Float64Array } } | null)?.value;
        const units = ["A", "B"].map((c) => `${o.of}/${c}`);
        const n = units.map((u) => ctx.set?.units.key.findIndex((k) => ctx.set?.strings.get(k) === u) ?? -1);
        const len = n.map((u) => (ctx.set && u >= 0 ? ctx.set.units.rangeStop[u] - ctx.set.units.rangeStart[u] : 0));
        if (!m || len[0] + len[1] !== m.size) return null;
        const columns: Record<string, unknown[]> = { target: [], b_target: [], value: [] };
        for (let p = 0; p < 2; p++)
          for (let q = 0; q < 2; q++) {
            const at = (i: number, j: number) => m.values[((q ? len[0] : 0) + j) * m.size + (p ? len[0] : 0) + i];
            columns.target.push(units[p]);
            columns.b_target.push(units[q]);
            columns.value.push(
              Float32Array.from({ length: len[p] * len[q] }, (_, x) => at(Math.floor(x / len[q]), x % len[q])),
            );
          }
        return { key: "", geometry: "grid", semanticType: "pae", label: "PAE", columns };
      },
    },
    {
      name: "molstar:secondary_structure",
      options: SimpleParamsSchema({ of: RequiredField(ref("complex"), "the complex") }),
      geometry: "ranges",
      async compute(ctx, options) {
        const of = (options as { of: string }).of;
        const rec = (ctx.source(of) as { value: Recorded } | null)?.value;
        const rows = rec?.secondary ?? [];
        const columns = {
          target: rows.map((r) => `${of}/${r.chain}`),
          start: rows.map((r) => r.start),
          stop: rows.map((r) => r.stop),
          category: rows.map((r) => r.kind),
        };
        return { key: "", geometry: "ranges", columns };
      },
    },
  ],
};

export function registries(molstar = true): SvsRegistries {
  const reg = createRegistries();
  register(reg, molstar ? [KernelRegistry, MvsRegistry, MolstarFixture] : [KernelRegistry, MvsRegistry]);
  return reg;
}

// --- Files of the other cases ------------------------------------------------------------------------------------

export const VH =
  "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
export const VL =
  "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";
const range = (a: number, b: number) => Array.from({ length: b - a + 1 }, (_, i) => a + i);
const sfx = (n: number, s: string) => [...s].map((c) => `${n}${c}`);
/** Kabat on VH (no 35A or 35B): 52A, 82A–C, 100A–C. */
export const VH_KABAT = [
  ...range(1, 52),
  "52A",
  ...range(53, 82),
  ...sfx(82, "ABC"),
  ...range(83, 100),
  ...sfx(100, "ABC"),
  ...range(101, 113),
];
export const VL_KABAT = range(1, 107);
const imgt = (missing: number[]) =>
  range(1, 128)
    .filter((n) => !missing.includes(n))
    .flatMap((n) => (n === 112 ? ["111.1", "112.1", 112] : [n]));
export const VH_IMGT = imgt([10, 31, 32, 33, 34, 60, 61, 62, 73, 84]);
export const VL_IMGT = imgt([10, ...range(29, 37), ...range(57, 64), 73, 81, 82]).filter(
  (x) => x !== "111.1" && x !== "112.1",
);

/** A tRNA-like family: SS_cons in WUSS, one row with a descending header. */
export const TRNA_STO = `# STOCKHOLM 1.0
#=GF ID   tRNA-like
#=GS T1/1-17  OS Homo sapiens
T1/1-17       GCGGAUUUAGCUCAGUU
T2/1-16       GCGGA-UUAGCUCAGUU
T3/20-5       GCCGAUUUAGCUCGGUU
T4/1-17       GCGGAUUUAGCACAGCU
#=GR T1/1-17 SS <<<<....>>>>.....
#=GC SS_cons   <<<<....>>>>.....
//
`;
/** RNAfold's structure for row 3's 17 residues. */
export const RNAFOLD_3 = "((((....))))(...)";

export const SH2_NWK = "((SRC_HUMAN,GRB2_HUMAN),(PTN11_HUMAN,P85A_HUMAN));";
export const SH2_SITES =
  "row\tposition\tcategory\taccession\n0\t153\tphosphosite\tP12931\n2\t62\tphosphosite\tP62993\n1\t20\tother\tQ06124\n";
export const SH2_ROWS = "row\tspecies\n0\tHomo sapiens\n1\tHomo sapiens\n2\tHomo sapiens\n3\tHomo sapiens\n";

export const INSERT = "ATGGAATTCAAAGGTGAAGAACTGTTCACCGGTGTTGTTCCGATCCTGGTTGAACTGGACGGTCTCGAGTAA";
export const INSERT_TSV =
  "start\tstop\tlabel\tcategory\tstrand\n0\t3\tstart codon\tstart_codon\t+\n3\t9\tEcoRI\trestriction_site\t+\n" +
  "9\t63\tGFP fragment\tCDS\t+\n63\t69\tXhoI\trestriction_site\t+\n";

/**
 * A synthetic 240-nucleotide cyclic plasmid (CY1's story; CG1b draws it as a ring): random letters but a start codon at
 * 40, EcoRI at 60 and HindIII at 84; its replication origin runs across the origin, [200, 30), written start > stop.
 */
export const PLASMID =
  "AAGCGAGACAGTGCCCATTAGGAATGCGCAACAAATAATCATGTAACGGTTCGGGCATGGGAATTCATAAGGACATGTTGGGGGAAGCTTCTGTGCAGACACGCGAGAAAC" +
  "GTTCATGGCTCAACGGGAATCTGGGATTGTCCCTATGAGTAGCGCTCGATTGCACAAAAAGAGTTGGCATTAGCATCTGCTTGTCATGGCGCGCCGGGACTTGGTCTGCAGG" +
  "GGTCGTGTGAACTTGCA";
export const PLASMID_TSV =
  "id\tstart\tstop\tlabel\tcategory\nori\t200\t30\trep origin\trep_origin\nlacZ\t40\t112\tlacZ\u03b1 fragment\tCDS\n" +
  "mcs\t58\t92\tMCS\tmisc_feature\nampR\t120\t192\tAmpR fragment\tCDS\n";

const KD = "VAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNL";
const mutate = (s: string, at: Record<number, string>) => [...s].map((c, i) => at[i] ?? c).join("");
/** Six kinase domain stretches, the gatekeeper at column 48 (ABL1 T315); synthetic around ABL1's. */
export const GK_COLUMN = 48;
export const KINASE_ROWS = [
  KD,
  mutate(KD, { 48: "I" }),
  mutate(KD, { 2: "I", 7: "P", 8: "G", 41: "-", 48: "T", 52: "S" }),
  mutate(KD, { 4: "M", 7: "P", 8: "S", 30: "M", 48: "T", 53: "C" }),
  mutate(KD, { 4: "M", 7: "P", 9: "A", 30: "I", 48: "T", 54: "F" }),
  mutate(KD, { 2: "I", 5: "E", 8: "A", 33: "R", 48: "M", 50: "P" }),
];

export const DESIGN_7 = "MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEV";
export const HOMOLOGS_STO = `# STOCKHOLM 1.0
#=GS design/1-55 DE Design 7
design/1-55   MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEV
hom1/3-57     MSEEKLRELAEKAIRLAKEGNPEEALRLLEEALRLAKESGDPELIARVLELLAEV
hom2/1-54     MSEEKLKELAEKAIRLAKEG-PEEALKLLEEALKLAKESGDPELIARVLELLAEV
#=GC SS_cons  ..HHHHHHHHHHHHHHHHH....HHHHHHHHHHHHHHH....HHHHHHHHHHHHH
//
`;

/**
 * Case 1's disorder profile (synthetic, deterministic; shaped like IUPred for p53): disordered in the N-terminal
 * transactivation and proline-rich parts (about 1-95), ordered across the DNA-binding domain (about 95-290),
 * disordered again in the C-terminus (about 290-393) with a dip in the tetramerization domain (about 325-356).
 * Smooth steps (logistic) plus a gentle ripple; a stand-in for real predictor output.
 */
function p53Disorder(p: number): number {
  const step = (x: number) => 1 / (1 + Math.exp(-x));
  const tetramer = step((p - 322) / 4) - step((p - 359) / 4);
  const d = 0.1 + 0.7 * step((95 - p) / 6) + 0.65 * step((p - 290) / 6) - 0.45 * tetramer + 0.05 * Math.sin(p / 7);
  return Math.min(1, Math.max(0, d));
}

/** Every file the cases read, by URI or host file name: a copy of the map, its texts made once. */
let FILES: Record<string, string> | undefined;
export function files(): Record<string, string> {
  FILES ??= makeFiles();
  return { ...FILES };
}
function makeFiles(): Record<string, string> {
  const barnaseFeats = BRS.barnase.features;
  const barstarFeats = BRS.barstar.features;
  const iupred = ["position,iupred", ...range(1, 393).map((p) => `${p},${p53Disorder(p).toFixed(2)}`)].join("\n");
  // Variant effects (synthetic): a score per state on every residue of mature barnase (precursor 47–156, the residues
  // 1BRS's chain A shows), so the heatmap lands on the drawn chain.
  const mature = range(47, 156);
  const effects = { index: mature, value: mature.map((i) => range(0, 19).map((s) => ((i + s) % 7) / 6 - 0.5)) };
  return {
    [`${UNIPROT}/P04637.fasta`]: fasta("sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens", P53.sequence),
    [`${UNIPROT}/P00648.fasta`]: fasta("sp|P00648|RNBR_BACAM Ribonuclease", BRS.barnase.sequence),
    [`${UNIPROT}/P11540.fasta`]: fasta("sp|P11540|BARS_BACAM Barstar", BRS.barstar.sequence),
    [`${UNIPROT}/P00519-1.fasta`]: fasta("sp|P00519-1|ABL1_HUMAN Isoform IA (synthetic)", ABL1A),
    [`${UNIPROT}/P00519-2.fasta`]: fasta("sp|P00519-2|ABL1_HUMAN Isoform IB (synthetic)", ABL1B),
    [`${UNIPROT}/P00648.json`]: features(barnaseFeats),
    [`${UNIPROT}/P11540.json`]: features(barstarFeats),
    [`${UNIPROT}/P04637.json`]: features(P53_FEATURES),
    "https://www.ebi.ac.uk/proteins/api/variation/P04637": variation(P53.sequence, P53_VARIANTS),
    "https://example.org/p53-iupred.csv": iupred,
    "https://www.ebi.ac.uk/pdbe/entry-files/download/1brs_updated.cif": json(brs1()),
    "https://www.ebi.ac.uk/pdbe/entry-files/download/1iep_updated.cif": json(iep1()),
    "https://www.ebi.ac.uk/pdbe/entry-files/download/1a3n_updated.cif": json(a3n1()),
    "https://alphafold.ebi.ac.uk/files/AF-BARNASE-BARSTAR-model_v4.cif": json(afModel()),
    "https://alphafold.ebi.ac.uk/files/AF-BARNASE-BARSTAR-predicted_aligned_error_v4.json": readFileSync(
      new URL("../fixtures/af-barnase-barstar-pae.json", import.meta.url),
      "utf8",
    ),
    "https://example.org/pfam/PF00042/full.nwk": pf00042Family().newick,
    "https://example.org/pfam/PF00042/rows.json": pf00042Family().rowsJson,
    "1iep_ligand_vina_poses.sdf": json(poses()),
    "trastuzumab.fasta": fasta("VH trastuzumab heavy chain variable domain", VH) + fasta("VL light chain", VL),
    "PF00017.seed.sto": SH2_STOCKHOLM,
    "PF00017.nwk": SH2_NWK,
    "sh2-sites.tsv": SH2_SITES,
    "sh2-rows.tsv": SH2_ROWS,
    "RF00005.seed.sto": TRNA_STO,
    "pET-insert.fasta": fasta("pET-insert", INSERT),
    "pET-insert.tsv": INSERT_TSV,
    "pDemo.fasta": fasta("pDemo (synthetic)", PLASMID),
    "pDemo.tsv": PLASMID_TSV,
    "barnase-effects.json": json(effects),
  };
}

/** Vertebrates with globins, the species the synthetic PF00042 rows run through. */
const GLOBIN_SPECIES = [
  "Homo sapiens",
  "Pan troglodytes",
  "Gorilla gorilla",
  "Macaca mulatta",
  "Mus musculus",
  "Rattus norvegicus",
  "Oryctolagus cuniculus",
  "Bos taurus",
  "Sus scrofa",
  "Equus caballus",
  "Canis lupus",
  "Felis catus",
  "Loxodonta africana",
  "Monodelphis domestica",
  "Ornithorhynchus anatinus",
  "Gallus gallus",
  "Taeniopygia guttata",
  "Anas platyrhynchos",
  "Anolis carolinensis",
  "Python bivittatus",
  "Xenopus laevis",
  "Xenopus tropicalis",
  "Danio rerio",
  "Takifugu rubripes",
  "Oryzias latipes",
  "Salmo salar",
  "Latimeria chalumnae",
  "Callorhinchus milii",
];

/**
 * Case 9's rows.json and full.nwk (synthetic, deterministic; built once). Rows run through species in runs of 6-45 that
 * share one (HBA's and HBB's runs are Homo sapiens); identity falls from 1 at HBA (0.95 at HBB) to about 0.3 far from
 * either, with a little wobble. The tree nests the same leaves in the same order: each run a balanced clade, the clades
 * balanced above them, an ultrametric tree: the leaves end together, every level of clade showing across the depth.
 */
let family: { rowsJson: string; newick: string } | undefined;
function pf00042Family() {
  if (family) return family;
  const { rows, hbaRow, hbbRow } = PF00042;
  const next = (x: number) => (Math.imul(x, 1664525) + 1013904223) >>> 0;
  const runs: { from: number; to: number; species: string }[] = [];
  for (let from = 0, state = 7, prev = -1; from < rows; ) {
    state = next(state);
    const to = Math.min(rows, from + 6 + ((state >>> 8) % 40));
    state = next(state);
    let sp = (state >>> 8) % GLOBIN_SPECIES.length;
    if (sp === prev) sp = (sp + 1) % GLOBIN_SPECIES.length;
    if ((from <= hbaRow && hbaRow < to) || (from <= hbbRow && hbbRow < to)) sp = 0;
    runs.push({ from, to, species: GLOBIN_SPECIES[sp] });
    prev = sp;
    from = to;
  }
  const meta: { row: number; species: string; identity: number }[] = [];
  for (const run of runs)
    for (let k = run.from; k < run.to; k++) {
      const near = Math.max(Math.exp(-Math.abs(k - hbaRow) / 2500), 0.95 * Math.exp(-Math.abs(k - hbbRow) / 2500));
      const wobble = ((Math.imul(k, 2654435761) >>> 0) % 1000) / 1000;
      const identity = k === hbaRow ? 1 : k === hbbRow ? 0.95 : 0.3 + 0.55 * near + 0.15 * wobble * (1 - 0.5 * near);
      meta.push({ row: k, species: run.species, identity: Math.round(identity * 100) / 100 });
    }
  // A balanced binary tree over items[lo, hi), ultrametric: a node sits above its taller child by a step (0.05-0.2, by a
  // hash; 0.2-0.6 over two leaves), a child's edge reaching up to it, so the leaves end together and every level shows.
  type Item = { text: string; h: number };
  const balanced = (items: Item[], lo: number, hi: number): Item => {
    if (hi - lo === 1) return items[lo];
    const mid = (lo + hi) >>> 1;
    const l = balanced(items, lo, mid);
    const r = balanced(items, mid, hi);
    const u = ((Math.imul(lo * 31 + hi, 2654435761) >>> 0) % 1000) / 1000;
    const h = Math.max(l.h, r.h) + (l.h === 0 && r.h === 0 ? 0.2 + 0.4 * u : 0.05 + 0.15 * u);
    return { text: `(${l.text}:${(h - l.h).toFixed(3)},${r.text}:${(h - r.h).toFixed(3)})`, h };
  };
  const leaves = range(0, rows - 1).map((k): Item => ({ text: `r${k}`, h: 0 }));
  const clades = runs.map((r) => balanced(leaves.slice(r.from, r.to), 0, r.to - r.from));
  family = { rowsJson: json(meta), newick: `${balanced(clades, 0, clades.length).text};` };
  return family;
}

/**
 * Case 9's PF00042 "full" alignment (synthetic): procedural rows, HBA_HUMAN at `hbaRow` and HBB_HUMAN at `hbbRow`
 * (1A3N's α and β, aligned over the same columns), conservation served.
 */
export const PF00042 = {
  uri: "https://example.org/pfam/PF00042/full",
  rows: 18000,
  columns: 320,
  hbaRow: 4711,
  hbbRow: 9133,
};
export const pf00042Server = () =>
  tileServer({
    rowCount: PF00042.rows,
    columnCount: PF00042.columns,
    row: (k) =>
      k === PF00042.hbaRow
        ? HBA.padStart(HBA.length + 40, "-").padEnd(PF00042.columns, "-")
        : k === PF00042.hbbRow
          ? HBB.padStart(HBB.length + 38, "-").padEnd(PF00042.columns, "-")
          : proceduralRow(k, 320),
    stats: { conservation: (c) => (c % 10) / 10 },
  });

/** Case 20's coding gene on chr2 of the streamed genome (stream-provider.ts): its CDS, its protein and transcript. */
export { CODING_CDS, TRANSCRIPT };
export const P1 = codingProtein();
/** Case 22's family: Pfam PF00042.29's 32 globins as ../fixtures/pf00042.json commits them (CC0), HBA_HUMAN first. */
export const GLOBINS = pf00042();

/**
 * The streamed genome of cases 19 and 20, by stream id: GRCh38's 24 chromosomes (synthetic letters, by region), its
 * genes (with transcript T1 on chr2, whose rows carry its codon `map`), variants and a phyloP-like signal.
 */
function genomeStreams(): (source: { stream: string | null }) => Scripted["provider"] | null {
  const s = hg38Streams({ letters: true });
  const keys = GRCH38.map((c) => c[0]);
  const genes = scripted(itemsManifest(GRCH38, GENE_LEVELS, VOCABULARY.genes), (r) => transcriptItems(r, keys));
  const by: Record<string, Scripted> = {
    hg38_genome: s.genome,
    hg38_genes: genes,
    hg38_variants: s.variants,
    hg38_phylop: s.phylop,
  };
  return (source) => (source.stream ? (by[source.stream]?.provider ?? null) : null);
}

/**
 * Fixture loaders: a URI or host file → its bytes, anything else failing as a missing file would; PF00042 on the tile
 * server, its requests in `tiles`; the streamed genome by its `hg38_*` stream ids.
 */
export function fixtureLoaders(
  map: Record<string, string> = files(),
  server = pf00042Server(),
): Loaders & { requests: string[]; tiles: string[] } {
  const requests: string[] = [];
  const tiles = tilesStream((uri) => (uri === PF00042.uri ? server.transport : null));
  const genome = genomeStreams();
  return {
    requests,
    tiles: server.calls,
    stream: (source) => tiles(source) ?? genome(source),
    async load(request) {
      const name = request.uri ?? request.hostFile ?? "";
      requests.push(name);
      const text = map[name];
      if (text === undefined) throw new Error(`no fixture for ${name}`);
      return new TextEncoder().encode(text);
    },
  };
}
