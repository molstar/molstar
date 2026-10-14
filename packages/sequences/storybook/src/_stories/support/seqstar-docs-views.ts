/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The documents of the Seq* stories that build one in code (lib/seqstar/storybook/src/_stories/*.stories.tsx), kept out of the
// story files so the smoke test (lib/seqstar/storybook/_test/smoke.test.ts) can open each on a headless viewer.
// No React here: a story imports its builders from this module, and nothing in it changes what a story shows.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import type { DocNode, Op, SectionStyle, SeqViewSpecDoc } from "@molstar/seq-core/format";
import { formatLink } from "@molstar/seq-core/links";
import { bookmarksOf } from "@molstar/seq-viewer";
import { CASES, CUT_ECORI, CUT_XHOI, case2, case10, contactsAD } from "../../../../core/_test/format/cases.ts";
import { BARNASE_MATURE } from "../../../../core/_test/format/load-fixtures.ts";
import { familyRows, familyWithTree } from "../../../../viewer/_test/families.ts";

// ---- Seq*/Viewer ----

/**
 * A synthetic family evolved down a known tree with unequal rates (deterministic, not real sequences): three clades off
 * a root, one with two subclades and a fast-evolving member (A5), and a distant outgroup. Rows come in a fixed shuffle,
 * so the data order is not the tree's.
 */
function evolvedFamily(): { labels: string[]; rows: string[] } {
  const AMINO = "ACDEFGHIKLMNPQRSTVWY";
  let seed = 11;
  const rnd = () => {
    seed = (Math.imul(seed, 1103515245) + 12345) & 0x7fffffff;
    return seed / 0x7fffffff;
  };
  const mutate = (from: string, rate: number) =>
    [...from]
      .map((c) => (c === "-" ? c : rnd() < rate * 0.1 ? "-" : rnd() < rate ? AMINO[Math.floor(rnd() * 20)] : c))
      .join("");
  const root = "WYFGKITRRESERLLLNAENPRGTFLVRESWFHPNITGVEKDAENLLLTRGVDGSFLARPWYWGDISREEVNEKLRDTADG";
  const leaves: [string, string][] = [];
  const clade = (name: string, rate: number, members: [string, number][]) => {
    const ancestor = mutate(root, rate);
    for (const [m, r] of members) leaves.push([`${name}${m}`, mutate(ancestor, r)]);
  };
  clade("A", 0.35, [
    ["1", 0.04],
    ["2", 0.04],
    ["3", 0.05],
    ["4", 0.05],
    ["5", 0.3],
  ]);
  const b = mutate(root, 0.35);
  for (const sub of ["1", "2"]) {
    const subAncestor = mutate(b, 0.12);
    for (const m of ["a", "b"]) leaves.push([`B${sub}${m}`, mutate(subAncestor, 0.03)]);
  }
  clade("C", 0.35, [
    ["1", 0.05],
    ["2", 0.05],
    ["3", 0.06],
  ]);
  leaves.push(["O1", mutate(root, 0.8)]);
  const shuffled = leaves.map((_, i) => leaves[(i * 5) % leaves.length]);
  return { labels: shuffled.map(([l]) => `SYN-${l}/1-${root.length}`), rows: shuffled.map(([, r]) => r) };
}
const EVOLVED = evolvedFamily();
/** The evolved family with a tree the viewer computes at load: `nj` or `upgma` rows in tree order. */
export function computedTreeDoc(provider: "nj" | "upgma"): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: `A guide tree by ${provider}` });
  const aln = b.data.alignment({ ref: "family", value: EVOLVED.rows, labels: EVOLVED.labels });
  aln.tree({ provider, options: { distance: "p" } } as never);
  const v = b.view({ ref: "fam", label: "Family", style: "msa", height: "box", height_rows: EVOLVED.rows.length + 1 });
  v.block({ alignment: aln });
  v.track({ label: "" }).ruler();
  v.track({ label: "Members" }).rows({ order: "tree" });
  return b.getDocument();
}
/**
 * A family with insertions and deletions (deterministic, not real sequences): four clades off a root, B with a deletion
 * of 5 residues, C an insertion of GSGS, D a deletion of 4 and a P; each clade's third member loses 3 more. Rows
 * are the family's own residues, padded at the end to one width: an alignment that was never made.
 */
const INDELS = (() => {
  const AMINO = "ACDEFGHIKLMNPQRSTVWY";
  let seed = 4;
  const rnd = () => {
    seed = (Math.imul(seed, 1103515245) + 12345) & 0x7fffffff;
    return seed / 0x7fffffff;
  };
  const sub = (from: string, rate: number) =>
    [...from].map((c) => (rnd() < rate ? AMINO[Math.floor(rnd() * 20)] : c)).join("");
  const edit = (from: string, at: number, drop: number, add: string) => from.slice(0, at) + add + from.slice(at + drop);
  const root = "WYFGKITRRESERLLLNAENPRGTFLVRESWFHPNITGVEKDAENLLLTRGVDGSFLARPWYWGDISREEVNEKLRDTADG";
  const clades = [
    sub(root, 0.08),
    edit(sub(root, 0.08), 20, 5, ""),
    edit(sub(root, 0.08), 45, 0, "GSGS"),
    edit(sub(root, 0.1), 60, 4, "P"),
  ];
  const leaves: [string, string][] = [];
  for (const [c, ancestor] of clades.entries()) {
    for (const m of [1, 2, 3]) {
      const own = sub(ancestor, 0.04);
      leaves.push([`${"ABCD"[c]}${m}`, m === 3 ? edit(own, 10 + c * 3, 3, "") : own]);
    }
  }
  const shuffled = leaves.map((_, i) => leaves[(i * 5) % leaves.length]);
  const width = Math.max(...shuffled.map(([, r]) => r.length));
  return {
    labels: shuffled.map(([l]) => `SYN-${l}/1-${root.length}`),
    rows: shuffled.map(([, r]) => r.padEnd(width, "-")),
  };
})();
/** `INDELS` as given, or realigned by the `realign` node provider (biomsa) with an NJ tree over the new rows. */
export function realignedDoc(realigned: boolean): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: realigned ? "Realigned, with a tree" : "Rows as given" });
  const given = b.data.alignment({ ref: "given", value: INDELS.rows, labels: INDELS.labels });
  const shown = realigned ? b.data.alignment({ ref: "re", provider: "realign", inputs: given } as never) : given;
  if (realigned) shown.tree({ provider: "nj" } as never);
  const v = b.view({ ref: "fam", label: "Family", style: "msa", height: "box", height_rows: INDELS.rows.length + 1 });
  v.block({ alignment: shown });
  v.track({ label: "" }).ruler();
  v.track({ label: "Members" }).rows({ order: "tree" });
  return b.getDocument();
}

/** A large synthetic alignment by value (deterministic letters and gaps, not real sequences) for frame times. */
export function largeAlignment(rows: number, columns: number) {
  const AA = "ACDEFGHIKLMNPQRSTVWY";
  const lines = ["# STOCKHOLM 1.0"];
  for (let k = 0; k < rows; k++) {
    let r = "";
    for (let c = 0; c < columns; c++) {
      const h = Math.imul((k % 97) * 73856093 + c * 19349663 + (k >> 5), 2654435761) >>> 0;
      r += c % 37 === 5 ? "W" : h % 7 === 0 ? "-" : AA[(c * 7 + (h % 5)) % 20];
    }
    lines.push(`R${k}/1-${columns} ${r}`);
  }
  lines.push("//", "");
  const b = createSeqViewSpec({ title: `A synthetic alignment of ${rows} × ${columns}` });
  const aln = b.data.alignment({
    ref: "big",
    source: { value: lines.join("\n"), format: "stockholm", max_bytes: 1 << 23 },
  });
  const v = b.view({ ref: "large", label: "Large", style: "msa", height: "box", height_rows: 24 });
  v.block({ alignment: aln });
  v.track({ label: "" }).ruler();
  v.track({ label: "Rows" }).rows();
  return b.getDocument();
}

/** A titin-sized sequence by value (synthetic), with features every 500 residues. */
export function longSequence(n: number) {
  const AA = "ACDEFGHIKLMNPQRSTVWY";
  const seq = Array.from({ length: n }, (_, i) => AA[(i * 7 + ((i * i) % 11)) % 20]).join("");
  const b = createSeqViewSpec({ title: `A synthetic sequence of ${n} residues` });
  const a = b.data.sequence({ ref: "long", label: "Titin-sized", value: seq });
  const domains = b.data.annotation({
    ref: "domains",
    geometry: "ranges",
    target: a,
    value: Array.from({ length: Math.floor(n / 500) }, (_, k) => ({
      start: k * 500 + 20,
      stop: k * 500 + 110 + (k % 5) * 40,
      label: `Ig ${k + 1}`,
    })),
  });
  const v = b.view({ ref: "long-view", label: "Tracks" }).block({ sequence: a });
  v.track({ label: "" }).ruler();
  v.track({ label: "Residues" }).letters();
  v.track({ label: "Domains" }).ranges({ data: domains });
  return b.getDocument();
}

export const trackRef = (doc: SeqViewSpecDoc, label: string) => {
  const walk = (n: DocNode): string | null =>
    n.kind === "track" && n.params?.label === label
      ? (n.ref ?? null)
      : ((n.children ?? []).map(walk).find(Boolean) ?? null);
  return walk(doc.root) as string;
};
/** Case 2 with snapshots as patches over it. */
export function antibodyWithSnapshots(doc: SeqViewSpecDoc): SeqViewSpecDoc {
  const residues = trackRef(doc, "Residues");
  const snapshot = (ref: string, title: string, description: string, ops: Op[]) => ({ ref, title, description, ops });
  return {
    ...doc,
    metadata: {
      ...doc.metadata,
      description:
        "Every change here is a patch: a header's menu (filter view, hide, remove), the section's switches, the chips. [The CDRs](!seq-highlight=cdrs) are in IMGT numbering.",
    },
    snapshots: [
      snapshot("only-cdrs", "Only the CDRs", "The axis cut to what [the CDRs](!seq-highlight=cdrs) cover.", [
        { op: "set", ref: "fv.main", params: { positions: "cdrs" } },
      ]),
      snapshot("imgt", "In IMGT", "The same residues numbered by IMGT: CDR-H3 is 105–117.", [
        { op: "set", ref: "fv.main", params: { numbering: "imgt" } },
      ]),
      snapshot("no-letters", "Residues hidden", "The letters hidden: a chip shows them again.", [
        { op: "set", ref: residues, params: { hidden: true } },
      ]),
    ],
  };
}

// ---- Seq*/React UI: the same case 2 with plain-text descriptions ----

/** Case 2 with snapshots as patches over it, described in plain text. */
export function antibodyPlain(): SeqViewSpecDoc {
  const doc = case2();
  const residues = trackRef(doc, "Residues");
  return {
    ...doc,
    metadata: {
      ...doc.metadata,
      description:
        "Every change here is a patch: a header's menu (filter view, hide, remove), the section's switches, the chips.",
      description_format: "plaintext",
    },
    snapshots: [
      {
        ref: "only-cdrs",
        title: "Only the CDRs",
        description: "The axis cut to what the CDRs cover.",
        ops: [{ op: "set", ref: "fv.main", params: { positions: "cdrs" } }],
      },
      {
        ref: "imgt",
        title: "In IMGT",
        description: "The same residues numbered by IMGT: CDR-H3 is 105–117.",
        ops: [{ op: "set", ref: "fv.main", params: { numbering: "imgt" } }],
      },
      {
        ref: "no-letters",
        title: "Residues hidden",
        description: "The letters hidden: a chip shows them again.",
        ops: [{ op: "set", ref: residues, params: { hidden: true } }],
      },
    ],
  } as SeqViewSpecDoc;
}

// ---- Seq*/Cyclic ----

/**
 * A synthetic 240-nucleotide plasmid (a committed fixture: the format's load-fixtures.ts) whose replication origin runs
 * across the origin, [200, 30), and three other features; EcoRI and HindIII cut sites. Four sections over it: from the
 * origin, rotated to 200, a piece across the origin, and rotated with a focus on the replication origin.
 */
export function plasmid(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "pDemo, a cyclic plasmid" });
  const p = b.data.sequence({
    ref: "pDemo",
    label: "pDemo",
    alphabet: "dna",
    cyclic: true,
    source: { host_file: "pDemo.fasta", format: "fasta" },
  });
  const feats = b.data.annotation({
    ref: "features",
    geometry: "ranges",
    target: p,
    source: { host_file: "pDemo.tsv", format: "tsv" },
  });
  const cuts = b.data.annotation({
    ref: "cuts",
    geometry: "ranges",
    target: p,
    label: "Restriction sites",
    value: { boundary: [61, 85], label: ["EcoRI", "HindIII"] },
  });
  const scissors = { theme: "uniform", params: { glyph: { char: "✂", label: "cut site" } } };
  const v = b.view({ ref: "map", label: "Linear maps" });
  const tracks = (sec: ReturnType<typeof v.section>) => {
    sec.track({ label: "" }).ruler();
    sec.track({ label: "pDemo" }).letters();
    sec.track({ label: "Features" }).ranges({ data: feats, labels: true });
    sec.track({ label: "Sites" }).ranges({ data: cuts, glyph: scissors });
    return sec;
  };
  tracks(v.section({ ref: "map.origin", label: "From the origin" })).block({ sequence: p });
  tracks(v.section({ ref: "map.rotated", label: "Rotated: origin 200" })).block({ sequence: p, origin: 200 });
  tracks(v.section({ ref: "map.piece", label: "A piece across the origin: range [180, 60]" })).block({
    sequence: p,
    range: [180, 60],
    label: "pDemo 181–60",
  });
  tracks(v.section({ ref: "map.focus", label: "Rotated, focused on the replication origin" }))
    .block({ sequence: p, origin: 150 })
    .focus({ annotation: feats, items: ["ori"] });
  // The same tracks on rings (the `circular` preset): the axis from the top, clockwise, the outer ring the ruler.
  const ring = b.view({ ref: "ring", label: "Ring", style: "circular" }).block({ sequence: p });
  ring.track({ label: "" }).ruler();
  ring.track({ label: "pDemo" }).letters();
  ring.track({ label: "Features" }).ranges({ data: feats, labels: true });
  ring.track({ label: "Sites" }).ranges({ data: cuts, glyph: scissors });
  return b.getDocument();
}

/** The pDemo unit's runs as coordinates (a host's form): [start, stop) pairs. */
export const pDemoRuns = (...xs: number[]) => ({ systems: ["seqstar:pDemo/unit"], runs: Int32Array.from(xs) });

// ---- Seq*/Hover ----

const HOVER_PROTEIN =
  "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNT";
const N = 160;
const at = <T>(n: number, f: (i: number) => T) => Array.from({ length: n }, (_, i) => f(i));

/** One chain with a track per type that has an emphasis; `style` on its section (`hover_emphasis`). */
export function tracksDoc(title: string, style?: Record<string, unknown>): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title });
  const p = b.data.sequence({ ref: "p53", label: "p53", value: HOVER_PROTEIN.slice(0, N) });
  const ann = (ref: string, geometry: string, value: unknown, more: object = {}) =>
    b.data.annotation({ ref, geometry, target: p, value, ...more } as never);
  const regions = ann("regions", "ranges", {
    id: ["tad", "tad", "dbd", "nls", "tet", "pro"],
    start: [0, 40, 60, 100, 120, 20],
    stop: [25, 55, 150, 110, 140, 70],
    label: ["TAD 1", "TAD 2", "DNA binding", "NLS", "Tetramer", "Proline-rich"],
    category: ["transactivation", "transactivation", "binding", "signal", "oligomer", "low complexity"],
  });
  const sites = ann("sites", "ranges", {
    start: [14, 19, 36, 72, 98, 119, 132],
    stop: [15, 20, 37, 73, 99, 120, 133],
    category: ["phospho", "phospho", "acetyl", "phospho", "ubiquitin", "acetyl", "phospho"],
  });
  const variants = ann("variants", "ranges", {
    start: [21, 21, 47, 72, 72, 72, 110, 151],
    stop: [22, 22, 48, 73, 73, 73, 111, 152],
    alt: ["A", "L", "S", "H", "W", "Q", "C", "V"],
  });
  const ss = ann("ss", "ranges", {
    start: [12, 30, 62, 84, 96, 112, 126, 146],
    stop: [24, 36, 70, 92, 104, 120, 136, 156],
    category: ["helix", "strand", "strand", "turn", "strand", "helix", "strand", "helix"],
  });
  const index = at(N, (i) => i);
  const score = ann("score", "values", {
    index,
    value: index.map((i) => Number((0.5 + 0.45 * Math.sin(i / 9)).toFixed(3))),
  });
  const contacts = ann(
    "contacts",
    "matrix",
    { index, value: index.map((i) => [(i % 7) / 7, ((i * 3) % 11) / 11, (i % 5) / 5, ((i * 7) % 13) / 13]) },
    { states: ["ATP", "Mg", "DNA", "Zn"] },
  );
  const v = b.view({ ref: "v", label: title, ...(style ? { style } : {}) } as never);
  v.block({ sequence: p });
  const sec = v;
  sec.track({ label: "" }).ruler();
  sec.track({ label: "Letters" }).letters();
  sec.track({ label: "Regions" }).ranges({ data: regions, labels: true });
  sec.track({ label: "Sites" }).points({ data: sites });
  sec.track({ label: "Variants" }).variants({ data: variants });
  sec.track({ label: "Secondary" }).secondary({ data: ss });
  sec.track({ label: "Score" }).values({ data: score });
  sec.track({ label: "Contacts" }).matrix({ data: contacts });
  return b.getDocument();
}

/** The synthetic family with its profile as a logo (no emphasis over its letters) and a matrix (its cell outlined). */
export function logoDoc(): SeqViewSpecDoc {
  const { labels, rows } = familyRows(12);
  const sto = ["# STOCKHOLM 1.0", ...rows.map((r, k) => `${labels[k].padEnd(14)}${r}`), "//", ""].join("\n");
  const b = createSeqViewSpec({ title: "Family" });
  const aln = b.data.alignment({ ref: "family", source: { value: sto, format: "stockholm" } });
  const profile = b.data.computed({ ref: "profile", provider: "profile", alignment: aln });
  const v = b.view({ ref: "fam", label: "Family", style: "msa", height: "box", height_rows: 6 });
  v.block({ alignment: aln });
  v.track({ label: "Logo" }).logo({ data: profile });
  v.track({ label: "Profile" }).matrix({ data: profile });
  return b.getDocument();
}

// ---- Seq*/Letter maps ----

const GENE =
  "GAATTCGCCACC" + // 5′ UTR
  "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC" + // exon 1's coding part
  "GTAAGTCTTTCCTAACTTCTGTTTCCACAG" + // the intron
  "ATTGGTACTTCAACAGCCAGGTGGAGTAA" + // exon 2's coding part and stop
  "GCGGCCGCATAAAGT"; // 3′ UTR
const CDS: [number, number][] = [
  [12, 46],
  [76, 105],
];

/** Case 14's gene fragment, its protein and exons in a builder: the pieces every card below shows. */
function gene(title: string) {
  const b = createSeqViewSpec({ title });
  const g = b.data.sequence({ ref: "gene", label: "Gene fragment", alphabet: "dna", value: GENE });
  const protein = translate(b.data, { ref: "protein", label: "Protein", from: g, cds: CDS });
  const exons = b.data.annotation({
    ref: "exons",
    geometry: "ranges",
    target: g,
    value: { start: [0, 76], stop: [46, 120], label: ["Exon 1", "Exon 2"] },
  });
  const start = b.data.annotation({
    ref: "start",
    geometry: "ranges",
    target: g,
    value: { start: [12], stop: [15], label: ["Start codon"] },
  });
  return { b, g, protein, exons, start };
}
type Sec = ReturnType<ReturnType<ReturnType<typeof createSeqViewSpec>["view"]>["section"]>;

/** The gene as stored, its reverse complement (the preset) and its mRNA (T shown as U) under its protein. */
export function letterMapStrands(): SeqViewSpecDoc {
  const { b, g, protein, exons, start } = gene("Case 14's gene: strand flip and codons through letter maps");
  const v = b.view({ ref: "maps", label: "Letter maps" });
  const stored = v.section({ ref: "maps.stored", label: "As stored" }).block({ sequence: g });
  stored.track({ label: "" }).ruler();
  stored.track({ label: "DNA" }).letters();
  stored.track({ label: "Exons" }).ranges({ data: exons, labels: true });
  const rc: Sec = v.section({ ref: "maps.rc", label: "Reverse complement", follow: stored.handle });
  rc.block({ sequence: g });
  rc.track({ label: "" }).ruler();
  rc.track({ label: "DNA" }).letters().emphasis({ data: start });
  rc.track({ label: "Exons" }).ranges({ data: exons, labels: true });
  b.preset.reverseComplement(rc);
  const rna = v.section({ ref: "maps.rna", label: "mRNA: T shown as U", follow: stored.handle }).block({ sequence: g });
  rna.track({ label: "mRNA" }).letters({ letter_map: { T: "U" } });
  rna.track({ label: "Protein" }).letters({ of: protein });
  return b.getDocument();
}

/** Only the exons (the preset): `positions` over the exon annotation, the intron collapsed to a break. */
export function exonsOnly(): SeqViewSpecDoc {
  const { b, g, protein, exons } = gene("Case 14's gene, exons only");
  const v = b.view({ ref: "ex", label: "Exons only" });
  const s = v.section({ ref: "ex.main", label: "Exons only (positions over the exons)" }).block({ sequence: g });
  s.track({ label: "" }).ruler();
  s.track({ label: "mRNA" }).letters({ letter_map: { T: "U" } });
  s.track({ label: "Protein" }).letters({ of: protein });
  b.preset.exonsOnly(s, exons);
  return b.getDocument();
}

/** The reverse complement as text, in groups: the strand a person would copy. */
export function asText(): SeqViewSpecDoc {
  const { b, g } = gene("Case 14's gene, the reverse complement as text");
  const v = b.view({ ref: "text", label: "Text", layout: "text" });
  const sec = v.section({ ref: "text.rc", label: "Reverse complement" }).block({ sequence: g });
  sec.track({ label: "DNA" }).letters();
  b.preset.reverseComplement(sec);
  return b.getDocument();
}

// ---- Seq*/Links ----

/** Case 14's gene fragment and its translation, as two sections: the protein follows the gene through the codon step. */
export function codingPair(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "A protein following its coding DNA" });
  const GENE =
    "GAATTCGCCACC" + // 5′ UTR
    "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC" + // exon 1's coding part
    "GTAAGTCTTTCCTAACTTCTGTTTCCACAG" + // the intron
    "ATTGGTACTTCAACAGCCAGGTGGAGTAA" + // exon 2's coding part and stop
    "GCGGCCGCATAAAGT"; // 3′ UTR
  const gene = b.data.sequence({ ref: "gene", label: "Gene fragment", alphabet: "dna", value: GENE });
  const protein = translate(b.data, {
    ref: "protein",
    label: "Protein",
    from: gene,
    cds: [
      [12, 46],
      [76, 105],
    ],
  });
  const exons = b.data.annotation({
    ref: "exons",
    geometry: "ranges",
    target: gene,
    value: { start: [0, 76], stop: [46, 120], label: ["Exon 1", "Exon 2"] },
  });
  const domain = b.data.annotation({
    ref: "domain",
    geometry: "ranges",
    target: protein,
    numbering: "ordinal",
    value: { beg: [5], end: [14], label: ["Domain"] },
  });
  const v = b.view({ ref: "pair", label: "Gene and protein" });
  const dna = v.section({ ref: "pair.dna", label: "Gene" }).block({ sequence: gene });
  dna.track({ label: "" }).ruler();
  dna.track({ label: "DNA" }).letters();
  dna.track({ label: "Exons" }).ranges({ data: exons, labels: true });
  const p = v.section({ ref: "pair.protein", label: "Protein", follow: dna.handle }).block({ sequence: protein });
  p.track({ label: "" }).ruler();
  p.track({ label: "Protein" }).letters();
  p.track({ label: "Domain" }).ranges({ data: domain, labels: true });
  return b.getDocument();
}

/** Case 5's pET insert: its reverse strand follows it, a reversed piece following a forward one. */
export function linkStrands(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "pET insert, the reverse strand following" });
  const ins = b.data.sequence({
    ref: "insert",
    label: "pET insert",
    alphabet: "dna",
    source: { host_file: "pET-insert.fasta", format: "fasta" },
  });
  const feats = b.data.annotation({
    ref: "features",
    geometry: "ranges",
    target: ins,
    source: { host_file: "pET-insert.tsv", format: "tsv" },
  });
  const cuts = b.data.annotation({
    ref: "cuts",
    geometry: "ranges",
    target: ins,
    label: "Restriction sites",
    value: { boundary: [CUT_ECORI, CUT_XHOI], label: ["EcoRI", "XhoI"] },
  });
  const scissors = { theme: "uniform", params: { glyph: { char: "✂", label: "cut site" } } };
  const v = b.view({ ref: "strands", label: "Both strands" });
  const tracks = (sec: ReturnType<typeof v.section>) => {
    sec.track({ label: "" }).ruler();
    sec.track({ label: "Insert" }).letters();
    sec.track({ label: "Features" }).ranges({ data: feats, labels: true });
    sec.track({ label: "Sites" }).ranges({ data: cuts, glyph: scissors });
    return sec;
  };
  const fwd = tracks(v.section({ ref: "strands.forward", label: "As written" }).block({ sequence: ins }));
  tracks(
    v
      .section({ ref: "strands.reverse", label: "Reverse strand", direction: -1, follow: fwd.handle })
      .block({ sequence: ins }),
  );
  return b.getDocument();
}

// ---- Seq*/Pieces ----

/**
 * Case 5's construct (a 72-nucleotide pET insert: a start codon, EcoRI, a GFP fragment, XhoI) on both strands: as
 * written, its reverse strand (`direction: -1`), and the GFP fragment alone, reversed (`range: [9, 63]`); and the
 * sequence reversed as text. The letters are the forward strand's, read right to left: complementing them is a letter
 * map's.
 */
export function flipped(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "pET insert, both strands" });
  const ins = b.data.sequence({
    ref: "insert",
    label: "pET insert",
    alphabet: "dna",
    source: { host_file: "pET-insert.fasta", format: "fasta" },
  });
  const feats = b.data.annotation({
    ref: "features",
    geometry: "ranges",
    target: ins,
    source: { host_file: "pET-insert.tsv", format: "tsv" },
  });
  const cuts = b.data.annotation({
    ref: "cuts",
    geometry: "ranges",
    target: ins,
    label: "Restriction sites",
    value: { boundary: [CUT_ECORI, CUT_XHOI], label: ["EcoRI", "XhoI"] },
  });
  const scissors = { theme: "uniform", params: { glyph: { char: "✂", label: "cut site" } } };
  const v = b.view({ ref: "map", label: "Both strands" });
  const tracks = (sec: ReturnType<typeof v.section>) => {
    sec.track({ label: "" }).ruler();
    sec.track({ label: "Insert" }).letters();
    sec.track({ label: "Features" }).ranges({ data: feats, labels: true });
    sec.track({ label: "Sites" }).ranges({ data: cuts, glyph: scissors });
  };
  tracks(v.section({ ref: "map.forward", label: "As written" }).block({ sequence: ins }));
  tracks(v.section({ ref: "map.reverse", label: "Reverse strand", direction: -1 }).block({ sequence: ins }));
  tracks(
    v
      .section({ ref: "map.gfp", label: "The GFP fragment, reversed" })
      .block({ sequence: ins, range: [9, 63], direction: -1, label: "GFP fragment" }),
  );
  const t = b.view({ ref: "seq", label: "Sequence, reverse strand", layout: "text", direction: -1 });
  t.block({ sequence: ins });
  t.track({ label: "Sequence" }).letters();
  t.track({ label: "Features" }).ranges({ data: feats, display: "underline" });
  return b.getDocument();
}

// ---- Seq*/Style ----

/** `doc` with the section `ref`'s `style` set (a copy). */
export function styled(doc: SeqViewSpecDoc, ref: string, style: string | SectionStyle | Record<string, unknown>) {
  const walk = (n: SeqViewSpecDoc["root"]): SeqViewSpecDoc["root"] => ({
    ...n,
    params: n.ref === ref ? { ...n.params, style } : n.params,
    children: n.children?.map(walk),
  });
  return { ...doc, root: walk(doc.root) };
}

export const HEADS: SectionStyle[] = [
  { preset: "msa", column_heads: "shared" },
  { preset: "msa", column_heads: "own_row" },
];

export const STYLE_FAMILY = familyWithTree(60, 12, { columns: true });
export const FLOORS: SectionStyle[] = [
  { preset: "msa", data_floor: 0 },
  { preset: "msa" },
  { preset: "msa", data_floor: 160 },
];

export const BAD = { preset: "msa", colum_heads: "own_row", data_floor: -10, row_height: 26 };

// ---- Seq*/Query: case 7 with a barnase family ----

/**
 * Case 7's 1BRS with a barnase family bound to chain A's entity (the eval's scenario 25; synthetic: barnase's mature
 * chain and eleven deterministic variants) and its conservation, in two views: the complex and the family. Results
 * `onto` the alignment lie in the Family view, so "Show me" switches to it.
 */
export function familyDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "Barnase–barstar with a barnase family" });
  const s = b.data.complex({
    ref: "1brs",
    provider: "molstar:structure",
    source: { service: "pdbe", resource: "structure", id: "1brs" },
  });
  const contacts = b.data.annotation({ ref: "contacts", geometry: "pairs", numbering: "auth", value: contactsAD() });
  const AMINO = "ACDEFGHIKLMNPQRSTVWY";
  const rows = Array.from({ length: 12 }, (_, k) =>
    [...BARNASE_MATURE].map((c, i) => (k > 0 && (i * 7 + k * 13) % 9 < k % 4 ? AMINO[(i + k) % 20] : c)).join(""),
  );
  const sto = ["# STOCKHOLM 1.0", ...rows.map((r, k) => `${k ? `HOM${String(k).padStart(2, "0")}` : "BARNASE"}  ${r}`)];
  const msa = b.data.alignment({ ref: "msa", source: { value: [...sto, "//", ""].join("\n"), format: "stockholm" } });
  msa.bind(0, s.sequence("entity:1"), { method: "align" });
  const cons = b.data.computed({ ref: "cons", provider: "column_stats", alignment: msa, measure: "conservation" });
  const v = b
    .view({ ref: "complex", label: "Complex" })
    .block({ unit: s.unit("A") })
    .block({ unit: s.unit("D") });
  v.track({ label: "" }).ruler();
  v.track({ label: "Residues" }).letters();
  v.track({ label: "Contacts" }).pairs({ data: contacts });
  const f = b.view({ ref: "family", label: "Family" }).block({ alignment: msa });
  f.track({ label: "Conservation" }).values({ data: cons });
  f.track({ label: "Rows" }).rows();
  return b.getDocument();
}

// ---- Seq*/DOM UI: case 10 with markdown descriptions ----

const T315 = { beg: "315", end: "315", numbering: "auth" };
/** Case 10's own descriptions, written for this story: command links of every key of the link table. */
const DOCUMENT = [
  "# The Stories drafts",
  "",
  `ABL kinase with imatinib. Hover [the gatekeeper, T315](${formatLink("seq-highlight", "1iep/A", T315)}) to light it in every view,`,
  `click [a stretch of the chain](${formatLink("seq-select", "1iep/A", { start: 10, stop: 20 })}) to select it, or`,
  `[focus the gatekeeper site](${formatLink("seq-focus", "site.gatekeeper", { items: "0" })}).`,
  "",
  `Next: [the alignment](${formatLink("seq-view", "alignment")}), or [start the story](!seq-snapshot=next).`,
  "",
  "A link to [an ordinary page](https://www.rcsb.org/structure/1IEP) opens as any link, and a [script link](javascript:alert(1)) is dropped by the sanitizer.",
].join("\n");
const SNAPSHOT = `This step shows [the pocket](${formatLink("seq-highlight", "1iep/A", T315)}) only; the descriptions change with the snapshot and release the links of the last one. [Next](!seq-snapshot=next).`;

/** Case 10 with the descriptions the story shows. */
export function describedCase10(): SeqViewSpecDoc {
  const doc = case10();
  return {
    ...doc,
    metadata: { ...doc.metadata, description: DOCUMENT, description_format: "markdown" },
    snapshots: doc.snapshots?.map((s) =>
      s.ref === "s2" ? { ...s, description: SNAPSHOT, description_format: "markdown" } : s,
    ),
  };
}

// ---- Seq*/Search: location sets ----

export const LOCATION_SECTION = "tracks.main";
export const BOOKMARKS = bookmarksOf(LOCATION_SECTION);
/** Case 18 (p53's queries, the Search stories' names, bookmarks and `highlights` entries): the Location sets story's. */
export const locationSetsDoc = (): SeqViewSpecDoc => CASES.case18();
