/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Model } from '@molstar/model/model/structure';

// Dynamic PDB's data (dynamicpdb.com, a JSON:API at /api/v1): an entry, its models (the deposited one, qFit's
// multiconformer one, …) with their metrics and per-residue data, and each model's files. Its API answers only its own
// origin (CORS) and its file bucket none, so the story reads a copy (`fixtures/`, the API's responses as they are, and
// the qFit model as BCIF); `DpdbSource` is where a live client would go once Dynamic PDB allows another origin.

export interface DpdbMutation {
  /** As written: L192F. */
  readonly text: string;
  readonly from: string;
  /** In the entry's UniProt numbering. */
  readonly position: number;
  readonly to: string;
}

export interface DpdbEntity {
  readonly description: string;
  readonly construct?: string;
  readonly organism?: string;
  readonly uniprot?: string;
  readonly mutations: readonly DpdbMutation[];
}

export interface DpdbEntry {
  readonly id: string;
  readonly title: string;
  readonly method: string;
  readonly resolution?: number;
  readonly spaceGroup?: string;
  readonly pdb?: string;
  readonly published?: string;
  readonly entities: readonly DpdbEntity[];
}

export interface DpdbResidue {
  readonly asym: string;
  readonly authAsym: string;
  readonly seq: number;
  readonly auth: number;
  readonly comp: string;
  readonly bIso: number;
  readonly occupancy: number;
  readonly rscc: number;
  /** UniProt accession and position, where SIFTS maps it. */
  readonly uniprot?: { readonly accession: string; readonly position: number };
}

export interface DpdbModel {
  readonly id: string;
  readonly title: string;
  /** Deposited, Multiconformer, … */
  readonly modelType: string;
  readonly purpose?: string;
  readonly atoms?: number;
  readonly residues?: number;
  readonly altlocFraction?: number;
  readonly ligands: readonly string[];
  /** By key: r_work, r_free, clashscore, ramachandran_outliers. */
  readonly metrics: Readonly<Record<string, number>>;
  readonly residueData: readonly DpdbResidue[];
  readonly primaryArtifact: string;
}

export interface DpdbArtifact {
  readonly id: string;
  readonly name: string;
  /** model, structure_factors, … */
  readonly type: string;
  readonly format: string;
  readonly size: number;
  readonly uri: string;
}

/** Where the data comes from: the API's three responses and a model file's address the browser may read. */
export interface DpdbSource {
  entry(id: string): Promise<unknown>;
  models(entry: string): Promise<unknown>;
  artifacts(model: string): Promise<unknown>;
  /** The address to read a model file from, and its format. */
  file(artifact: DpdbArtifact): { readonly url: string; readonly format: string };
}

/** The entry page on dynamicpdb.com. */
export const dpdbPage = (entry: string) => `https://dynamicpdb.com/entries/${entry}`;

type Json = Record<string, any>;

const num = (v: unknown) => (typeof v === 'number' && Number.isFinite(v) ? v : undefined);
const str = (v: unknown) => (typeof v === 'string' && v ? v : undefined);

/** `L192F, F194Y` (UniProt numbering); what does not read so is left out. */
export function parseMutations(text: string | undefined): DpdbMutation[] {
  const out: DpdbMutation[] = [];
  for (const part of (text ?? '').split(/[\s,;]+/)) {
    const m = /^([A-Z])(\d+)([A-Z])$/.exec(part);
    if (m) out.push({ text: part, from: m[1], position: +m[2], to: m[3] });
  }
  return out;
}

export function parseEntry(doc: unknown): DpdbEntry {
  const d = (doc as Json)?.data;
  if (d?.type !== 'entries') throw new Error('Not a Dynamic PDB entry');
  const a: Json = d.attributes ?? {};
  return {
    id: d.id,
    title: str(a.title)?.replace(/'$/, '') ?? d.id,
    method: str(a.method) ?? '',
    resolution: num(a.resolution),
    spaceGroup: str(a.space_group),
    pdb: str(a.external_refs?.pdb),
    published: str(a.published_at),
    entities: (a.polymer_entities ?? []).map(
      (e: Json): DpdbEntity => ({
        description: str(e.description) ?? '',
        construct: str(e.construct),
        organism: str(e.source_organisms?.[0]?.scientific_name),
        uniprot: str(e.uniprot_mappings?.[0]?.accession),
        mutations: parseMutations(e.mutations),
      }),
    ),
  };
}

function parseResidue(r: Json): DpdbResidue {
  const [accession, position] = String(r.uniprot_position ?? '').split(':');
  return {
    asym: r.label_asym_id,
    authAsym: r.auth_asym_id,
    seq: r.label_seq_id,
    auth: r.auth_seq_id,
    comp: r.label_comp_id,
    bIso: r.b_iso,
    occupancy: r.occupancy,
    rscc: r.rscc,
    ...(accession && position ? { uniprot: { accession, position: +position } } : {}),
  };
}

/** The models, the deposited one first. */
export function parseModels(doc: unknown): DpdbModel[] {
  const list: Json[] = (doc as Json)?.data ?? [];
  const models = list.map((d): DpdbModel => {
    const a: Json = d.attributes ?? {};
    const m: Json = a.metadata ?? {};
    return {
      id: d.id,
      title: str(a.title) ?? d.id,
      modelType: str(m.model_type) ?? '',
      purpose: str(m.purpose),
      atoms: num(m.atom_count),
      residues: num(m.modeled_residues),
      altlocFraction: num(m.altloc_fraction),
      ligands: m.ligands ?? [],
      metrics: Object.fromEntries((a.metrics ?? []).map((x: Json) => [x.key, x.value])),
      residueData: (a.residue_data ?? []).map(parseResidue),
      primaryArtifact: a.primary_artifact_id,
    };
  });
  const deposited = (x: DpdbModel) => (x.modelType === 'Deposited' ? 0 : 1);
  return models.sort((x, y) => deposited(x) - deposited(y) || x.id.localeCompare(y.id));
}

export function parseArtifacts(doc: unknown): DpdbArtifact[] {
  return ((doc as Json)?.data ?? []).map((d: Json): DpdbArtifact => {
    const a: Json = d.attributes ?? {};
    return { id: d.id, name: a.name, type: a.type, format: a.format, size: a.size_bytes, uri: a.uri };
  });
}

/** A model file's format by its name. */
export const formatOf = (name: string) => (/\.bcif$/i.test(name) ? 'bcif' : /\.pdb$/i.test(name) ? 'pdb' : 'mmcif');

const fixture = (name: string) => new URL(`./fixtures/${name}`, import.meta.url).href;

async function json(url: string) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return r.json();
}

/**
 * The copy in `fixtures/`: entry dpdb_gijqn2no (3K7F, FAAH). The deposited model is read from RCSB, as the API says;
 * qFit's from the copy, as its bucket has no CORS.
 */
export const DpdbFixtures: DpdbSource = {
  entry: () => json(fixture('entry.json')),
  models: () => json(fixture('models.json')),
  artifacts: (model) => json(fixture(`artifacts-m${model.slice(-3)}.json`)),
  file: (a) =>
    a.uri.startsWith('https://files.dynamicpdb.com/')
      ? { url: fixture('3k7f-qfit.bcif'), format: 'bcif' }
      : { url: a.uri, format: formatOf(a.name) },
};

/** The live API: needs Dynamic PDB to allow the page's origin (the API) and a CORS rule on its file bucket. */
export function dpdbLive(base = 'https://dynamicpdb.com/api/v1'): DpdbSource {
  return {
    entry: (id) => json(`${base}/entries/${id}`),
    models: (id) => json(`${base}/entries/${id}/models`),
    artifacts: (model) => json(`${base}/models/${model}/artifacts`),
    file: (a) => ({ url: a.uri, format: formatOf(a.name) }),
  };
}

// --- tracks ----------------------------------------------------------------------------------------------------------

/**
 * A residue's place in a model: its chain (label_asym_id, the sequence's unit) and its author number (the units' `auth`
 * numbering). By author number, as the models do not share label_seq_id (qFit's starts at 1 where the deposited's has 6
 * unmodeled residues before it).
 */
export interface ResidueKey {
  readonly authAsym: string;
  readonly auth: number;
}
const keyOf = (k: ResidueKey) => `${k.authAsym}:${k.auth}`;

/** A model's polymer residues: author chain and number to label chain (the sequence's unit). */
export function labelChains(model: Model): Map<string, string> {
  const h = model.atomicHierarchy;
  const out = new Map<string, string>();
  for (let c = 0; c < h.chains._rowCount; c++) {
    const start = h.chainAtomSegments.offsets[c];
    const end = h.chainAtomSegments.offsets[c + 1];
    const asym = h.chains.label_asym_id.value(c);
    const authAsym = h.chains.auth_asym_id.value(c);
    for (let r = h.residueAtomSegments.index[start]; r <= h.residueAtomSegments.index[end - 1]; r++) {
      out.set(keyOf({ authAsym, auth: h.residues.auth_seq_id.value(r) }), asym);
    }
  }
  return out;
}

/** A values table (Seq*'s `values` geometry) in the units' `auth` numbering. */
export interface ValuesTable {
  readonly target: string[];
  readonly pos: string[];
  readonly value: number[];
}

/** Values by residue as a table over the model's units: the residues the model has. */
export function valuesTable(
  model: Map<string, string>,
  rows: readonly (ResidueKey & { readonly value: number })[],
): ValuesTable {
  const t: ValuesTable = { target: [], pos: [], value: [] };
  for (const r of rows) {
    const asym = model.get(keyOf(r));
    if (asym === undefined || !Number.isFinite(r.value)) continue;
    t.target.push(asym);
    t.pos.push(`${r.auth}`);
    t.value.push(r.value);
  }
  return t;
}

export type ResidueField = 'bIso' | 'occupancy' | 'rscc';

export const residueValues = (residues: readonly DpdbResidue[], field: ResidueField) =>
  residues.map((r) => ({ authAsym: r.authAsym, auth: r.auth, value: r[field] }));

/** A ranges table (Seq*'s `ranges` geometry) for the `variants` layer: the mutations where the residues are. */
export interface VariantsTable {
  readonly target: string[];
  readonly beg: string[];
  readonly end: string[];
  readonly alt: string[];
  readonly label: string[];
}

/** The entity's mutations on every chain that has the residue (by its UniProt position in the residue data). */
export function mutationsTable(
  model: Map<string, string>,
  residues: readonly DpdbResidue[],
  mutations: readonly DpdbMutation[],
): VariantsTable {
  const t: VariantsTable = { target: [], beg: [], end: [], alt: [], label: [] };
  const at = new Map(mutations.map((m) => [m.position, m]));
  for (const r of residues) {
    const m = r.uniprot && at.get(r.uniprot.position);
    const asym = m && model.get(keyOf(r));
    if (!m || asym === undefined) continue;
    t.target.push(asym);
    t.beg.push(`${r.auth}`);
    t.end.push(`${r.auth}`);
    t.alt.push(m.to);
    t.label.push(m.text);
  }
  return t;
}

/** A mutation's residues as a target (author numbering, every chain that has it): `A:192, B:192`; empty: none. */
export function mutationTarget(residues: readonly DpdbResidue[], m: DpdbMutation): string {
  return residues
    .filter((r) => r.uniprot?.position === m.position)
    .map((r) => `${r.authAsym}:${r.auth}`)
    .join(', ');
}

/** A residue's alternate conformations in a model, and how far its Cα moves from another model's (Å, the most). */
export interface Conformers extends ResidueKey {
  readonly count: number;
  readonly spread: number;
}

/** Cα positions by residue key and altloc ('' for none). */
function alphaCarbons(
  model: Model,
): Map<string, { readonly alts: Set<string>; readonly ca: [number, number, number][] }> {
  const h = model.atomicHierarchy;
  const { x, y, z } = model.atomicConformation;
  const out = new Map<string, { alts: Set<string>; ca: [number, number, number][] }>();
  for (let c = 0; c < h.chains._rowCount; c++) {
    const authAsym = h.chains.auth_asym_id.value(c);
    for (let i = h.chainAtomSegments.offsets[c]; i < h.chainAtomSegments.offsets[c + 1]; i++) {
      const r = h.residueAtomSegments.index[i];
      const key = keyOf({ authAsym, auth: h.residues.auth_seq_id.value(r) });
      let e = out.get(key);
      if (!e) out.set(key, (e = { alts: new Set(), ca: [] }));
      const alt = h.atoms.label_alt_id.value(i);
      if (alt) e.alts.add(alt);
      if (h.atoms.label_atom_id.value(i) === 'CA') e.ca.push([x[i], y[i], z[i]]);
    }
  }
  return out;
}

/** Per residue of `multi`: its conformations, and its Cα's largest distance from `reference`'s (the same crystal frame). */
export function conformers(multi: Model, reference: Model): Conformers[] {
  const ref = alphaCarbons(reference);
  const out: Conformers[] = [];
  for (const [key, e] of alphaCarbons(multi)) {
    const at = ref.get(key)?.ca[0];
    let spread = 0;
    if (at) {
      for (const p of e.ca) spread = Math.max(spread, Math.hypot(p[0] - at[0], p[1] - at[1], p[2] - at[2]));
    }
    const sep = key.lastIndexOf(':');
    out.push({
      authAsym: key.slice(0, sep),
      auth: +key.slice(sep + 1),
      count: Math.max(1, e.alts.size),
      spread,
    });
  }
  return out;
}

/** A metric's place on its bar, and which side is better: lower for every metric here. */
export const MetricScales: Readonly<
  Record<string, { readonly label: string; readonly max: number; readonly digits: number }>
> = {
  r_work: { label: 'R-work', max: 0.35, digits: 3 },
  r_free: { label: 'R-free', max: 0.35, digits: 3 },
  clashscore: { label: 'Clashscore', max: 20, digits: 2 },
  ramachandran_outliers: { label: 'Rama. outliers', max: 2, digits: 1 },
};
