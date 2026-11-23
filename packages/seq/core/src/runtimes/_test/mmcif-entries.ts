/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { validateAndNormalize } from '@molstar/core/repr/schema';
import type { CategoryColumn, Node } from '@molstar/core/repr/tree';
import { parseCifText } from '@molstar/io/reader/cif/text/parser';
import type { CifFrame } from '@molstar/io/reader/cif/data-model';
import { checkDataNode, SeqLanguage, SeqLanguageVersion, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { expect } from 'vitest';
import { DefaultOpenLimits, type OpenedSource } from '../../data/protocol.js';
import { type MmcifOptions, openMmcif } from '../mmcif.js';

// The mmCIF runtime's conformance suite's entries (_spec/README.md "Data", mmCIF): the PDB files recorded by
// `examples/scripts/record-structures.ts`, a summary of the nodes the runtime reads from each, the counts recorded from
// each file (read by hand from its categories, 2026-11-20), and the same counts counted again from the file's own
// categories by `@molstar/io`'s parser alone (no code of the runtime's), so a count is held twice.

export type EntryId = '1EHZ' | '1FIN' | '4HHB' | '3SGJ';

export const fixtureBytes = (id: EntryId) =>
  new Uint8Array(
    readFileSync(fileURLToPath(new URL(`../../../../examples/fixtures/structures/${id}.cif.gz`, import.meta.url))),
  );

export const openEntry = (id: EntryId, options: MmcifOptions = {}) =>
  openMmcif(
    fixtureBytes(id),
    { signal: new AbortController().signal, limits: DefaultOpenLimits },
    { address: `/data/${id.toLowerCase()}`, format: 'mmcif' },
    options,
  );

const values = (c: CategoryColumn) => Array.from({ length: c.length }, (_, i) => c.get(i));
const runs = (o: { start: Int32Array; end: Int32Array } | undefined) =>
  o ? Array.from(o.start, (s, k) => `${s}-${o.end[k]}`).join(',') : 'all';

/** What an entry reads as, in a few words a node. */
export interface Summary {
  readonly sequences: Record<string, string>;
  readonly molecules: Record<string, string>;
  readonly compounds: number;
  /** A unit: `of <entity> auth <chain> observed <runs> labels <first>..<last> (<n>)`. */
  readonly units: Record<string, string>;
  readonly links: Record<string, number>;
  readonly issues: readonly string[];
}

export function summary(source: OpenedSource): Summary {
  const out = { sequences: {}, molecules: {}, compounds: 0, units: {}, links: {} } as {
    -readonly [K in keyof Summary]: Summary[K];
  };
  for (const { node } of source.coverage.nodes ?? []) {
    const p = (node.params ?? {}) as Record<string, any>;
    if (node.kind === 'sequence')
      out.sequences[node.id!] = `${p.alphabet} ${p.letters.length} modified ${p.modified?.position.length ?? 0}`;
    else if (node.kind === 'molecule') {
      const own = node.children?.find((c) => c.kind === 'pairs')?.params?.rows as { a_start: Int32Array } | undefined;
      out.molecules[node.id!] =
        `${p.class} [${values(p.residues).join(' ')}] atoms ${p.atoms?.residue.length ?? 0} ids ${
          p.identifiers ? values(p.identifiers.namespace).join(',') : '-'
        }${own ? ` links ${own.a_start.length}` : ''}`;
    } else if (node.kind === 'compounds') out.compounds = node.children?.length ?? 0;
    else if (node.kind === 'complex')
      for (const c of node.children ?? []) {
        const q = (c.params ?? {}) as Record<string, any>;
        if (c.kind === 'unit') {
          const labels = (c.children?.[0]?.params?.labels as CategoryColumn | undefined)
            ? values(c.children![0].params!.labels as CategoryColumn).filter((x) => x !== null)
            : [];
          out.units[c.id!] =
            `of ${q.of.split(' ').pop().slice(0, -1).replace(/"/g, '')} auth ${q.auth_chain} observed ${runs(q.observed)} labels ${labels[0]}..${labels[labels.length - 1]} (${labels.length})`;
        } else if (c.kind === 'pairs') for (const t of values(q.rows.type)) out.links[t!] = (out.links[t!] ?? 0) + 1;
      }
  }
  return { ...out, issues: (source.issues ?? []).map((i) => i.code) };
}

/** The nodes pass the language's description and its data rules, as the overlay checks a source's. */
export function expectValid(nodes: readonly { node: Node }[]) {
  for (const { node } of nodes) {
    const spec: Node = {
      kind: 'spec',
      params: { version: SeqLanguageVersion },
      children: [{ kind: 'data', id: 'data', children: [node] }],
    };
    expect(validateAndNormalize(SeqLanguage, spec, SeqValidateOptions).issues).toEqual([]);
    const said: string[] = [];
    checkDataNode(node, { target: () => undefined }, (code) => said.push(code));
    expect(said).toEqual([]);
  }
}

/** The counts recorded from each file, by hand from its categories (the deposited model, no waters). */
export const Recorded: Readonly<Record<EntryId, Summary>> = {
  '1EHZ': {
    // tRNA-Phe: 76 nucleotides, 14 of them modified (`pdbx_struct_mod_residue`, Y37 is YYG); 11 metal ions on the
    // chain (its other 36 `metalc` rows bind waters), its 25 `covale` rows the backbone at modified residues
    sequences: { '1': 'rna 76 modified 14' },
    molecules: { '2': 'ion [MG] atoms 1 ids ccd', '3': 'ion [MN] atoms 1 ids ccd' },
    compounds: 17,
    units: {
      A: 'of 1 auth A observed all labels 1..76 (76)',
      B: 'of 2 auth A observed all labels 590..590 (1)',
      C: 'of 3 auth A observed all labels 530..530 (1)',
      D: 'of 2 auth A observed all labels 510..510 (1)',
      E: 'of 3 auth A observed all labels 520..520 (1)',
      F: 'of 2 auth A observed all labels 540..540 (1)',
      G: 'of 2 auth A observed all labels 560..560 (1)',
      H: 'of 2 auth A observed all labels 570..570 (1)',
      I: 'of 2 auth A observed all labels 580..580 (1)',
      J: 'of 3 auth A observed all labels 550..550 (1)',
    },
    links: { metal: 11 },
    issues: ['seq.io.mmcif.links'],
  },
  '1FIN': {
    // CDK2 (298) and cyclin A (260, numbered from 173) twice, an ATP on each CDK2; no `struct_conn`
    sequences: { '1': 'protein 298 modified 0', '2': 'protein 260 modified 0' },
    molecules: { '3': 'ligand [ATP] atoms 31 ids ccd' },
    compounds: 21,
    units: {
      A: 'of 1 auth A observed all labels 1..298 (298)',
      B: 'of 2 auth B observed all labels 173..432 (260)',
      C: 'of 1 auth C observed all labels 1..298 (298)',
      D: 'of 2 auth D observed all labels 173..432 (260)',
      E: 'of 3 auth A observed all labels 299..299 (1)',
      F: 'of 3 auth C observed all labels 299..299 (1)',
    },
    links: {},
    issues: [],
  },
  '4HHB': {
    // α (141) and β (146) twice, a heme in each (43 heavy atoms, each Fe on its proximal His), two phosphates (the
    // P alone modeled: a ligand still, by its formula `O4 P -3`)
    sequences: { '1': 'protein 141 modified 0', '2': 'protein 146 modified 0' },
    molecules: { '3': 'ligand [HEM] atoms 43 ids ccd', '4': 'ligand [PO4] atoms 1 ids ccd' },
    compounds: 21,
    units: {
      A: 'of 1 auth A observed all labels 1..141 (141)',
      B: 'of 2 auth B observed all labels 1..146 (146)',
      C: 'of 1 auth C observed all labels 1..141 (141)',
      D: 'of 2 auth D observed all labels 1..146 (146)',
      E: 'of 3 auth A observed all labels 142..142 (1)',
      F: 'of 4 auth B observed all labels 147..147 (1)',
      G: 'of 3 auth B observed all labels 148..148 (1)',
      H: 'of 3 auth C observed all labels 142..142 (1)',
      I: 'of 4 auth D observed all labels 147..147 (1)',
      J: 'of 3 auth D observed all labels 148..148 (1)',
    },
    links: { metal: 4 },
    issues: [],
  },
  '3SGJ': {
    // an IgG1 Fc (two chains, ends unmodeled) and FcγRIIIa, three N-glycans as branched entities (8, 5 and 2
    // residues, their 12 `pdbx_entity_branch_link` rows each molecule's links), 4 N-glycosylation and 7 disulfide rows
    sequences: { '1': 'protein 225 modified 0', '2': 'protein 204 modified 0' },
    molecules: {
      '3': 'branched [NAG NAG BMA MAN NAG MAN NAG FUC] atoms 99 ids linucs,wurcs links 7',
      '4': 'branched [NAG NAG BMA MAN MAN] atoms 61 ids linucs,wurcs links 4',
      '5': 'branched [NAG NAG] atoms 28 ids linucs,wurcs links 1',
      '6': 'ligand [MLI] atoms 7 ids ccd',
    },
    compounds: 25,
    units: {
      A: 'of 1 auth A observed 3-222 labels 226..444 (219)',
      B: 'of 1 auth B observed 4-221 labels 227..443 (217)',
      C: 'of 2 auth C observed 5-174 labels 6..174 (169)',
      D: 'of 3 auth D observed all labels 1..8 (8)',
      E: 'of 3 auth E observed all labels 1..8 (8)',
      F: 'of 4 auth F observed all labels 1..5 (5)',
      G: 'of 5 auth G observed all labels 1..2 (2)',
      H: 'of 6 auth A observed all labels 509..509 (1)',
    },
    links: { disulfide: 7, glycosidic: 4 },
    issues: [],
  },
};

/** The first data block of an entry, by `@molstar/io`'s parser alone. */
export async function frameOf(id: EntryId): Promise<CifFrame> {
  const parsed = await parseCifText(gunzipSync(fixtureBytes(id))).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result.blocks[0];
}

/**
 * The counts read again from the file's categories by the plainest rule each states, none of the converter's: the
 * entities and `struct_asym` chains not water, each polymer's `entity_poly_seq` rows, the `struct_conn` rows of each
 * `conn_type_id` (with how many have an end on a water chain and how many both ends on one chain), each branched
 * entity's `pdbx_entity_branch_link` rows, and each chain's author numbers from its scheme
 * (`pdbx_poly_seq_scheme.auth_seq_num` and `pdb_ins_code` where modeled, `pdbx_branch_scheme`'s and
 * `pdbx_nonpoly_scheme`'s `pdb_seq_num`) as `<first>..<last> (<n>)`. Where they differ from what the runtime reads,
 * the tests say why.
 */
export function countedFromFile(frame: CifFrame) {
  const col = (n: string, f: string) => {
    const c = frame.categories[n];
    const field = c?.getField(f);
    return c ? Array.from({ length: c.rowCount }, (_, i) => (field?.valueKind(i) === 0 ? field.str(i) : '')) : [];
  };
  const type = new Map(col('entity', 'id').map((id, i) => [id, col('entity', 'type')[i]]));
  const asymEntity = new Map(col('struct_asym', 'id').map((id, i) => [id, col('struct_asym', 'entity_id')[i]]));
  const water = (asym: string) => type.get(asymEntity.get(asym) ?? '') === 'water';
  const chains = [...asymEntity.keys()].filter((a) => !water(a));
  const lengths: Record<string, number> = {};
  for (const e of col('entity_poly_seq', 'entity_id')) lengths[e] = (lengths[e] ?? 0) + 1;
  const conn: Record<string, { rows: number; onWater: number; oneChain: number }> = {};
  const [a, b] = ['ptnr1_label_asym_id', 'ptnr2_label_asym_id'].map((f) => col('struct_conn', f));
  col('struct_conn', 'conn_type_id').forEach((t, i) => {
    const c = (conn[t] ??= { rows: 0, onWater: 0, oneChain: 0 });
    c.rows++;
    if (water(a[i]) || water(b[i])) c.onWater++;
    if (a[i] === b[i]) c.oneChain++;
  });
  const branchLinks: Record<string, number> = {};
  for (const e of col('pdbx_entity_branch_link', 'entity_id')) branchLinks[e] = (branchLinks[e] ?? 0) + 1;
  const auth = new Map<string, string[]>();
  const scheme = (cat: string, asymF: string, numF: string, insF?: string) => {
    const [asym, num] = [col(cat, asymF), col(cat, numF)];
    const ins = insF ? col(cat, insF) : [];
    asym.forEach((x, i) => {
      if (!num[i] || water(x)) return;
      if (!auth.has(x)) auth.set(x, []);
      auth.get(x)!.push(`${num[i]}${ins[i] ?? ''}`);
    });
  };
  scheme('pdbx_poly_seq_scheme', 'asym_id', 'auth_seq_num', 'pdb_ins_code');
  scheme('pdbx_branch_scheme', 'asym_id', 'pdb_seq_num');
  scheme('pdbx_nonpoly_scheme', 'asym_id', 'pdb_seq_num', 'pdb_ins_code');
  const labels: Record<string, string> = {};
  for (const [x, l] of auth) labels[x] = `${l[0]}..${l[l.length - 1]} (${l.length})`;
  return {
    entities: [...type.values()].filter((t) => t !== 'water').length,
    chains: chains.length,
    lengths,
    conn,
    branchLinks,
    labels,
  };
}

/** The author numbers a summary's units read: `<first>..<last> (<n>)` by chain. */
export const unitLabels = (s: Summary) =>
  Object.fromEntries(Object.entries(s.units).map(([u, d]) => [u, /labels (.*)$/.exec(d)![1]]));

/**
 * The links the runtime reads, against the file's plain counts, each difference said: a `metalc` row with an end on a
 * water is not read (waters are not read); a `covale` row with both ends on one chain is the chain itself (a polymer's
 * backbone at a modified residue) or a glycan's own bond, which is read from `pdbx_entity_branch_link` instead, as
 * the molecule's links; `hydrog` and the other types are not links (_spec/README.md "Data", mmCIF).
 */
export function expectLinksAsCounted(s: Summary, counted: ReturnType<typeof countedFromFile>) {
  const c = (t: string) => counted.conn[t] ?? { rows: 0, onWater: 0, oneChain: 0 };
  expect(s.links.metal ?? 0).toBe(c('metalc').rows - c('metalc').onWater);
  expect(s.links.disulfide ?? 0).toBe(c('disulf').rows - c('disulf').onWater);
  expect((s.links.covalent ?? 0) + (s.links.glycosidic ?? 0)).toBe(
    c('covale').rows - c('covale').oneChain - c('covale').onWater,
  );
  for (const [e, d] of Object.entries(s.molecules))
    expect(/ links (\d+)$/.exec(d)?.[1] ?? '0', `molecule ${e}`).toBe(String(counted.branchLinks[e] ?? 0));
}
