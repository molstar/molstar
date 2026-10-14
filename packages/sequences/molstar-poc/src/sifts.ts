/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SIFTS from the file's own categories (seqspec:libraries#the-mol-adapter): PDBe's updated mmCIF and BinaryCIF carry
// `pdbx_sifts_unp_segments`, each row an entity instance's segment onto a UniProt accession. An entity gets one
// `uniprot:<acc>` numbering per accession, a `segments` numbering with holes, from its first instance's rows.

import type { NumberingInput } from '@molstar/seq-core/model';
import type { CifFrame } from 'molstar/lib/mol-io/reader/cif';

export interface SiftsSegment {
  readonly entityId: string;
  readonly asymId: string;
  readonly accession: string;
  readonly best: boolean;
  /** label_seq_id range, inclusive. */
  readonly seqStart: number;
  readonly seqEnd: number;
  /** UniProt range, inclusive. */
  readonly unpStart: number;
  readonly unpEnd: number;
}

/** The file's SIFTS segments; [] when it has none. */
export function siftsSegments(frame: CifFrame | null): SiftsSegment[] {
  const cat = frame?.categories.pdbx_sifts_unp_segments;
  if (!cat) return [];
  const f = (n: string) => cat.getField(n);
  const [entity, asym, acc, best, s0, s1, u0, u1] = [
    'entity_id',
    'asym_id',
    'unp_acc',
    'best_mapping',
    'seq_id_start',
    'seq_id_end',
    'unp_start',
    'unp_end',
  ].map(f);
  if (!entity || !asym || !acc || !s0 || !s1 || !u0 || !u1) return [];
  const out: SiftsSegment[] = [];
  for (let i = 0; i < cat.rowCount; i++) {
    out.push({
      entityId: entity.str(i),
      asymId: asym.str(i),
      accession: acc.str(i),
      best: best ? best.str(i) === 'y' : true,
      seqStart: s0.int(i),
      seqEnd: s1.int(i),
      unpStart: u0.int(i),
      unpEnd: u1.int(i),
    });
  }
  return out;
}

/**
 * An entity's `uniprot:<acc>` numberings (`[index_start, label_start, length]` segments; index 0 is label_seq_id
 * `first`), from its first instance's rows: the best mapping's accessions, or every accession (`all`: isoforms too).
 */
export function siftsNumberings(
  frame: CifFrame | null,
  entityId: string,
  first: number,
  all: boolean,
): { numberings: NumberingInput[]; accessions: string[] } {
  const rows = siftsSegments(frame).filter((r) => r.entityId === entityId && (all || r.best));
  const asym = rows[0]?.asymId;
  const byAcc = new Map<string, [number, number, number][]>();
  for (const r of rows) {
    if (r.asymId !== asym) continue;
    const segs = byAcc.get(r.accession) ?? [];
    segs.push([r.seqStart - first, r.unpStart, r.seqEnd - r.seqStart + 1]);
    byAcc.set(r.accession, segs);
  }
  const numberings: NumberingInput[] = [...byAcc].map(([acc, segments]) => ({
    name: `uniprot:${acc}`,
    method: 'segments',
    segments: segments.sort((a, b) => a[0] - b[0]),
  }));
  return { numberings, accessions: [...byAcc.keys()] };
}
