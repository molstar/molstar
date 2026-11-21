/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqSpecBuilder } from '@molstar/seqx-core/lang';

// PDBe's SIFTS as an application-layer runtime (examples.md, Showcases): the mappings response of an entry
// (`/pdbe/api/mappings/uniprot/<entry>`) read into chain mappings, and those into `structure` maps (maps.md: segments with
// a `provenance` word, onto the entity of the structure a source describes) the way a host with its own SIFTS would write
// them. The response is read as published: `unp_start` and `unp_end` are UniProt's numbers (from 1, the last included),
// `start.residue_number` and `end.residue_number` the entity's positions (from 1, `label_seq_id`), and the author's
// numbers beside them (`author_residue_number`, null where SIFTS has none).

/** One mapping of the response: a chain's segment of the entity onto a UniProt entry. */
export interface SiftsMapping {
  readonly entry: string;
  readonly accession: string;
  /** The entity (`entity.id`) and the chain instance (`struct_asym_id`) and the author's chain. */
  readonly entity: string;
  readonly asym: string;
  readonly chain: string;
  /** UniProt's numbers, from 1, both included. */
  readonly unpStart: number;
  readonly unpEnd: number;
  /** The entity's positions, from 1, both included. */
  readonly start: number;
  readonly end: number;
  /** The author's numbers of the segment's ends (null: SIFTS names none). */
  readonly authStart: number | null;
  readonly authEnd: number | null;
  readonly identity: number;
  readonly coverage: number;
}

interface Json {
  [key: string]: any;
}

/** Reads a mappings response: every UniProt entry's mappings, in the order the response lists them. */
export function parseSifts(text: string): SiftsMapping[] {
  const root = JSON.parse(text) as Json;
  const out: SiftsMapping[] = [];
  for (const entry of Object.keys(root)) {
    const uniprot = root[entry]?.UniProt as Json | undefined;
    for (const accession of Object.keys(uniprot ?? {})) {
      for (const m of (uniprot![accession].mappings ?? []) as Json[]) {
        out.push({
          entry: entry.toUpperCase(),
          accession,
          entity: String(m.entity_id),
          asym: String(m.struct_asym_id),
          chain: String(m.chain_id),
          unpStart: m.unp_start,
          unpEnd: m.unp_end,
          start: m.start.residue_number,
          end: m.end.residue_number,
          authStart: m.start.author_residue_number ?? null,
          authEnd: m.end.author_residue_number ?? null,
          identity: m.identity,
          coverage: m.coverage,
        });
      }
    }
  }
  return out;
}

/** The UniProt entry's chains: what a coverage track lists (`1FIN A, C`). */
export function chainsOf(mappings: readonly SiftsMapping[], accession: string): string[] {
  return [...new Set(mappings.filter((m) => m.accession === accession).map((m) => m.chain))];
}

/**
 * Whether a mapping pairs its two runs one to one: the entity's residues (`end − start + 1`) as many as UniProt's
 * (`unpEnd − unpStart + 1`). SIFTS also writes segments with an insertion or a deletion inside; `from + i ↦ to + i` cannot
 * say those, so they are left out of the maps and said (`siftsIssues`), never mapped wrong.
 */
export const isOneToOne = (m: SiftsMapping) => m.end - m.start === m.unpEnd - m.unpStart && m.end >= m.start;

/** What the mappings left out, one sentence each: the application shows them. */
export function siftsIssues(mappings: readonly SiftsMapping[]): string[] {
  return mappings
    .filter((m) => !isOneToOne(m))
    .map(
      (m) =>
        `${m.entry} chain ${m.chain}: the entity's ${m.start}–${m.end} (${m.end - m.start + 1}) and ${m.accession} ${m.unpStart}–${m.unpEnd} (${m.unpEnd - m.unpStart + 1}) differ in length: not mapped`,
    );
}

/** One `structure` map: a chain's segment onto an accession (SIFTS per chain). */
export interface SiftsSegment {
  readonly accession: string;
  /** The chain instance (`struct_asym_id`: a unit's id in the described complex) and the author's chain (`chain_id`). */
  readonly asym: string;
  readonly chain: string;
  /** UniProt's positions, 0-based, and the unit's, 0-based (the entity's positions: `label_seq_id` − 1). */
  readonly from: number;
  readonly to: number;
  readonly length: number;
}

/**
 * The segments of the mappings, one per one-to-one mapping (a chain's own), `from + i ↦ to + i`. Which residues of the
 * unit were observed is the structure source's (its unit's `observed` runs), never SIFTS's: a segment covers the entity.
 */
export function siftsSegments(mappings: readonly SiftsMapping[]): SiftsSegment[] {
  return mappings.filter(isOneToOne).map((m) => ({
    accession: m.accession,
    asym: m.asym,
    chain: m.chain,
    from: m.unpStart - 1,
    to: m.start - 1,
    length: m.end - m.start + 1,
  }));
}

/** The id of a segment's map: `sifts-<entry>-<accession>-<asym>` (a second segment of the pair is numbered). */
export const siftsMapId = (entry: string, segment: SiftsSegment, n = 0) =>
  `sifts-${entry.toLowerCase()}-${segment.accession.toLowerCase()}-${segment.asym.toLowerCase()}${n ? `-${n + 1}` : ''}`;

/**
 * Writes the `structure` maps of the mappings to `b`: one per chain's segment, from the sequence `proteins` names for its
 * accession onto the unit `unitAddress` gives for the chain instance (`struct_asym_id`, what the structure source calls
 * its units); `provenance` `sifts`, ranked first (`priority` 1) so a section's paths take it before an alignment's. A
 * mapping of an accession no sequence is given for is left out. `suffix` tells a second set of maps (onto the spec's own
 * sequences of the units) from the first. Returns the maps' addresses by id.
 */
export function addSiftsMaps(
  b: SeqSpecBuilder,
  entry: string,
  mappings: readonly SiftsMapping[],
  proteins: Readonly<Record<string, string>>,
  unitAddress: (asym: string) => string,
  suffix = '',
): Map<string, string> {
  const out = new Map<string, string>();
  const seen = new Map<string, number>();
  for (const s of siftsSegments(mappings)) {
    const from = proteins[s.accession];
    if (from === undefined) continue;
    const pair = `${s.accession}|${s.asym}`;
    const n = seen.get(pair) ?? 0;
    seen.set(pair, n + 1);
    const id = `${siftsMapId(entry, s, n)}${suffix}`;
    out.set(
      id,
      b.data.map(id, {
        from: from as never,
        to: unitAddress(s.asym) as never,
        method: 'structure',
        segments: { from: [s.from], to: [s.to], length: [s.length], provenance: ['sifts'] },
        priority: 1,
      }),
    );
  }
  return out;
}
