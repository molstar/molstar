/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The system registry with no viewer (seqspec:libraries#the-system-registry, seqspec:libraries#coordinates): a set's
// own systems by key, a host's system (a structure's chain, here a stand-in) registered beside them, and targets mapped
// between the host's system, the set's units and an alignment's columns, both ways; changes reported through the
// callback.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { type CoordinateSystem, CoordinateSystemKind, Extents, KeyKind, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { compile } from '../../src/format/index.js';
import { Coordinates, type SystemDescriptor, SystemRegistry, TranslateMaps } from '../../src/translate/index.js';
import { registries } from '../format/load-fixtures.js';

const K = CoordinateSystemKind;

/** P (10 residues) as unit `u` of complex `c`; an alignment whose row 0 is P with a gap at column 3. */
async function setOf(): Promise<SeqSet> {
  const b = createSeqViewSpec({ title: 'Registry' });
  b.data.sequence({ ref: 'p', value: 'MKVLATGHEK' });
  b.data.complex({ ref: 'c' }).unit({ ref: 'u', sequence: 'p' });
  b.data.alignment({ ref: 'msa', alphabet: 'protein', value: ['MKV-LATGHEK', 'MKVQLAT-HEK'] });
  b.data.mapping({ ref: 'm0', from: 'msa/0', to: 'p', method: 'identity' } as never);
  b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler();
  const { set, issues } = await compile(b.getDocument(), registries(false), {
    load: () => Promise.reject(new Error('no files')),
  });
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  return set;
}

const chain = (key: string, stop: number, more: Partial<SystemDescriptor> = {}): SystemDescriptor => ({
  key,
  kind: 'external',
  start: 0,
  stop,
  refs: {},
  members: null,
  structure: null,
  numberings: [],
  owner: 'test',
  ...more,
});
function unitOf(set: SeqSet, ref: string): CoordinateSystem {
  const t = { kind: 0, id: 0, row: -1 };
  if (!set.byName(ref, t) || t.kind !== KeyKind.unit) throw new Error(`no unit ${ref}`);
  return { kind: K.unit, id: t.id, row: -1 };
}
/** Coordinates as [key, start, stop] triples. */
const named = (c: Coordinates) =>
  Array.from({ length: c.runs.length / 3 }, (_, i) => [c.systems[c.runs[3 * i]], c.runs[3 * i + 1], c.runs[3 * i + 2]]);
const coords = (key: string, start: number, stop: number): Coordinates => ({
  systems: [key],
  runs: Int32Array.of(0, start, stop),
});

/** The registry over the set, a chain `3d/A` over its residues 1..10 mapped onto unit `u` (residue i → index i - 1). */
async function registry() {
  const set = await setOf();
  let changes = 0;
  const r = new SystemRegistry({ onChange: () => changes++ });
  r.update(set);
  const u = unitOf(set, 'u');
  r.register(chain('3d/A', 11), (s) => [TranslateMaps.bySegments(s, u, Int32Array.of(1, 0, 10))]);
  r.register(chain('3d', 0, { members: ['3d/A'] }));
  return { set, r, changes: () => changes };
}

describe('SystemRegistry, headless', () => {
  it("names the set's own systems by key, a host's beside them, and reports changes", async () => {
    const { r, changes } = await registry();
    expect(changes()).toBeGreaterThan(0);
    expect(r.descriptor('seqstar:u')?.kind).toBe('unit');
    expect(r.descriptor('seqstar:msa')?.kind).toBe('columns');
    expect(r.descriptor('seqstar:msa/1')?.kind).toBe('row');
    expect(r.descriptor('seqstar:c')?.members).toEqual(['seqstar:u']);
    expect(r.descriptor('3d/A')?.owner).toBe('test');
    expect(r.system('seqstar:nothing')).toBeNull();
    const before = changes();
    r.unregister('3d/A');
    expect(changes()).toBe(before + 1);
  });

  it("maps a structure's residues onto the alignment's columns, and columns back onto the structure", async () => {
    const { r } = await registry();
    // Residues 4 and 5 (L, A): P's indices 3 and 4, row 0's columns 4 and 5 (the gap at column 3 between).
    expect(named(r.map(coords('3d/A', 4, 6), 'seqstar:msa'))).toEqual([['seqstar:msa', 4, 6]]);
    // Columns 0..5 hold P's 0..4 (column 3 a gap): residues 1..5.
    expect(named(r.map(coords('seqstar:msa', 0, 6), '3d/A'))).toEqual([['3d/A', 1, 6]]);
    // Onto a composite system: through each member.
    expect(named(r.map(coords('seqstar:msa', 4, 5), '3d'))).toEqual([['3d/A', 4, 5]]);
    // The other row: column 4 is its L too, on no structure (row 1 is mapped onto nothing).
    expect(named(r.map(coords('seqstar:msa', 4, 5), 'seqstar:msa/1'))).toEqual([['seqstar:msa/1', 4, 5]]);
  });

  it('maps a unit onto a row, and loci of the set as targets', async () => {
    const { set, r } = await registry();
    // A row's own system counts its residues: P's 2..5 are row 0's.
    expect(named(r.map(coords('seqstar:u', 2, 5), 'seqstar:msa/0'))).toEqual([['seqstar:msa/0', 2, 5]]);
    const [loci] = r.lociOf(coords('seqstar:msa', 0, 2));
    expect(loci.kind).toBe('column-loci');
    expect(named(r.coordinatesOf(loci))).toEqual([['seqstar:msa', 0, 2]]);
    expect(named(r.map(loci, '3d/A'))).toEqual([['3d/A', 1, 3]]);
    expect(Coordinates.toSeqLoci(set, Extents.create(1))).toBeNull();
    // A unit's loci list its complex too (its descriptor's members name the units).
    const [unitLoci] = r.lociOf(coords('3d/A', 1, 3));
    expect(r.coordinatesOf(unitLoci).systems).toEqual(['seqstar:u', 'seqstar:c']);
  });

  it('gives nothing for an unknown target or a set not given', async () => {
    const { r } = await registry();
    expect(named(r.map(coords('3d/A', 1, 3), 'seqstar:nothing'))).toEqual([]);
    expect(named(r.map(coords('nothing', 1, 3), 'seqstar:u'))).toEqual([]);
    expect(named(new SystemRegistry().map(coords('3d/A', 1, 3), 'seqstar:u'))).toEqual([]);
  });
});
