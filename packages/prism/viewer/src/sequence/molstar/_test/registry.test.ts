/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The link on core's system registry, no viewer (seqspec:libraries#the-system-registry,
// seqspec:libraries#the-mol-adapter): 1TUP's chains registered beside P04637 and an alignment whose first row is
// P04637, and targets mapped between the structure, the sequence set and the alignment's columns both ways.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { compile } from '@molstar/seq-core/format';
import { Extents } from '@molstar/seq-core/model';
import { type Coordinates, SystemRegistry } from '@molstar/seq-core/translate';
import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { readStructure, StructureLink, structureOf } from '../index.js';
import { bytesOf, loaders, P53, registries } from './fixtures.js';

/** P04637 with 1TUP's p53 chains mapped onto it by SIFTS, and an alignment: P04637 with two gap columns after 100. */
async function linked() {
  const b = createSeqViewSpec({ title: 'p53' });
  const p53 = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
  b.data.complex({
    ref: '1tup',
    provider: 'molstar:structure',
    source: { service: 'pdbe', resource: 'structure', id: '1TUP' },
  });
  b.data.mapping({
    ref: '1tup-map',
    from: p53,
    to: '1tup/entity:3',
    method: 'by_numbering',
    numbering: 'uniprot:P04637',
  } as never);
  const row = `${P53.sequence.slice(0, 100)}--${P53.sequence.slice(100)}`;
  b.data.alignment({ ref: 'msa', alphabet: 'protein', value: [row, row.replace(/-/g, 'A')] });
  b.data.mapping({ ref: 'msa-map', from: 'msa/0', to: 'p53', method: 'identity' } as never);
  b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler();
  const { set, issues } = await compile(b.getDocument(), registries(), loaders());
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
  let changes = 0;
  const r = new SystemRegistry({ onChange: () => changes++ });
  r.update(set);
  const link = StructureLink.attach(r, '3d-1', structure, set, '1tup');
  return { r, link, changes: () => changes };
}

const coords = (key: string, start: number, stop: number): Coordinates => ({
  systems: [key],
  runs: Int32Array.of(0, start, stop),
});
const named = (c: Coordinates) =>
  Array.from({ length: c.runs.length / 3 }, (_, i) => [c.systems[c.runs[3 * i]], c.runs[3 * i + 1], c.runs[3 * i + 2]]);

describe('StructureLink on a SystemRegistry, no viewer', () => {
  it("maps a structure's residue onto the alignment's columns, and a column back onto every copy", async () => {
    const { r, link, changes } = await linked();
    expect(changes()).toBeGreaterThan(0);
    // Chain C's R175 (its unit's index 81) is P04637's 174, row 0's column 176 (past the two gap columns).
    const c = link.chains.find((x) => x.chain.labelAsymId === 'C');
    expect(c?.key).toBe('molstar:3d-1/C@1_555');
    expect(named(r.map(coords('seqstar:1tup/C', 81, 82), 'seqstar:p53'))).toEqual([['seqstar:p53', 174, 175]]);
    expect(named(r.map(coords('seqstar:1tup/C', 81, 82), 'seqstar:msa'))).toEqual([['seqstar:msa', 176, 177]]);
    // Back as Mol*'s loci: the registry's runs, then the link's.
    const ext = Extents.create(1);
    r.extentsOf(r.map(coords('seqstar:msa', 176, 177), 'seqstar:1tup/C'), ext);
    expect(StructureElement.Loci.size(link.toLoci(ext))).toBeGreaterThan(0);
    // Column 176 back onto the structure: R175 of chains C, D and E (the entity's copies).
    const back = named(r.map(coords('seqstar:msa', 176, 177), link.key));
    expect(back.map((x) => x[0])).toEqual(['molstar:3d-1/C@1_555', 'molstar:3d-1/D@1_555', 'molstar:3d-1/E@1_555']);
    expect(back.every((x) => x[2] === (x[1] as number) + 1)).toBe(true);
    // A gap column is no residue.
    expect(named(r.map(coords('seqstar:msa', 100, 102), link.key))).toEqual([]);
  });

  it('detach unregisters the structure: its keys map to nothing, the set keeps its own', async () => {
    const { r, link, changes } = await linked();
    const before = changes();
    link.detach();
    expect(changes()).toBeGreaterThan(before);
    expect(named(r.map(coords('seqstar:msa', 176, 177), link.key))).toEqual([]);
    expect(named(r.map(coords('seqstar:1tup/C', 81, 82), 'seqstar:msa'))).toEqual([['seqstar:msa', 176, 177]]);
  });
});
