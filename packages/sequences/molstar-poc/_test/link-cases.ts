/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The case of the link's 3D-hover guard (link.heap.ts) and its quick twin (link.quick.ts): a Mol* hover, a loci of
// chain C of 1TUP made once (Mol* makes its hover loci; the link only reads them), to member runs and on to Seq* runs
// in reused `Extents`. The heap guard and the quick twin share it.

import { compile } from '@molstar/seq-core/format';
import { CoordinateSystemKind, Extents, KeyKind } from '@molstar/seq-core/model';
import { Translator } from '@molstar/seq-core/translate';
import { expect } from 'vitest';
import { readStructure, StructureLink, structureOf } from '../src/index.js';
import { bytesOf, loaders, p53Doc, registries } from './fixtures.js';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 50_000;

/** The step: a Mol* hover to member runs and Seq* runs; `sink` reads what it wrote. */
export async function molHover(): Promise<{ step: (i: number) => void; sink: () => number }> {
  const { set } = await compile(p53Doc([{ ref: '1tup', pdb: '1TUP', entity: '3' }]), registries(), loaders());
  const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
  const link = StructureLink.attach(new Translator(), '3d', structure, set, '1tup');
  const t = { kind: 0, id: 0, row: -1 };
  set.byName('1tup/C', t);
  expect(t.kind).toBe(KeyKind.unit);
  const at = (i: number) => {
    const e = Extents.create(1);
    Extents.add(e, { kind: CoordinateSystemKind.unit, id: t.id, row: -1 }, i, i + 1);
    return link.toLoci(e);
  };
  const hovers = Array.from({ length: 64 }, (_, i) => at(i * 3));
  const members = Extents.create(8);
  const seq = Extents.create(8);
  let sink = 0;
  const step = (i: number) => {
    members.count = 0;
    seq.count = 0;
    link.fromLoci(hovers[i & 63], members);
    link.toSeq(members, seq);
    sink = (sink + seq.start[0]) | 0;
  };
  return { step, sink: () => sink };
}
