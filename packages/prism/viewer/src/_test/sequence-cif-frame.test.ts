/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import { headlessViewer, until } from './harness.js';

// The followed structure's document reads the categories of the CIF block it came from (sequence.md §2): SIFTS
// numberings, a predicted model's origin.

const fixtures = path.resolve(import.meta.dirname, '../sequence/molstar/_test/fixtures/structures');
const fixture = (name: string) => new File([fs.readFileSync(path.join(fixtures, name))], name);

async function assembled(name: string) {
  const v = await headlessViewer({ files: [fixture(name)] });
  const viewer = (await v.views.sequence!.load()).viewers.get();
  await until(() => !!viewer.shown.value && !!viewer.seq.graph.assembled.value);
  return viewer.seq.graph.assembled.value!.set;
}

describe("the followed structure's CIF block", () => {
  it('gives 1TUP its SIFTS numbering onto P04637', async () => {
    const set = await assembled('1tup.bcif');
    const uniprot = set.strings.find('uniprot:P04637');
    expect(uniprot).toBeGreaterThanOrEqual(0);
    expect(set.molecules.some((m) => !!m?.numberings.byName(uniprot))).toBe(true);
  });

  it('marks the AlphaFold model predicted', async () => {
    const set = await assembled('AF-P04637-F1-model_v6.bcif');
    const origins = new Set<string>();
    for (let u = 0; u < set.units.count; u++) origins.add(set.strings.get(set.units.origin[u]) ?? '');
    expect([...origins]).toEqual(['predicted']);
  });
});
