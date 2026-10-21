/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import { parseArtiatomiEm } from '@molstar/io/reader/artiatomi/em';
import { CIF } from '@molstar/io/reader/cif';
import { parseCryoEtDataPortalNdjson } from '@molstar/io/reader/cryoet/ndjson';
import { parseDynamoTbl } from '@molstar/io/reader/dynamo/tbl';
import { getRelionStarTomogramNames, parseRelionStar } from '@molstar/io/reader/relion/star';
import { parseSimularium } from '@molstar/io/reader/simularium/parser';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { createParticleListFromArtiatomiEm, getArtiatomiMotivelistTomogramIds } from '../em.js';
import { createParticleListFromCryoEtDataPortalNdjson } from '../ndjson.js';
import { createParticleListFromSimularium, getSimulariumFrameCount } from '../simularium.js';
import { createParticleListFromRelionStar } from '../star.js';
import { createParticleListFromDynamoTbl, getDynamoTblTomogramIds } from '../tbl.js';

// The fixtures and their sources: data/examples/particles/README.md.
const dataDir = path.resolve(import.meta.dirname, '../../../../../../data/examples/particles');

function readText(name: string) {
  return fs.readFileSync(path.join(dataDir, name), 'utf8');
}

function readBinary(name: string) {
  return new Uint8Array(fs.readFileSync(path.join(dataDir, name)));
}

function unwrap<T>(r: { isError: true; message: string } | { isError: false; result: T }): T {
  if (r.isError) throw new Error(r.message);
  return r.result;
}

function expectFinite(list: ParticleList) {
  for (let i = 0; i < list.coordinates.length; ++i) expect(Number.isFinite(list.coordinates[i])).toBe(true);
}

function expectUnitRotations(list: ParticleList) {
  const q = list.rotations!;
  expect(q.length).toBe(list.count * 4);
  for (let i = 0; i < list.count; ++i) {
    const o = i * 4;
    expect(Math.hypot(q[o], q[o + 1], q[o + 2], q[o + 3])).toBeCloseTo(1, 4);
  }
}

describe('particle fixtures in data/examples/particles', () => {
  it('reads the CryoET Data Portal oriented points', async () => {
    const data = unwrap(await parseCryoEtDataPortalNdjson(readText('cryoet-ribosomes.ndjson')).run());
    expect(data.records.length).toBe(200);
    expect(data.records.every((r) => r.type === 'orientedPoint')).toBe(true);

    const list = createParticleListFromCryoEtDataPortalNdjson(data, { pixelSize: 6.802 });
    expect(list.count).toBe(200);
    expectFinite(list);
    expectUnitRotations(list);
  });

  it('reads the RELION STAR file', async () => {
    const cif = unwrap(await CIF.parseText(readText('relion-subtomo.star')).run());
    const data = unwrap(parseRelionStar(cif));
    expect(data.opticsBlock).toBeDefined();
    expect(getRelionStarTomogramNames(cif)).toEqual([]);

    const list = createParticleListFromRelionStar(data);
    expect(list.count).toBe(300);
    expectFinite(list);
    expectUnitRotations(list);
    // pixel coordinates scaled by the optics block's 10 Å/px
    expect(list.coordinates[0]).toBeCloseTo(917.987, 2);
  });

  it('reads the Dynamo table', async () => {
    const data = unwrap(await parseDynamoTbl(readText('dynamo.tbl')).run());
    expect(data.rowCount).toBe(300);
    expect(getDynamoTblTomogramIds(data)).toEqual([1, 2]);

    const list = createParticleListFromDynamoTbl(data);
    expect(list.count).toBe(300);
    expectFinite(list);
    expectUnitRotations(list);

    expect(createParticleListFromDynamoTbl(data, { tomos: [2] }).count).toBe(150);
  });

  it('reads the Artiatomi motivelist', async () => {
    const data = unwrap(await parseArtiatomiEm(readBinary('artiatomi.em')).run());
    expect(data.header.dimX).toBe(20);
    expect(data.header.dimY).toBe(225);
    expect(getArtiatomiMotivelistTomogramIds(data).length).toBeGreaterThan(0);

    const list = createParticleListFromArtiatomiEm(data, { pixelSize: 13 });
    expect(list.count).toBe(225);
    expectFinite(list);
    expectUnitRotations(list);
  });

  it('reads the Simularium trajectory', async () => {
    const data = unwrap(await parseSimularium(readBinary('small.simularium')).run());
    expect(getSimulariumFrameCount(data)).toBe(5);

    for (let frameIndex = 0; frameIndex < 5; ++frameIndex) {
      const list = createParticleListFromSimularium(data, { frameIndex });
      expect(list.count).toBe(5);
      expectFinite(list);
    }
  });
});
