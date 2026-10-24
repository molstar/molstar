/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import { computeInteractions } from '../interactions.js';
import { HydrophobicParams, InteractionsParams, type InteractionsProps } from '../params.js';

const ctx = { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() };

async function load1cbs() {
  const file = path.resolve(import.meta.dirname, '../../../../../../../data/examples/1cbs_updated.cif');
  const parsed = await CIF.parse(fs.readFileSync(file, 'utf8')).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  return Structure.ofModel(trajectory.representative);
}

function featureKeys(structure: Structure) {
  return structure.units.map((u) => [...u.transientCache.keys()].filter((k) => `${k}`.startsWith('features-')).length);
}

// The features depend on the feature providers only: a change of a contact's params (a slider) reuses them, and the
// units' caches do not grow by an entry per change (perf-audit-2 06/P2).
describe('interaction features cache', () => {
  it('is keyed by the feature providers, not by every param', async () => {
    const structure = await load1cbs();
    const defaults = PD.getDefaultValues(InteractionsParams);
    const a = await computeInteractions(ctx, structure, defaults);
    const keys = featureKeys(structure);

    const hb = defaults.providers['hydrogen-bonds'];
    if (hb.name !== 'on') throw new Error('expected H-bonds on by default');
    const moved: InteractionsProps = {
      ...defaults,
      providers: {
        ...defaults.providers,
        'hydrogen-bonds': { name: 'on', params: { ...hb.params, distanceMax: 3.2 } },
      },
    };
    const b = await computeInteractions(ctx, structure, moved);
    for (const u of structure.units) expect(b.unitsFeatures.get(u.id)).toBe(a.unitsFeatures.get(u.id));
    expect(featureKeys(structure)).toEqual(keys);

    // other providers need other features: hydrophobic contacts alone
    const providers = Object.fromEntries(
      Object.keys(defaults.providers).map((k) => [k, { name: 'off', params: {} }]),
    ) as InteractionsProps['providers'];
    providers.hydrophobic = { name: 'on', params: PD.getDefaultValues(HydrophobicParams) };
    const bridges = Object.fromEntries(
      Object.keys(defaults.bridges).map((k) => [k, { name: 'off', params: {} }]),
    ) as InteractionsProps['bridges'];
    const c = await computeInteractions(ctx, structure, { ...defaults, providers, bridges });
    const unit = structure.units[0].id;
    expect(c.unitsFeatures.get(unit)).not.toBe(a.unitsFeatures.get(unit));
    expect(c.unitsFeatures.get(unit)!.count).toBeGreaterThan(0);
  });
});

// A tester's subset (index list and grid) is made once per type set for the features' lifetime, and the grids get an
// exact box without the 98-direction boundary (perf-audit-2 06/P5).
describe('interaction feature grids', () => {
  it('reuse a subset per type set; the full grid has the exact box and a sphere holding every feature', async () => {
    const structure = await load1cbs();
    const r = await computeInteractions(ctx, structure, PD.getDefaultValues(InteractionsParams));
    const unit = structure.units.reduce((a, b) => (a.elements.length >= b.elements.length ? a : b));
    const features = r.unitsFeatures.get(unit.id);
    expect(features.count).toBeGreaterThan(100);

    const types = new Set(Array.from(features.types).slice(0, 3));
    const subset = features.subset(types);
    expect(features.subset(new Set(Array.from(types).reverse()))).toBe(subset);
    expect(features.subset(new Set([...types, 255 as any]))).not.toBe(subset);

    const { box, sphere } = features.lookup3d.boundary;
    const { x, y, z, count } = features;
    const min = [Infinity, Infinity, Infinity],
      max = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < count; i++) {
      const p = [x[i], y[i], z[i]];
      for (let k = 0; k < 3; k++) {
        min[k] = Math.min(min[k], p[k]);
        max[k] = Math.max(max[k], p[k]);
      }
      const d = Math.hypot(p[0] - sphere.center[0], p[1] - sphere.center[1], p[2] - sphere.center[2]);
      expect(d).toBeLessThanOrEqual(sphere.radius * (1 + 1e-12));
    }
    expect(Array.from(box.min)).toEqual(min);
    expect(Array.from(box.max)).toEqual(max);

    // the subset grid finds what a scan finds
    const i0 = OrderedSet.getAt(subset.indices, 0);
    const found = subset.lookup3d.find(x[i0], y[i0], z[i0], 6);
    const hits = new Set<number>();
    for (let k = 0; k < found.count; k++) hits.add(OrderedSet.getAt(subset.indices, found.indices[k]));
    const scan = new Set<number>();
    for (let t = 0; t < OrderedSet.size(subset.indices); t++) {
      const j = OrderedSet.getAt(subset.indices, t);
      if (Math.hypot(x[j] - x[i0], y[j] - y[i0], z[j] - z[i0]) <= 6) scan.add(j);
    }
    expect(hits).toEqual(scan);
  });
});
