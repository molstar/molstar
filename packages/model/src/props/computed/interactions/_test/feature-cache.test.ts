/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
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
