/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { parseTrr } from '@molstar/io/reader/trr/parser';
import { coordinatesFromTrr } from '@molstar/model/formats/structure/trr';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { CoordinatesFormatCategory } from './category.js';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

export { CoordinatesFromTrr };
type CoordinatesFromTrr = typeof CoordinatesFromTrr;
const CoordinatesFromTrr = PluginStateTransform.BuiltIn({
  name: 'coordinates-from-trr',
  display: { name: 'Parse TRR', description: 'Parse TRR binary data.' },
  from: [SO.Data.Binary],
  to: SO.Molecule.Coordinates,
})({
  apply({ a }) {
    return Task.create('Parse TRR', async (ctx) => {
      const parsed = await parseTrr(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const coordinates = await coordinatesFromTrr(parsed.result).runInContext(ctx);
      return new SO.Molecule.Coordinates(coordinates, { label: a.label, description: 'Coordinates' });
    });
  },
});

export { TrrProvider };
const TrrProvider = DataFormatProvider({
  name: 'trr',
  label: 'TRR',
  description: 'TRR',
  category: CoordinatesFormatCategory,
  binaryExtensions: ['trr'],
  parse: (plugin, data) => {
    const coordinates = plugin.state.data.build().to(data).apply(CoordinatesFromTrr);

    return coordinates.commit();
  },
});

type TrrProvider = typeof TrrProvider;

/** The Trr data format. */
export const Trr: PluginRegistryEntry = {
  formats: [TrrProvider],
};
