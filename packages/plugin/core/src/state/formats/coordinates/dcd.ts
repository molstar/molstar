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
import { parseDcd } from '@molstar/io/reader/dcd/parser';
import { coordinatesFromDcd } from '@molstar/model/formats/structure/dcd';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { CoordinatesFormatCategory } from './category.js';

export { CoordinatesFromDcd };
type CoordinatesFromDcd = typeof CoordinatesFromDcd;
const CoordinatesFromDcd = PluginStateTransform.BuiltIn({
  name: 'coordinates-from-dcd',
  display: { name: 'Parse DCD', description: 'Parse DCD binary data.' },
  from: [SO.Data.Binary],
  to: SO.Molecule.Coordinates,
})({
  apply({ a }) {
    return Task.create('Parse DCD', async (ctx) => {
      const parsed = await parseDcd(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const coordinates = await coordinatesFromDcd(parsed.result).runInContext(ctx);
      return new SO.Molecule.Coordinates(coordinates, { label: a.label, description: 'Coordinates' });
    });
  },
});

export { DcdProvider };
const DcdProvider = DataFormatProvider({
  label: 'DCD',
  description: 'DCD',
  category: CoordinatesFormatCategory,
  binaryExtensions: ['dcd'],
  parse: (plugin, data) => {
    const coordinates = plugin.state.data.build().to(data).apply(CoordinatesFromDcd);

    return coordinates.commit();
  },
});

type DcdProvider = typeof DcdProvider;
