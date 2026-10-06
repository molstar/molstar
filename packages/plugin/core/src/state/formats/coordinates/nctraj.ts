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
import { parseNctraj } from '@molstar/io/reader/nctraj/parser';
import { coordinatesFromNctraj } from '@molstar/model/formats/structure/nctraj';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { CoordinatesFormatCategory } from './category.js';

export { CoordinatesFromNctraj };
type CoordinatesFromNctraj = typeof CoordinatesFromNctraj;
const CoordinatesFromNctraj = PluginStateTransform.BuiltIn({
  name: 'coordinates-from-nctraj',
  display: { name: 'Parse NCTRAJ', description: 'Parse NCTRAJ binary data.' },
  from: [SO.Data.Binary],
  to: SO.Molecule.Coordinates,
})({
  apply({ a }) {
    return Task.create('Parse NCTRAJ', async (ctx) => {
      const parsed = await parseNctraj(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const coordinates = await coordinatesFromNctraj(parsed.result).runInContext(ctx);
      return new SO.Molecule.Coordinates(coordinates, { label: a.label, description: 'Coordinates' });
    });
  },
});

export { NctrajProvider };
const NctrajProvider = DataFormatProvider({
  name: 'nctraj',
  label: 'NCTRAJ',
  description: 'NCTRAJ',
  category: CoordinatesFormatCategory,
  binaryExtensions: ['nc', 'nctraj'],
  parse: (plugin, data) => {
    const coordinates = plugin.state.data.build().to(data).apply(CoordinatesFromNctraj);

    return coordinates.commit();
  },
});

type NctrajProvider = typeof NctrajProvider;
