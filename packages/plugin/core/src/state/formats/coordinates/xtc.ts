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
import { parseXtc } from '@molstar/io/reader/xtc/parser';
import { coordinatesFromXtc } from '@molstar/model/formats/structure/xtc';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { CoordinatesFormatCategory } from './category.js';

export { CoordinatesFromXtc };
type CoordinatesFromXtc = typeof CoordinatesFromXtc;
const CoordinatesFromXtc = PluginStateTransform.BuiltIn({
  name: 'coordinates-from-xtc',
  display: { name: 'Parse XTC', description: 'Parse XTC binary data.' },
  from: [SO.Data.Binary],
  to: SO.Molecule.Coordinates,
})({
  apply({ a }) {
    return Task.create('Parse XTC', async (ctx) => {
      const parsed = await parseXtc(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const coordinates = await coordinatesFromXtc(parsed.result).runInContext(ctx);
      return new SO.Molecule.Coordinates(coordinates, { label: a.label, description: 'Coordinates' });
    });
  },
});

export { XtcProvider };
const XtcProvider = DataFormatProvider({
  name: 'xtc',
  label: 'XTC',
  description: 'XTC',
  category: CoordinatesFormatCategory,
  binaryExtensions: ['xtc'],
  parse: (plugin, data) => {
    const coordinates = plugin.state.data.build().to(data).apply(CoordinatesFromXtc);

    return coordinates.commit();
  },
});

type XtcProvider = typeof XtcProvider;
