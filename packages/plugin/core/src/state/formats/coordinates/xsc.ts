/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection, type StateTransform } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { parseXsc } from '@molstar/io/reader/xsc/parser';
import type { PluginContext } from '@molstar/plugin/context';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import {
  CoordinatesWithBoxes,
  TrajectoryFromModelAndCoordinates,
} from '@molstar/plugin/state/transforms/structure/hierarchy';

// A NAMD extended system file (.xsc, .xst) opened on its own (a dropped file): the periodic boxes of a DCD's frames,
// with their origins. It is parsed (`ParseXsc`) and given to the last trajectory built from
// coordinates: its coordinates get the boxes (`CoordinatesWithBoxes`), replacing boxes they had.

export { ParseXsc };
type ParseXsc = typeof ParseXsc;
const ParseXsc = PluginStateTransform.BuiltIn({
  name: 'parse-xsc',
  display: { name: 'Parse NAMD XSC', description: 'Parse a NAMD extended system file (.xsc, .xst): periodic boxes.' },
  from: SO.Data.String,
  to: SO.Format.Xsc,
})({
  apply({ a }) {
    return Task.create('Parse NAMD XSC', async (ctx) => {
      const parsed = await parseXsc(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      const n = parsed.result.entries.length;
      return new SO.Format.Xsc(parsed.result, { label: a.label, description: `${n} box${n === 1 ? '' : 'es'}` });
    });
  },
});

/** Files that give coordinates their periodic boxes; not coordinates (`load.trajectory` does not offer them). */
export const PeriodicBoxesFormatCategory = 'Periodic boxes';

/**
 * Gives the last trajectory built from coordinates the boxes in `boxesRef` (a `ParseXsc` cell): `CoordinatesWithBoxes` on its
 * coordinates, the trajectory pointed at it; boxes it had are replaced. Returns whether one was found.
 */
export async function attachPeriodicBoxes(plugin: PluginContext, boxesRef: StateTransform.Ref): Promise<boolean> {
  const state = plugin.state.data;
  const trajectories = state.select(StateSelection.Generators.ofTransformer(TrajectoryFromModelAndCoordinates));
  const trajectory = trajectories[trajectories.length - 1];
  if (!trajectory) {
    plugin.log.warn('NAMD box file: no trajectory built from coordinates to give its boxes to.');
    return false;
  }
  const params = trajectory.transform.params as { modelRef: string; coordinatesRef: string };
  let coordinatesRef = params.coordinatesRef;
  const current = state.cells.get(coordinatesRef);
  // boxes it had: from the coordinates under them, the old ones deleted
  const replaced = current?.transform.transformer === CoordinatesWithBoxes ? coordinatesRef : undefined;
  if (replaced) coordinatesRef = state.tree.transforms.get(replaced)!.parent;
  const b = state.build();
  const withBoxes = b.to(coordinatesRef).apply(CoordinatesWithBoxes, { boxesRef }, { dependsOn: [boxesRef] });
  b.to(trajectory.transform.ref).update(TrajectoryFromModelAndCoordinates, (p) => ({
    ...p,
    coordinatesRef: withBoxes.ref,
  }));
  if (replaced) b.delete(replaced);
  await b.commit();
  return true;
}

export { XscProvider };
const XscProvider = DataFormatProvider({
  name: 'xsc',
  label: 'NAMD XSC',
  description: 'NAMD extended system file (periodic boxes)',
  category: PeriodicBoxesFormatCategory,
  stringExtensions: ['xsc', 'xst'],
  parse: async (plugin, data) => ({ boxes: await plugin.state.data.build().to(data).apply(ParseXsc).commit() }),
  visuals: async (plugin, parsed: { boxes: { ref: string } }) => {
    await attachPeriodicBoxes(plugin, parsed.boxes.ref as StateTransform.Ref);
  },
});

type XscProvider = typeof XscProvider;

/** The NAMD XSC data format. */
export const Xsc: PluginRegistryEntry = {
  formats: [XscProvider],
};
