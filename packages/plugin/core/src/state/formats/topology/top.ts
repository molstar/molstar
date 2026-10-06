/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Neli Fonseca <neli@ebi.ac.uk>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { parseTop } from '@molstar/io/reader/top/parser';
import { topologyFromTop } from '@molstar/model/formats/structure/top';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { TopologyFormatCategory } from './category.js';

export { ParseTop };
type ParseTop = typeof ParseTop;
const ParseTop = PluginStateTransform.BuiltIn({
  name: 'parse-top',
  display: { name: 'Parse TOP', description: 'Parse TOP from String data' },
  from: [SO.Data.String],
  to: SO.Format.Top,
})({
  apply({ a }) {
    return Task.create('Parse TOP', async (ctx) => {
      const parsed = await parseTop(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Top(parsed.result);
    });
  },
});

export { TopologyFromTop };
type TopologyFromTop = typeof TopologyFromTop;
const TopologyFromTop = PluginStateTransform.BuiltIn({
  name: 'topology-from-top',
  display: { name: 'TOP Topology', description: 'Create topology from TOP.' },
  from: [SO.Format.Top],
  to: SO.Molecule.Topology,
})({
  apply({ a }) {
    return Task.create('Create Topology', async (ctx) => {
      const topology = await topologyFromTop(a.data).runInContext(ctx);
      return new SO.Molecule.Topology(topology, { label: topology.label || a.label, description: 'Topology' });
    });
  },
});

export { TopProvider };
const TopProvider = DataFormatProvider({
  name: 'top',
  label: 'TOP',
  description: 'TOP',
  category: TopologyFormatCategory,
  stringExtensions: ['top'],
  parse: async (plugin, data) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParseTop, {}, { state: { isGhost: true } });
    const topology = format.apply(TopologyFromTop);

    await format.commit();

    return { format: format.selector, topology: topology.selector };
  },
});

type TopProvider = typeof TopProvider;
