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
import { parsePrmtop } from '@molstar/io/reader/prmtop/parser';
import { topologyFromPrmtop } from '@molstar/model/formats/structure/prmtop';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { TopologyFormatCategory } from './category.js';

export { ParsePrmtop };
type ParsePrmtop = typeof ParsePrmtop;
const ParsePrmtop = PluginStateTransform.BuiltIn({
  name: 'parse-prmtop',
  display: { name: 'Parse PRMTOP', description: 'Parse PRMTOP from String data' },
  from: [SO.Data.String],
  to: SO.Format.Prmtop,
})({
  apply({ a }) {
    return Task.create('Parse PRMTOP', async (ctx) => {
      const parsed = await parsePrmtop(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Prmtop(parsed.result);
    });
  },
});

export { TopologyFromPrmtop };
type TopologyFromPrmtop = typeof TopologyFromPrmtop;
const TopologyFromPrmtop = PluginStateTransform.BuiltIn({
  name: 'topology-from-prmtop',
  display: { name: 'PRMTOP Topology', description: 'Create topology from PRMTOP.' },
  from: [SO.Format.Prmtop],
  to: SO.Molecule.Topology,
})({
  apply({ a }) {
    return Task.create('Create Topology', async (ctx) => {
      const topology = await topologyFromPrmtop(a.data).runInContext(ctx);
      return new SO.Molecule.Topology(topology, { label: topology.label || a.label, description: 'Topology' });
    });
  },
});

export { PrmtopProvider };
const PrmtopProvider = DataFormatProvider({
  label: 'PRMTOP',
  description: 'PRMTOP',
  category: TopologyFormatCategory,
  stringExtensions: ['prmtop', 'parm7'],
  parse: async (plugin, data) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParsePrmtop, {}, { state: { isGhost: true } });
    const topology = format.apply(TopologyFromPrmtop);

    await format.commit();

    return { format: format.selector, topology: topology.selector };
  },
});

type PrmtopProvider = typeof PrmtopProvider;
