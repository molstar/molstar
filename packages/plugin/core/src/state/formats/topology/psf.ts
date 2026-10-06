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
import { parsePsf } from '@molstar/io/reader/psf/parser';
import { topologyFromPsf } from '@molstar/model/formats/structure/psf';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { TopologyFormatCategory } from './category.js';

export { ParsePsf };
type ParsePsf = typeof ParsePsf;
const ParsePsf = PluginStateTransform.BuiltIn({
  name: 'parse-psf',
  display: { name: 'Parse PSF', description: 'Parse PSF from String data' },
  from: [SO.Data.String],
  to: SO.Format.Psf,
})({
  apply({ a }) {
    return Task.create('Parse PSF', async (ctx) => {
      const parsed = await parsePsf(a.data).runInContext(ctx);
      if (parsed.isError) throw new Error(parsed.message);
      return new SO.Format.Psf(parsed.result);
    });
  },
});

export { TopologyFromPsf };
type TopologyFromPsf = typeof TopologyFromPsf;
const TopologyFromPsf = PluginStateTransform.BuiltIn({
  name: 'topology-from-psf',
  display: { name: 'PSF Topology', description: 'Create topology from PSF.' },
  from: [SO.Format.Psf],
  to: SO.Molecule.Topology,
})({
  apply({ a }) {
    return Task.create('Create Topology', async (ctx) => {
      const topology = await topologyFromPsf(a.data).runInContext(ctx);
      return new SO.Molecule.Topology(topology, { label: topology.label || a.label, description: 'Topology' });
    });
  },
});

export { PsfProvider };
const PsfProvider = DataFormatProvider({
  label: 'PSF',
  description: 'PSF',
  category: TopologyFormatCategory,
  stringExtensions: ['psf'],
  parse: async (plugin, data) => {
    const format = plugin.state.data
      .build()
      .to(data)
      .apply(ParsePsf, {}, { state: { isGhost: true } });
    const topology = format.apply(TopologyFromPsf);

    await format.commit();

    return { format: format.selector, topology: topology.selector };
  },
});

type PsfProvider = typeof PsfProvider;
