/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { InteractionsParams } from '@molstar/model/props/computed/interactions';
import { Structure } from '@molstar/model/model/structure';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { StateTransformer } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { computeContacts } from '@molstar/interactions-extension/compute';
import { getCustomInteractionData } from '@molstar/interactions-extension/custom';
import type { InteractionElementSchema, StructureInteractions } from '@molstar/interactions-extension/model';
import { buildInteractionsShape, InteractionVisualParams } from '@molstar/interactions-extension/visuals';

const Factory = StateTransformer.builderFactory('interactions-extension');

export class InteractionData extends SO.Create<{ interactions: StructureInteractions }>({
  name: 'Interactions',
  typeClass: 'Data',
}) {}

export const ComputeContacts = Factory({
  name: 'compute-contacts',
  display: 'Compute Contacts',
  from: SO.Molecule.Structure.Selections,
  to: InteractionData,
  params: {
    interactions: PD.Group(InteractionsParams),
  },
})({
  apply({ params, a }) {
    return Task.create('Compute Contacts', async (ctx) => {
      const interactions = await computeContacts(ctx, a.data, { interactions: params.interactions });
      return new InteractionData({ interactions }, { label: 'Interactions' });
    });
  },
});

export const CustomInteractions = Factory({
  name: 'custom-interactions',
  display: 'Custom Interactions',
  from: SO.Root,
  to: InteractionData,
  params: {
    interactions: PD.Value<InteractionElementSchema[]>([], { isHidden: true }),
  },
})({
  apply({ params, dependencies }) {
    return Task.create('Custom Interactions', async (ctx) => {
      const structures: { [ref: string]: Structure } = {};
      for (const [k, v] of Object.entries(dependencies ?? {})) {
        structures[k] = v.data as Structure;
      }
      const interactions = getCustomInteractionData(params.interactions, structures);
      return new InteractionData({ interactions }, { label: 'Custom Interactions' });
    });
  },
});

export const InteractionsShape = Factory({
  name: 'interactions-shape',
  display: { name: 'Interactions Shape' },
  from: InteractionData,
  to: SO.Shape.Provider,
  params: InteractionVisualParams,
})({
  canAutoUpdate: () => true,
  apply({ a, params }) {
    return new SO.Shape.Provider(
      {
        label: 'Interactions Shape Provider',
        data: { interactions: a.data.interactions, params },
        params: PD.withDefaults(Mesh.Params, {}),
        getShape: (_, data: { interactions: StructureInteractions; params: InteractionVisualParams }, __, prev: any) =>
          buildInteractionsShape(data.interactions, data.params, prev?.geometry),
        geometryUtils: Mesh.Utils,
      },
      { label: 'Interactions Shape Provider' },
    );
  },
});
