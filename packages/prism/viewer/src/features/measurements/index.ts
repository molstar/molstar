/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { operationCommand, runAsCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import { MeasureApi } from './operations.js';

declare module '../../api/area.js' {
  // The feature's area (typed as present, as in the default spec: api/area.ts `BuiltInApis`).
  interface BuiltInApis {
    measure: MeasureApi;
  }
}

// The `measurements` feature (viewer-spec.md §3.1): the `measure` area over Mol*'s measurement manager, and
// the Measure and Label section of the Selection tab. Its views are a chunk of their own.

const SelectionTab = 'core.tab.selection';

export const MeasurementsFeature: ViewerExtension = {
  id: 'measurements',
  title: 'Measurements',
  builtIn: true,
  api: { measure: (ctx) => new MeasureApi(ctx) },
  commands: (api) => {
    const measure = api.measure;
    if (!measure) return [];
    return [
      operationCommand(api, measure.distance, {
        category: 'Measure',
        view: { title: 'Selection › Measure and Label', tab: 'selection', section: 'measurements.section.measure' },
        icon: 'measure-distance',
        title: 'Measure a distance…',
        alias: ['distance', 'measure distance'],
        action: 'Measure',
      }),
      operationCommand(api, measure.angle, {
        category: 'Measure',
        view: { title: 'Selection › Measure and Label', tab: 'selection', section: 'measurements.section.measure' },
        icon: 'measure-angle',
        title: 'Measure an angle…',
        alias: ['angle', 'measure angle'],
        action: 'Measure',
      }),
      operationCommand(api, measure.dihedral, {
        category: 'Measure',
        view: { title: 'Selection › Measure and Label', tab: 'selection', section: 'measurements.section.measure' },
        icon: 'measure-angle',
        title: 'Measure a dihedral…',
        alias: ['dihedral', 'torsion'],
        action: 'Measure',
      }),
      operationCommand(api, measure.label, {
        category: 'Measure',
        view: { title: 'Selection › Measure and Label', tab: 'selection', section: 'measurements.section.measure' },
        icon: 'tag',
        title: 'Label an atom or residue…',
        alias: ['label', 'annotate'],
        action: 'Add',
      }),
      operationCommand(api, measure.orientation, {
        category: 'Measure',
        view: { title: 'Selection › Measure and Label', tab: 'selection', section: 'measurements.section.measure' },
        icon: 'move-3d',
        title: 'Show the principal axes…',
        alias: ['orientation', 'principal axes', 'box'],
        action: 'Show',
      }),
      operationCommand(api, measure.plane, {
        category: 'Measure',
        view: { title: 'Selection › Measure and Label', tab: 'selection', section: 'measurements.section.measure' },
        icon: 'square',
        title: 'Show the best-fit plane…',
        alias: ['plane', 'best fit plane'],
        action: 'Show',
      }),
      operationCommand(api, measure.remove, {
        category: 'Measure',
        icon: 'trash-2',
        title: 'Remove a measurement…',
        alias: ['remove measurement', 'delete measurement'],
        action: 'Remove',
      }),
      operationCommand(api, measure.visibility, {
        category: 'Measure',
        icon: 'eye',
        title: 'Show or hide a measurement…',
        alias: ['hide measurement', 'show measurement'],
      }),
      operationCommand(api, measure.options, {
        category: 'Measure',
        icon: 'sliders-horizontal',
        title: 'Measurement options…',
        alias: ['unit', 'text color'],
        action: 'Set',
      }),
      operationCommand(api, measure.get, {
        category: 'Measure',
        icon: 'list',
        title: 'List the measurements',
        alias: ['measurements'],
      }),
      {
        id: 'measure.open',
        title: 'Measure…',
        category: 'Measure',
        icon: 'ruler',
        alias: ['measure', 'measurements'],
        // The section is in the Selection panel: without the panel there is nowhere to open.
        when: () => (api.view.hasTab(SelectionTab) ? true : 'There is no Selection panel.'),
        async run() {
          const r = await runAsCommand(api, 'view.open', {
            tab: SelectionTab,
            section: 'measurements.section.measure',
          });
          return { ...r, undoable: false };
        },
      },
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).measurementsUI(ctx),
};
