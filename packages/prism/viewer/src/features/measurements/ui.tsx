/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { MeasurementsAdd } from '../../views/structure/measurements-add.js';

// The Measure section of the Selection tab (parity.md §2; the mock's `registerMockMeasurements`): the picks, what can be
// measured from them, the measurements made and their options, over the `measure` area as the person. The measurement
// nodes of the scene tree are the scene source's (not this feature's).

export function measurementsUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: 'measurements.section.measure',
        slot: 'selection',
        title: 'Measure',
        icon: 'ruler',
        keywords: ['distance', 'angle', 'dihedral', 'label', 'orientation', 'plane'],
        order: 2,
        render: domView(() => {
          const measure = ctx.api.measure;
          return measure ? <MeasurementsAdd model={measure.model(ctx.api)} /> : <></>;
        }),
      },
    ],
  };
}
