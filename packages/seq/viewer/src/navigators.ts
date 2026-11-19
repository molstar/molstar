/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type NavigatorExport, NavigatorSvgExport } from '@molstar/seqx-core/render/navigator-export';
import { type NavigatorDef, NavigatorModel } from './navigator.js';
import type { SeqViewer } from './viewer.js';

// The viewer's navigators part (_spec/README.md "Composition", Navigators): their models and their part of an SVG
// export, listed by the spec (`navigators`), so a viewer that lists none bundles none of them. Without it a declared
// navigator is not drawn, the ruler rail it would replace stays, and one issue says so
// (`seq.viewer.navigators.none`); the gestures are the `navigator` behavior's (without it, drawn with no targets).

export interface SeqNavigatorsPart {
  /** A navigator's model (a declared one by its id, or a host's definition). */
  model(viewer: SeqViewer, id: string, def?: NavigatorDef): NavigatorModel;
  /** How an export draws the navigators it includes. */
  readonly svgExport: NavigatorExport;
}

/** The viewer's navigators: what the default spec lists. */
export const NavigatorsPart: SeqNavigatorsPart = {
  model: (viewer, id, def) => new NavigatorModel(viewer, id, def),
  svgExport: NavigatorSvgExport,
};
