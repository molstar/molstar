/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type NavigatorExport, NavigatorSvgExport } from '@molstar/seqx-core/render/navigator-export';
import type { Node } from '@molstar/core/repr/tree';
import { isSquareNavigator } from '@molstar/seqx-core/render/navigator-square-load';
import { withViewKind } from '@molstar/seqx-core/resolve';
import type { ResolvedNavigator, ResolvedView } from '@molstar/seqx-core/view';
import { loadSquare, type NavigatorDef, NavigatorModel, squareLoaded } from './navigator.js';
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
  /**
   * What a draw of these navigators waits for (`SeqViewer.settle`): the square's chunk when one is a square and it has
   * not landed (`loadSquare`), else nothing.
   */
  loaded?(view: ResolvedView, navigators: readonly ResolvedNavigator[]): Promise<void> | undefined;
}

/** Whether a view's `navigators` may hold a square (a navigator on both axes bound to two sections). */
function maySquare(navigators: Node): boolean {
  return !!navigators.children?.some((n) => {
    const of = n.params?.of;
    return n.kind === 'navigator' && n.params?.axes === 'both' && Array.isArray(of) && of.length === 2;
  });
}

// a model of a spec whose view may show a square waits for the square's chunk (core's "Resolve", Lazy kinds), so a
// declared square is drawn at its first frame; a host's square (`viewer.navigator(def)`) waits in `settle`
withViewKind('navigators', loadSquare, maySquare);

/** The viewer's navigators: what the default spec lists. */
export const NavigatorsPart: SeqNavigatorsPart = {
  model: (viewer, id, def) => new NavigatorModel(viewer, id, def),
  svgExport: NavigatorSvgExport,
  loaded: (view, navigators) =>
    !squareLoaded() && navigators.some((n) => !n.hidden && isSquareNavigator(n, view)) ? loadSquare() : undefined,
};
