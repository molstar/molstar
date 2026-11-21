/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInLayerTypes, BuiltInLayouts, BuiltInNavigatorPictures } from '@molstar/seqx-core/render/built-ins';
import { MoleculeLettersType } from '@molstar/seqx-core/render/letters-molecule';
import { NavigatorSvgExport } from '@molstar/seqx-core/render/navigator-export';
import { composeSeqContext, type SeqContext, type SeqContextSpec, type SeqRegistryEntry } from './context.js';

// The default headless composition (_spec/README.md "Context"): the core's layer types and layouts listed before the
// spec's own registry, and its navigators in the SVG export. A smaller composition calls `composeSeqContext` with its
// own list and bundles only that.

/** The core's layer types, section layouts and navigator pictures as registry entries: what the default compositions list. */
export const BuiltInRenderRegistry: readonly SeqRegistryEntry[] = [
  // the letters with molecules (their labels, modified residues' marks) in the core letters' place
  ...Object.values({ ...BuiltInLayerTypes, letters: MoleculeLettersType }).map((type) => ({
    kind: 'layer-type' as const,
    type,
  })),
  ...Object.values(BuiltInLayouts).map((layout) => ({ kind: 'layout' as const, layout })),
  // the navigators' pictures (_spec/README.md "Render", Navigators): the axis, the tracks small, the bands, the rows
  ...Object.values(BuiltInNavigatorPictures).map((picture) => ({ kind: 'navigator-picture' as const, picture })),
];

/**
 * The headless context of the default composition: `composeSeqContext` with the built-in layer types and layouts
 * registered before the spec's own entries (the same object again counts, another under a built-in's name is an error)
 * and the built-in navigators drawn in its exports.
 */
export function createSeqContext(spec: SeqContextSpec = {}): SeqContext {
  const context = composeSeqContext({ ...spec, registry: [...BuiltInRenderRegistry, ...(spec.registry ?? [])] });
  // in place: the context's tables are getters (`extend` swaps them), which a spread would freeze
  const render = context.renderSvg;
  return Object.assign(context, {
    renderSvg: (s, options) => render(s, { navigatorExport: NavigatorSvgExport, ...options }),
  } satisfies Pick<SeqContext, 'renderSvg'>);
}
