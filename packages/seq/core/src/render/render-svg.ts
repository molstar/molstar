/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { BuiltInLayerTypes, BuiltInLayouts } from './built-ins.js';
import { exportSvg, type RenderSvgOptions, type RenderSvgResult } from './export.js';

// The headless default (_spec/README.md "Render", Export): `exportSvg` with everything the core ships, the built-in
// layer types, layouts and navigators (their drawing and pictures a lazy chunk, loaded for an export that includes
// them), unless the options name their own. A composition that lists fewer calls
// `exportSvg` (a viewer, `createSeqContext`'s base) and bundles fewer.

/**
 * Whether a spec holds molecules, complexes or modified residues: what the letters with molecules draw (a chunk of
 * their own, loaded only then: a spec of plain letters never reads it).
 */
function molecular(spec: Node): boolean {
  for (const d of spec.children ?? [])
    if (d.kind === 'data')
      for (const n of d.children ?? [])
        if (n.kind === 'molecule' || n.kind === 'complex' || n.params?.modified !== undefined) return true;
  return false;
}

/** A spec drawn to an SVG string once the data its view needs has settled (or the timeout passed): the built-ins. */
/** Whether an export may include navigators: asked to, or a view's `export.navigators` is `include` and not refused. */
function navigators(spec: Node, options: RenderSvgOptions): boolean {
  if (options.navigators) return options.navigators === 'include';
  for (const c of spec.children ?? [])
    if (c.kind === 'views')
      for (const v of c.children ?? [])
        if ((v.params?.export as { navigators?: string } | undefined)?.navigators === 'include') return true;
  return false;
}

export async function renderSvg(spec: Node, options: RenderSvgOptions): Promise<RenderSvgResult> {
  const types = molecular(spec)
    ? { ...BuiltInLayerTypes, letters: (await import('./letters-molecule.js')).MoleculeLettersType }
    : BuiltInLayerTypes;
  // the navigators' drawing and pictures, loaded only for an export that includes them
  const navs = navigators(spec, options)
    ? await Promise.all([import('./navigator-export.js'), import('./navigator-picture.js')])
    : undefined;
  return exportSvg(spec, {
    types,
    layouts: BuiltInLayouts,
    ...(navs ? { pictures: navs[1].BuiltInNavigatorPictures, navigatorExport: navs[0].NavigatorSvgExport } : {}),
    ...options,
  });
}
