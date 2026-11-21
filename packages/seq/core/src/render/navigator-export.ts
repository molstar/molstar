/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { DrawLists } from './draw-lists.js';
import type { Overlay } from './layer-type.js';
import {
  type NavigatorFrame,
  NavigatorRenderer,
  type NavigatorRendererOptions,
  navigatorBox,
  windowText,
} from './navigator.js';
import { drawListsToSvg, escapeXml, type SvgOptions } from './svg.js';
import type { ViewLayout } from './tracks.js';

// The navigators of an export (_spec/README.md "Render", Export): an export's part, given to `exportSvg` by the
// composition that lists navigators (`renderSvg`, `createSeqContext`, a viewer's `navigator` behavior), so a composition
// that lists none bundles none of their drawing.

/** What the navigators of an export draw with: a renderer's options but the navigator, and the frames' common input. */
export interface NavigatorExportInput extends Omit<NavigatorRendererOptions, 'navigator' | 'view'> {
  readonly view: NavigatorRendererOptions['view'];
  readonly dpr: number;
  readonly cvdSafe: boolean;
  /** The store's dynamic nodes now (asked each frame). */
  overlay(): Overlay;
}

/** The included navigators of one export: their frames over the view's layout, their issues, and the SVG around it. */
export interface ExportedNavigators {
  readonly count: number;
  /** Each navigator's frame over the view's layout, at its box around the view, and what they want. */
  frames(layout: ViewLayout, settled: boolean): { readonly frames: readonly NavigatorFrame[]; readonly wants: Want[] };
  /** Their issues, after the last frames. */
  issues(): readonly Issue[];
  /** The view's SVG with the navigators stacked by `place` around it. */
  wrap(
    view: string,
    lists: DrawLists,
    model: StaticModel,
    frames: readonly NavigatorFrame[],
    svg: Pick<SvgOptions, 'fonts' | 'fade' | 'idPrefix'>,
  ): string;
}

/** Makes an export's navigators: the view's shown ones. */
export type NavigatorExport = (input: NavigatorExportInput) => ExportedNavigators;

/** The core's navigators of an export: outlines and wash, the windows' text as titles. */
export const NavigatorSvgExport: NavigatorExport = (input) => {
  const { dpr, cvdSafe, overlay, ...options } = input;
  const navs = input.view.navigators
    .filter((n) => !n.hidden)
    .map((navigator) => new NavigatorRenderer({ ...options, navigator }));
  return {
    count: navs.length,
    frames(layout, settled) {
      const frames = navs.map((n) =>
        n.frame({
          size: navigatorBox(n.navigator, layout.width, layout.height, options.pictures),
          dpr,
          settled,
          overlay: overlay(),
          view: layout,
          ...(cvdSafe ? { cvdSafe: true } : {}),
        }),
      );
      return { frames, wants: frames.flatMap((f) => f.wants) };
    },
    issues: () => navs.flatMap((n) => n.allIssues),
    wrap: (view, lists, model, frames, svg) => withNavigators(view, lists, model, navs, frames, svg),
  };
};

/**
 * The view's SVG with its navigators stacked by `place` around it (top and bottom across its width, left and right
 * along its height, a miniature at its own size): each a group with its windows' `<title>` and `data-navigator` its
 * id; the view's metadata on the outer root.
 */
function withNavigators(
  view: string,
  lists: DrawLists,
  model: StaticModel,
  navs: readonly NavigatorRenderer[],
  frames: readonly NavigatorFrame[],
  svg: Pick<SvgOptions, 'fonts' | 'fade' | 'idPrefix'>,
): string {
  const vw = lists.width;
  const vh = lists.height;
  // each navigator's box, a zoom strip's height below it
  const boxes = navs.map((n) => {
    const b = navigatorBox(n.navigator, vw, vh, n.options.pictures);
    return { width: b.width, height: b.height + n.zoomHeight };
  });
  const sum = (p: string, k: 'width' | 'height') =>
    boxes.reduce((t, b, i) => t + (navs[i].navigator.place === p ? b[k] : 0), 0);
  const top = sum('top', 'height');
  const left = sum('left', 'width');
  const width = left + vw + sum('right', 'width');
  const height = top + vh + sum('bottom', 'height');
  const at = { top: 0, bottom: top + vh, left: 0, right: left + vw };
  const nest = (inner: string, x: number, y: number) => inner.replace(/^<svg /, `<svg x="${x}" y="${y}" `);
  const meta = /^<svg [^>]*?viewBox="[^"]*"([^>]*)>/.exec(view)?.[1] ?? '';
  const parts = [nest(view, left, top)];
  navs.forEach((n, i) => {
    const b = boxes[i];
    const place = n.navigator.place;
    const x = place === 'left' ? at.left : place === 'right' ? at.right : left;
    const y = place === 'top' ? at.top : place === 'bottom' ? at.bottom : top;
    if (place === 'top') at.top += b.height;
    else if (place === 'bottom') at.bottom += b.height;
    else if (place === 'left') at.left += b.width;
    else at.right += b.width;
    const titles = frames[i].windows
      .filter((w) => w.shown)
      .map((w) => `Window: ${windowText(model, n.navigator, w)}`)
      .join('; ');
    const inner = drawListsToSvg(frames[i].lists, { ...svg, idPrefix: `${svg.idPrefix ?? ''}n${i}-` });
    parts.push(
      `<g data-navigator="${escapeXml(n.navigator.id)}"><title>${escapeXml(titles || n.navigator.id)}</title>${nest(inner, x, y)}</g>`,
    );
  });
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}"${meta}>${parts.join('')}</svg>`;
}
