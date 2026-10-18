/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `view.exportSvg` action (seqspec:viewer#actions; the default UI's Export part): `renderSvg` over the live
// context, its string handed to the host's `offer` (a download in the UI). The host registers it with the platform's
// pieces, so this entry stays platform-free: `{ actions: [exportSvgAction({ ports: { deflate }, offer })] }` in its
// registry. The live context is left as it was, its sections' size and view and its measurer, families and current view
// (`renderSvg` borrows it and puts them back before it writes the file), and exports queue per context, since they
// share its viewports. Run as data it is
// `{ command: "seqv:runAction", params: { name: "view.exportSvg", target: { view, width } } }`.

import type { ActionEntry } from '../actions.js';
import type { SeqViewerContext } from '../context.js';
import type { SvgPorts } from './backend.js';
import { type RenderSvgOptions, renderSvg } from './render.js';

export interface ExportSvgTarget {
  /** The view to export; default the current one. */
  readonly view?: string;
  /** CSS px, every section's width. */
  readonly width: number;
  /** The file's pixels per CSS px; default 1. */
  readonly scale?: number;
}
export interface ExportSvgOptions {
  /** The browser's pieces: `deflate` over `CompressionStream("deflate")` for the images; absent: stored PNGs. */
  readonly ports?: SvgPorts;
  /** What the figure shows beyond the tracks (`overlay`, `nav`, `background`, `title`, `gap`), and the em boxes for `dy`. */
  readonly render?: Pick<RenderSvgOptions, 'overlay' | 'nav' | 'background' | 'title' | 'gap' | 'metrics'>;
  /** The finished file, named after its view. */
  offer(svg: string, name: string): void | Promise<void>;
}

const queues = new WeakMap<object, Promise<void>>();
const positive = (n: unknown) => typeof n === 'number' && Number.isFinite(n) && n > 0;

export function exportSvgAction(o: ExportSvgOptions): ActionEntry {
  const entry: ActionEntry<ExportSvgTarget> = {
    name: 'view.exportSvg',
    title: 'Export as SVG',
    when: (ctx, t) =>
      t !== null &&
      typeof t === 'object' &&
      positive(t.width) &&
      (t.scale === undefined || positive(t.scale)) &&
      (t.view === undefined || ctx.views.list.value.some((v) => v.ref === t.view)),
    run(ctx: SeqViewerContext, t) {
      const name = `${t.view ?? ctx.views.current.value ?? 'view'}.svg`;
      const job = (queues.get(ctx) ?? Promise.resolve())
        .then(() => renderSvg(ctx, { ...o.render, width: t.width, scale: t.scale, view: t.view, ports: o.ports }))
        .then((svg) => o.offer(svg, name))
        .catch((e) => ctx.report('seqv.command.error', `view.exportSvg failed: ${e instanceof Error ? e.message : e}`));
      queues.set(ctx, job);
    },
  };
  return entry as ActionEntry;
}
