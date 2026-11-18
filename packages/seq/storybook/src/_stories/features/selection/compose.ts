/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqCoreConfig, type SeqPorts } from '@molstar/seqx-core/context';
import type { RenderSvgResult } from '@molstar/seqx-core/render';
import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';

// The selection story's composition, DOM-free: the story and its headless twin call these, so the twin acts on the
// viewer as the story's pointer and keys do (through `viewer.input`, the bindings and the commands, never around them)
// and draws the same SVG beside it.

/** How the selection is drawn: its shape and the tint's place (the core's `emphasis` render items). */
export interface SelectionLook {
  readonly shape?: 'column' | 'tracks';
  readonly place?: 'over' | 'under';
}

/**
 * A viewer of an example view, composed from the default spec with the view's registry and the selection's look as the
 * core's render items (read once: another look is another viewer).
 */
export function composeViewer(
  view: ExampleView,
  options: { target?: HTMLElement; delayMs?: number; ports?: SeqPorts; look?: SelectionLook } = {},
): SeqViewer {
  const base = DefaultSeqViewerSpec();
  const extra = view.registry?.(options.delayMs !== undefined ? { delayMs: options.delayMs } : {}) ?? [];
  const look = options.look ?? {};
  const emphasis = {
    ...(look.shape ? { selectionShape: look.shape } : {}),
    ...(look.place ? { selectionPlace: look.place } : {}),
  };
  return createSeqViewer(
    {
      ...base,
      ...(options.target ? {} : { ui: undefined }),
      registry: [...(base.registry ?? []), ...extra],
      ...(options.ports ? { ports: options.ports } : {}),
      ...(Object.keys(emphasis).length ? { config: [[SeqCoreConfig.render, { emphasis }]] } : {}),
    },
    options.target,
  );
}

/** The middle of a position's cell on a track (sections and tracks by index), in the canvas's CSS px. */
export function cellPoint(viewer: SeqViewer, section: number, track: number, position: number): [number, number] {
  const s = viewer.frame?.layout.sections[section];
  const t = s?.tracks[track];
  if (!s?.axis || !t) throw new Error(`No track ${track} in section ${section} laid out`);
  return [(s.axis.x(position) + s.axis.x(position + 1)) / 2, t.box.y + Math.min(6, t.box.height / 2)];
}

/** A click as the canvas gets it: the pointer over the point, pressed and released there. */
export function clickAt(
  viewer: SeqViewer,
  [x, y]: [number, number],
  keys: { shiftKey?: boolean; ctrlKey?: boolean; metaKey?: boolean } = {},
): void {
  viewer.input.move(x, y, keys);
  viewer.input.down(x, y, keys);
  viewer.input.up();
  viewer.input.click(x, y, keys);
}

/** A drag from one point to another in `steps` moves, as the canvas gets it. */
export function dragFrom(viewer: SeqViewer, from: [number, number], to: [number, number], steps = 4): void {
  viewer.input.move(from[0], from[1]);
  viewer.input.down(from[0], from[1]);
  for (let i = 1; i <= steps; i++)
    viewer.input.move(from[0] + ((to[0] - from[0]) * i) / steps, from[1] + ((to[1] - from[1]) * i) / steps);
  viewer.input.up();
}

/** What the status line says: the selection as an address, the cursor, and what the view shows. */
export function describeSelection(viewer: SeqViewer): string {
  const c = viewer.selection.cursor.value;
  const cursor = c ? `${c.system} ${c.position + 1}` : '—';
  return `selection: ${viewer.addresses.selectionAddress() ?? '—'} · cursor: ${cursor}`;
}

/** The SVG beside the canvas: the viewer's export of the same spec, the selection drawn and projected as the canvas does. */
export function selectionSvg(viewer: SeqViewer, idPrefix: string): Promise<RenderSvgResult> {
  return viewer.renderSvg({ idPrefix });
}
