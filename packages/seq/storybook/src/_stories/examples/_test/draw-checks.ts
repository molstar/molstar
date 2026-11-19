/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Frame, type RangesLay, recordDrawLists } from '@molstar/seqx-core/render';
import type { ExampleView } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { uiHistory } from '@molstar/seqx-viewer/ui-history';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { readFileSync } from 'node:fs';
import { vi } from 'vitest';

// The checks the example views' drawing tests share: every drawn mark of a ranges layer is hit (and outlined when
// hovered), and every item of the data in view has its marks.

export const HoverOwner = 'test';

export interface OpenOptions {
  /** The gene file's text (the views that read it by location). */
  readonly gene: string;
  /** Each answer's delay, or `staggered` (the story's own: later regions later). */
  readonly delayMs?: number | 'staggered';
  /** How long the fake clock runs after the load. */
  readonly wait?: number;
  /** Keep the real `Date` and `performance` (a bound on time reads them): only the timers are faked. */
  readonly realClock?: boolean;
  /** The host's history, for edits that are undone (the spec takes part in it). */
  readonly history?: boolean;
}

/** A viewer of an example view on a fake clock, composed as the stories compose it: loaded, the clock run, drawn. */
export async function openExample(view: ExampleView, options: OpenOptions): Promise<SeqViewer> {
  const { gene, delayMs = 300, wait = 6000, realClock = false, history = false } = options;
  vi.useFakeTimers({
    toFake: realClock ? ['setTimeout', 'clearTimeout'] : ['setTimeout', 'clearTimeout', 'Date', 'performance'],
  });
  const text = view.id === 'gene-by-region' || view.id === 'gene-region-ideogram' ? regionSeqText() : gene;
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    ...(history ? { history: uiHistory() } : {}),
    registry: [...(base.registry ?? []), ...(view.registry?.(delayMs === 'staggered' ? {} : { delayMs }) ?? [])],
    ports: {
      assets: () => {
        // every dataset’s files (`fixtures/…`) are read as they are committed, at once (a read of the disk is not on the
        // fake clock: the layout waits for the letters); every other location answers the gene's text
        const read = (uri: string) => readFileSync(`${PackageFolder}${uri}`);
        const fixture = (uri: string) => uri.startsWith('fixtures/');
        return {
          reader: {
            readText: async (uri) => (fixture(uri) ? read(uri).toString('utf8') : text),
            readBytes: async (uri) => (fixture(uri) ? new Uint8Array(read(uri)) : new TextEncoder().encode(text)),
          },
        };
      },
    },
  });
  v.resize(720);
  void v.load(view.build());
  await vi.advanceTimersByTimeAsync(wait);
  v.draw();
  return v;
}

/** The hover pass's rects. */
const hoverRects = (frame: Frame) =>
  recordDrawLists(frame.lists, ['hover'])
    .split('\n')
    .filter((l) => l.startsWith('rect ')).length;

/** Each drawn mark of every ranges layer hit at its centre: what is wrong, and how many marks were checked. */
export function checkMarks(v: SeqViewer) {
  const wrong: string[] = [];
  let marks = 0;
  let merged = 0;
  const layers = v.frame!.layout.sections.flatMap((s) => s.tracks.flatMap((t) => t.layers));
  for (const l of layers) {
    if (l.type.name !== 'ranges') continue;
    const lay = l.lay as RangesLay;
    for (let j = 0; j < lay.m; j++) {
      marks++;
      const x = (lay.mx0[j] + lay.mx1[j]) / 2;
      const y = lay.box.y + (lay.mlane[j] + 0.5) * lay.laneHeight;
      // the items the mark draws: its own, or those of the run's marks
      const drawn = new Set<string | number>();
      if (lay.mitem[j] >= 0) drawn.add(lay.keys[lay.mitem[j]]);
      else {
        merged++;
        for (let k = 0; k < (lay.hn ?? 0); k++)
          if (lay.hlane[k] === lay.mlane[j] && lay.hx1[k] >= lay.mx0[j] && lay.hx0[k] <= lay.mx1[j])
            drawn.add(lay.keys[lay.hitem[k]]);
      }
      const hit = v.hit(x, y);
      const where = `${l.key} mark ${j} (lane ${lay.mlane[j]}, x ${lay.mx0[j]}–${lay.mx1[j]})`;
      if (hit?.layer !== l.key || hit.element === undefined || !drawn.has(hit.element.key)) {
        wrong.push(`${where}: ${hit?.layer} ${JSON.stringify(hit?.element)}`);
        continue;
      }
      // hovered, the item is outlined
      v.hover.set(hit, HoverOwner);
      const f = v.draw() ?? v.frame!;
      if (!hoverRects(f)) wrong.push(`${where}: hovered ${hit.element.key}, nothing outlined`);
      v.hover.set(undefined, HoverOwner);
    }
  }
  return { wrong, marks, merged };
}
