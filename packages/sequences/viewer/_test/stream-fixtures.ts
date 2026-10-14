/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1's streamed genome (core/_test/format/stream-provider.ts): a complex from a letters stream (GRCh38's 24 chromosomes,
// 3.1 G positions, or a few given ones), a values stream over it, a ruler and a values track, on a context whose
// loaders hold the two scripted providers. Shared by the stream port's tests, its heap guards and the hg38 bench.

import type { DataPatch, Loaders, SeqViewSpecDoc, StreamNode, StreamRequest } from '@molstar/seq-core/format';
import { GRCH38, hg38Streams, sampledValues, streamedDoc } from '../../core/_test/format/stream-provider.js';
import { type LayerRenderer, SeqViewerConfig } from '../src/index.js';
import type { Streams } from '../src/streams.js';
import { type Opened, opened, preparedRenderers } from './layout-fixtures.js';

export { streamedDoc };
export const GENOME = 'browser.main';

export interface StreamedGenome extends Opened {
  readonly streams: ReturnType<typeof hg38Streams>;
  readonly manager: Streams;
  /** The values layer's ref. */
  readonly layer: string;
}

/** The streamed genome opened at `size` (default 1,000 × 200), the values provider gated or not. */
export async function openGenome(
  o: {
    gated?: boolean;
    chromosomes?: readonly (readonly [string, number])[];
    size?: readonly [number, number];
    renderers?: LayerRenderer[];
    answer?: NonNullable<Parameters<typeof hg38Streams>[0]>['answer'];
    nodes?: NonNullable<Parameters<typeof hg38Streams>[0]>['nodes'];
    doc?: SeqViewSpecDoc;
    /** The genome answers letters (ST3a). */
    letters?: boolean;
  } = {},
): Promise<StreamedGenome> {
  const streams = hg38Streams({
    gated: o.gated,
    chromosomes: o.chromosomes,
    answer: o.answer,
    nodes: o.nodes,
    letters: o.letters,
  });
  const loaders: Loaders = { load: (r) => Promise.reject(new Error(`no file ${r.uri}`)), stream: streams.stream };
  const t = await opened(o.doc ?? streamedDoc(), {
    renderers: o.renderers ?? preparedRenderers(),
    loaders,
    size: { [GENOME]: o.size ?? [1000, 200] },
  });
  const m = t.ctx.graph.assembled.value?.view;
  let layer = '';
  for (let k = 0; m && k < m.resolved.layers.count; k++)
    if (m.set.strings.get(m.resolved.layers.type[k]) === 'values') layer = m.set.strings.get(m.resolved.layers.ref[k]);
  return { ...t, streams, manager: t.ctx.streams as Streams, layer };
}

export const { MaxInFlight, CacheBytes, Margin } = SeqViewerConfig.Stream;

/**
 * A composer frame over the streamed genome (every chromosome; DPR 2, 1,000 × 200) with its regions held: at the whole
 * view, 10 kb or 100 positions (letters drawn) inside chr2, every request answered first; `frame` pans ±9 px among held
 * regions, one pass and one draw into the recording backend. `requests`: what the providers have been asked so far.
 * `letters` (ST3a): with the sequence track, its letters streamed.
 */
export async function genomeFrame(at: 'whole' | '10 kb' | '100 b', draw = true, dynamic = false, letters = false) {
  const { composed, recording } = await import('./composer-fixtures.js');
  const { Dirty } = await import('../src/draw/composer.js');
  // ST1b-a: with the dynamic tracks, the manifest's genes and every data region's peak held.
  const streams = dynamic ? hg38Streams({ nodes: [GENES], answer: withPeaks() }) : hg38Streams({ letters });
  const loaders: Loaders = { load: (r) => Promise.reject(new Error(`no file ${r.uri}`)), stream: streams.stream };
  const doc = dynamic ? dynamicDoc() : streamedDoc({ sequence: letters });
  const t = await composed(doc, { width: 1000, height: 200, dpr: 2, loaders });
  const manager = t.ctx.streams as Streams;
  const vp = t.vp;
  if (!vp) throw new Error('no viewport');
  if (at === 'whole') vp.fit();
  else if (at === '10 kb') vp.zoomTo(1, 100_000_000, 100_010_000);
  else vp.zoomTo(1, 100_040_100, 100_040_200);
  const be = recording();
  const flush = () => new Promise((ok) => setTimeout(ok, 0));
  for (let k = 0; k < 8; k++) {
    for (const dx of [9, -9]) {
      vp.panBy(dx, 0, 'pointer');
      t.layout(GENOME);
      manager.pass();
      await flush();
    }
    // The nodes' frame and its graph update (a new view), then the frame drawn.
    t.scheduler.advance(0);
    await t.ctx.graph.update();
    t.frame(be, Dirty.all);
  }
  let i = 0;
  return {
    manager,
    backend: be,
    requests: () => streams.phylop.requests.length + streams.genome.requests.length,
    frame: () => {
      // A pan (small ints: a double passed to a call is boxed), its pass, and the frame drawn.
      vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
      manager.pass();
      if (draw) t.frame(be, Dirty.all);
    },
  };
}

// --- ST1b-a: dynamic nodes ----------------------------------------------------------------------------------------

/** A region of data (level -1) is this many positions (streams.ts). */
export const DATA_SPAN = 4096;
/** The peak node a data region of `system` carries: a ranges annotation on its chromosome, keyed by the region. */
export function peakNode(system: string, start: number, over: Record<string, unknown> = {}): StreamNode {
  const tx = Math.floor(start / DATA_SPAN);
  const at = tx * DATA_SPAN + 100;
  return {
    kind: 'annotation',
    key: `peak-${system}-${tx}`,
    params: {
      geometry: 'ranges',
      target: `hg38/${system}`,
      label: `peak ${tx}`,
      value: { start: [at], stop: [at + 300] },
      ...over,
    },
  };
}
/** The manifest's node: genes on chr1, known up front. */
export const GENES: StreamNode = {
  kind: 'annotation',
  key: 'genes',
  params: {
    geometry: 'ranges',
    target: 'hg38/chr1',
    value: { start: [1_000_000, 5_000_000], stop: [1_050_000, 5_100_000] },
  },
};
/** The peak a 10 kb window at 10 Mb on chr1 brings (its region's), which `dynamicDoc` names. */
export const PEAK_10MB = 'phylop_src/peak-chr1-2441';

/** The streamed genome with two tracks of dynamic nodes: the manifest's genes, and the peak at 10 Mb on chr1. */
export function dynamicDoc(): SeqViewSpecDoc {
  const doc = streamedDoc();
  const view = doc.root.children?.find((c) => c.kind === 'view');
  const section = view?.children?.find((c) => c.kind === 'section');
  const track = (ref: string, data: string) => ({
    kind: 'track',
    ref,
    params: { label: ref },
    children: [{ kind: 'layer', ref: `${ref}.l`, params: { type: 'ranges', data } }],
  });
  section?.children?.push(track('genes', 'phylop_src/genes'), track('peak', PEAK_10MB));
  return doc;
}

/** Each data region's answer: sampled values, with its peak node (and `extra`'s). */
export const withPeaks =
  (extra?: (r: StreamRequest) => StreamNode[]) =>
  (r: StreamRequest): DataPatch | null => {
    if (r.level >= 0) return null;
    const p = sampledValues(
      r,
      GRCH38.map((c) => c[0]),
    );
    return { ...p, nodes: [peakNode(r.system as string, r.start), ...(extra?.(r) ?? [])] };
  };
