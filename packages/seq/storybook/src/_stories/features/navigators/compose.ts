/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { NavigatorPicture } from '@molstar/seqx-core/extension';
import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import type { ExampleView } from '@molstar/seqx-examples';
import { addDomains, conservation, Domains, P53 } from '@molstar/seqx-examples/p53/data';
import type { NavigatorModel } from '@molstar/seqx-viewer/navigator';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The navigators' stories' composition, DOM-free: the Strip story and its headless twin call these; the coverage strip
// (a host's picture, core's "Render", Navigators, Pictures) and the host-drawn navigator (the viewer's NavigatorModel as
// a host's contract) with their headless twins.

/** Where the Strip story opens: the DNA-binding domain's first half (a window smaller than the strip: it can be dragged). */
export const StripOpening = '/data/p53/(range :from 94 :to 200)';

/** Moves the viewer to where the Strip story opens. */
export function openStrip(viewer: SeqViewer): void {
  viewer.runCommand('address.focus', { address: StripOpening });
}

/** p53's conservation, as the spec holds it and the readout reads it. */
const Conservation = conservation();

/**
 * A host's picture, `host:coverage`: the values its `data` param names as one track filling the navigator (its zoom
 * strip the axis), its readout the position and the value there.
 */
export const CoveragePicture: NavigatorPicture = {
  name: 'host:coverage',
  axes: ['positions'],
  size: 32,
  fills: true,
  picture: (c) => {
    const data = c.navigator.params.data;
    if (c.zoom || typeof data !== 'string') return undefined;
    const track: Node = {
      kind: 'track',
      id: 'coverage',
      params: { label: '', height: c.height },
      children: [{ kind: 'layer', id: 'coverage', params: { type: 'values', data } }],
    };
    return { tracks: [track] };
  },
  // the host holds its data: the value under the pointer read from it
  readout: (p) => {
    const v = Conservation[p];
    return `residue ${p + 1}${typeof v === 'number' && !Number.isNaN(v) ? ` · conservation ${v.toFixed(2)}` : ''}`;
  },
};

/** p53's tracks under a coverage strip of its conservation. */
function coverageSpec(): Node {
  const b = createSeqSpec({ title: 'p53 with a coverage strip' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  const values = b.data.values('conservation', { target: p53, values: Conservation });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    s.track('conservation', { label: 'Conservation' }).layer('values', { type: 'values', data: values });
    // the host's picture, its param plain (a reader without it draws the axis, the param kept)
    v.navigator('coverage', { ...stripNavigator(s), content: 'host:coverage', data: values } as never);
  });
  return b.build();
}

/** The coverage strip as an example view: its picture registered beside the default registry. */
export const CoverageView: ExampleView = {
  dataset: 'p53',
  id: 'coverage-strip',
  title: 'Coverage strip',
  summary:
    'A host’s navigator picture, `host:coverage`: p53’s conservation filling the strip above its tracks, the window the view’s.',
  build: coverageSpec,
  registry: () => [{ kind: 'navigator-picture', picture: CoveragePicture }],
};

/** The host-drawn navigator's spec: p53's letters and domains, no declared navigator (the host adds its own). */
export function hostDrawnSpec(): Node {
  const b = createSeqSpec({ title: 'p53, a host-drawn navigator' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = addDomains(b, p53);
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}

export const HostSection = '/views/main/protein';

/** What a host draws of its navigator: its domains as blocks and the windows, x in the navigator's px. */
export interface HostPicture {
  /** A block's lane: the first free one at its x, so blocks that overlap (a domain inside another) never share one. */
  readonly blocks: readonly {
    readonly x: number;
    readonly width: number;
    readonly label: string;
    /** What it is called in full, and the residues it covers (1-based, inclusive). */
    readonly name: string;
    readonly first: number;
    readonly last: number;
    readonly lane: number;
  }[];
  /** The lanes the blocks use (the host divides its height by them). */
  readonly lanes: number;
  readonly windows: readonly {
    readonly section: string;
    readonly x: number;
    readonly width: number;
    readonly idle: boolean;
  }[];
}

/** The block under a point of the host's picture (x, y in the navigator's px, `height` its height), or undefined. */
export function blockAt(
  pic: HostPicture,
  x: number,
  y: number,
  height: number,
): HostPicture['blocks'][number] | undefined {
  const lane = (height - 16) / pic.lanes;
  return pic.blocks.find((b) => x >= b.x && x < b.x + b.width && y >= 8 + b.lane * lane && y < 8 + (b.lane + 1) * lane);
}

/** The readout of a hovered block: its name and residues, then the model's own words about the pointer's place. */
export function blockReadout(
  block: HostPicture['blocks'][number] | undefined,
  modelReadout: string | undefined,
): string {
  if (!block) return modelReadout ?? 'Hover the navigator';
  return `${block.name} (${block.label}), residues ${block.first}–${block.last}${modelReadout ? ` · ${modelReadout}` : ''}`;
}

/**
 * The host's own picture of the navigator after a `draw`: p53's domains placed through the model (`xRange`: the
 * navigator's display order, whatever the direction or the pieces) and each window as the model's last frame placed it
 * (its drawn ends: the floor applied), idle as its slider says.
 */
export function hostPicture(model: NavigatorModel): HostPicture {
  const placed = Domains.id.flatMap((label, i) => {
    const r = model.xRange(Domains.start[i], Domains.end[i]);
    return r
      ? [
          {
            x: r[0],
            width: r[1] - r[0],
            label,
            name: Domains.label[i],
            first: Domains.start[i] + 1,
            last: Domains.end[i],
          },
        ]
      : [];
  });
  // each block the first lane whose last block ended before it: overlapping labels never overprint (F27)
  const ends: number[] = [];
  const blocks = [...placed]
    .sort((a, b) => a.x - b.x)
    .map((b) => {
      let lane = ends.findIndex((end) => end <= b.x);
      if (lane < 0) lane = ends.length;
      ends[lane] = b.x + b.width;
      return { ...b, lane };
    });
  const windows = (model.frame?.windows ?? [])
    .filter((w) => w.shown)
    .map((w) => ({
      section: w.section,
      x: w.dx0,
      width: w.dx1 - w.dx0,
      idle: model.values.value.find((v) => v.section === w.section)?.disabled ?? false,
    }));
  return { blocks, lanes: Math.max(1, ends.length), windows };
}
