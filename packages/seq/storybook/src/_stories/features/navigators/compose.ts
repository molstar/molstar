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
  readonly blocks: readonly { readonly x: number; readonly width: number; readonly label: string }[];
  readonly windows: readonly { readonly x: number; readonly width: number; readonly idle: boolean }[];
}

/**
 * The host's own picture of the navigator: p53's domains placed by the host (the navigator's fit: the whole protein
 * across its width) and each window as the model's last frame placed it (its drawn ends: the floor applied).
 */
export function hostPicture(model: NavigatorModel): HostPicture {
  const px = model.width / P53.length;
  const blocks = Domains.id.map((label, i) => ({
    x: Domains.start[i] * px,
    width: (Domains.end[i] - Domains.start[i]) * px,
    label,
  }));
  const windows = (model.frame?.windows ?? [])
    .filter((w) => w.shown)
    .map((w, i) => ({ x: w.dx0, width: w.dx1 - w.dx0, idle: model.values.value[i]?.disabled ?? false }));
  return { blocks, windows };
}
