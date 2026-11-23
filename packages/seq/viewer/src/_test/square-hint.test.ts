// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { mountNavigator, SquareHint } from '../mount-navigator.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

// The square's hint (_spec/README.md "Composition", Navigators, The square's hint): one muted line under the square
// naming its gestures; `hint: false` (or `navigator.squareHint: false`) leaves it out.

let frames: Frames;
let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});

beforeEach(() => {
  frames = installFrames();
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames.restore();
  document.body.replaceChildren();
});

const Ring = '/views/main/ring';

async function make(squareHint?: boolean) {
  const b = createSeqSpec();
  const g = b.data.assembly(
    'g',
    assembly(
      'Four',
      ['chr1', 'chr2', 'chr3', 'chr4'].map((name) => ({ name, length: 1_000_000 })),
    ),
  );
  b.views.view('main', (v) => {
    genome(v, 'ring', g, { layout: 'rings' } as never);
    v.navigator('square', { of: [Ring, Ring], axes: 'both', system: g, switch: 'include', size: [200, 200] } as never);
  });
  viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    ...(squareHint === undefined ? {} : { config: [[viewerItem('navigator.squareHint'), squareHint]] }),
  });
  await viewer.load(b.build());
  await viewer.settle();
  return viewer;
}

describe('the square’s hint', () => {
  it('one muted line under the square, taken away with it', async () => {
    const v = await make();
    const host = document.createElement('div');
    document.body.append(host);
    const dispose = mountNavigator(host, v, 'square');
    const hint = host.querySelector('.seqx-navigator-hint')!;
    expect(hint.textContent).toBe('Drag to move · drag an edge to resize · double-click for the whole genome');
    expect(hint.textContent).toBe(SquareHint);
    expect(host.lastElementChild).toBe(hint);
    dispose();
    expect(host.querySelector('.seqx-navigator-hint')).toBeNull();
  });

  it('`hint: false`, or the configuration item, leaves it out', async () => {
    const v = await make();
    const host = document.createElement('div');
    document.body.append(host);
    mountNavigator(host, v, 'square', { hint: false });
    expect(host.querySelector('.seqx-navigator-hint')).toBeNull();
    expect(host.querySelector('[data-navigator]')).not.toBeNull();
    v.retire();
    const off = await make(false);
    const other = document.createElement('div');
    document.body.append(other);
    mountNavigator(other, off, 'square');
    expect(other.querySelector('.seqx-navigator-hint')).toBeNull();
  });
});
