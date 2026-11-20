/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A narrow frame (core's _spec "Render", A narrow frame; _spec/README.md "Composition", Navigators): the vertical
// overview beside the canvas, given the row's width as the default `ui.navigators` gives it (`setRoom`), keeps its 96 px
// in a wide frame, narrows toward 48 px to keep the plot 160 px wide, and hides below that, the rows rail kept; headless
// on a fake clock, the canvas resized to what the overview leaves as the mount's row does.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, overviewNavigator } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Letters = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAA';

function spec(): Node {
  const b = createSeqSpec({ title: 'narrow' });
  const n = 30;
  const label = Array.from({ length: n }, (_, i) => `r${i}`);
  const letters = Array.from({ length: n }, (_, i) => Letters.slice(i % 40, (i % 40) + 40));
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label, letters } });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    s.track('rows', { label: '', height: 100 }).layer('r', { type: 'rows', data: aln });
    v.navigator('nav', overviewNavigator(s));
  });
  return b.build();
}

/** A viewer in a row `frame` px wide: the overview's room the row, the canvas what the overview leaves. */
async function inFrame(frame: number) {
  frames?.restore();
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('view.jumpMs'), 0]] as never,
  });
  viewers.push(viewer);
  await viewer.load(spec());
  viewer.resize(frame - 96);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  const nav = viewer.navigatorModel('nav')!;
  nav.setRoom(frame);
  // the mount's loop: the part sized, the canvas resized to the rest, drawn again (twice: it is stable)
  for (let i = 0; i < 2; i++) {
    nav.draw();
    viewer.resize(frame - nav.width);
    viewer.draw();
  }
  nav.draw();
  return { viewer, nav, columns: viewer.frame!.layout.narrow!.columnsPx };
}

describe('a narrow frame: the vertical overview beside the canvas', () => {
  it('at 1280 px: its own 96 px, the rows rail beside it', async () => {
    const { viewer, nav } = await inFrame(1280);
    expect(nav.width).toBe(96);
    expect(nav.frame).toBeDefined();
    expect(viewer.width).toBe(1280 - 96);
    expect(viewer.frame!.layout.sections[0].rowsRails?.length).toBe(1);
  });

  it('narrowed toward 48 px to keep the plot 160 px wide, then hidden at 480 px; the rows rail kept', async () => {
    const { columns } = await inFrame(1280);
    const narrowed = await inFrame(columns + 160 + 70);
    expect(narrowed.nav.width).toBe(70);
    expect(narrowed.viewer.frame!.layout.narrow!.plotPx).toBeGreaterThanOrEqual(160);
    expect(narrowed.viewer.frame!.layout.narrow!.plotPx).toBeLessThan(161);
    const hidden = await inFrame(columns + 160 + 40);
    expect(hidden.nav.width).toBe(0);
    expect(hidden.nav.draw()).toBeUndefined();
    expect(hidden.viewer.frame!.layout.sections[0].rowsRails?.length).toBe(1);
    // at 480 px: the plot beside it is what decides
    const at480 = await inFrame(480);
    expect(at480.nav.width).toBe(480 - columns - 48 < 160 ? 0 : Math.min(96, 480 - columns - 160));
    expect(at480.viewer.frame!.layout.sections[0].rowsRails?.length).toBe(1);
  });

  it('no room given, or a host’s own size, keeps its width', async () => {
    const { viewer, nav, columns } = await inFrame(columns0());
    expect(nav.width).toBe(0);
    nav.setRoom(0);
    nav.draw();
    expect(nav.width).toBe(96);
    nav.resize(120, 0);
    nav.setRoom(columns + 100);
    nav.draw();
    expect(nav.width).toBe(120);
    expect(viewer.width).toBeGreaterThan(0);
  });
});

/** A frame too narrow for the overview whatever the columns (the plot under 160 px beside its floor). */
const columns0 = () => 300;
