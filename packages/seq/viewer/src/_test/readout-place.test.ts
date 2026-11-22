/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A navigator part's readout is never over the control in use (_spec/README.md "Composition", Navigators): beside the
// window on the side with room, else above or below the part; during a drag its box touches neither the window (its
// grips included) nor the pointer.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { GripReach, placeReadout, PointerReach, pointerBox, type Box } from '../readout-place.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { installFrames, type Frames } from './frames.js';

const Letters =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD';

const P = '/views/v/protein';

/** p53's letters (a ruler with them, or none), a strip over the section or over two independent sections. */
function spec(
  options: { ruler?: boolean; twin?: boolean; group?: boolean; strip?: boolean; follow?: boolean; same?: boolean } = {},
) {
  const b = createSeqSpec({ title: 'navigators' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters });
  const q = b.data.sequence('q', { alphabet: 'protein', letters: Letters });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p });
    if (options.ruler) s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('letters').layer('letters', { type: 'letters' });
    let t: { address: string } | undefined;
    if (options.twin) {
      // independent: another system; in a group: the same system (one range); a follower: follows the protein
      t = v.section('twin', {
        of: options.group || options.follow || options.same ? p : q,
        ...(options.follow ? { follow: s } : {}),
      });
      (t as ReturnType<typeof v.section>).track('letters').layer('letters', { type: 'letters' });
    }
    if (options.strip !== false) v.navigator('strip', stripNavigator(t && !options.group ? [s, t.address] : s));
  });
  return b.build();
}

async function make(root: Node, config: readonly (readonly [unknown, unknown])[] = []) {
  frames = installFrames();
  const viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [[viewerItem('viewport.maxPxPerPosition'), 32], ...config] as never,
  });
  viewers.push(viewer);
  await viewer.load(root);
  viewer.resize(800);
  await viewer.settle();
  frames.flush();
  viewer.draw();
  return viewer;
}

/** The navigator drawn over the viewer's last frame, its wants answered. */
async function drawNav(viewer: SeqViewer, id = 'strip') {
  const nav = viewer.navigatorModel(id)!;
  nav.draw();
  await viewer.settle();
  frames!.flush();
  viewer.draw();
  const f = nav.draw()!;
  return { nav, f };
}

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const hit = (a: Box, b: Box) => a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1;

describe('placeReadout', () => {
  const room = [600, 40] as const;
  const win = (x0: number, x1: number): Box => ({ x0: x0 - GripReach, x1: x1 + GripReach, y0: 0, y1: 40 });

  it('rests at the top right where nothing is under it', () => {
    expect(placeReadout(120, 20, ...room, [win(10, 60)])).toEqual([476, 4]);
  });

  it('a window under the corner: the top left, then beside it, then above', () => {
    expect(placeReadout(120, 20, ...room, [win(500, 590)])).toEqual([4, 4]);
    expect(placeReadout(120, 20, ...room, [win(0, 600)])[1]).toBeLessThan(0);
    // both corners taken, room between: beside the window on its right side
    expect(placeReadout(100, 20, ...room, [win(10, 50), win(520, 590)])[0]).toBeGreaterThan(50);
  });

  it("the pointer's box is avoided as a window is", () => {
    const p = pointerBox(540, 20)!;
    const [x, y] = placeReadout(120, 20, ...room, [], p);
    expect(hit({ x0: x, x1: x + 120, y0: y, y1: y + 20 }, p)).toBe(false);
    expect(p.x1 - p.x0).toBe(2 * PointerReach);
  });

  it('a pointer away has no box', () => {
    expect(pointerBox(Number.NaN, Number.NaN)).toBeUndefined();
  });
});

describe('a strip’s readout during a drag', () => {
  it('its box touches neither the window, its grips nor the pointer, at any place of the window', async () => {
    const v = await make(spec());
    for (const [a, b] of [
      [10, 40],
      [100, 150],
      [330, 393],
      [0, 393],
    ]) {
      v.view.setRange(P, a, b);
      v.draw();
      const { nav, f } = await drawNav(v);
      const w = f.windows[0];
      const mid = (w.dx0 + w.dx1) / 2;
      expect(nav.down(mid, 10)).toBe(true);
      nav.move(mid + 1, 10);
      const x = nav.width * 0.5;
      for (const [rw, rh] of [
        [120, 20],
        [200, 20],
      ]) {
        const [left, top] = nav.readoutAt(rw, rh);
        const box: Box = { x0: left, x1: left + rw, y0: top, y1: top + rh };
        const win: Box = { x0: w.dx0 - GripReach, x1: w.dx1 + GripReach, y0: 0, y1: nav.height };
        expect(hit(box, win), `window ${a}-${b}`).toBe(false);
        expect(hit(box, pointerBox(mid + 1, 10)!), `pointer ${a}-${b}`).toBe(false);
        expect(x).toBeGreaterThan(0);
      }
      nav.up();
    }
  });
});
