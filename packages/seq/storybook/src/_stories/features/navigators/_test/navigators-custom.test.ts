/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The custom navigators' stories, headless on a fake clock, composed as the stories compose them: the coverage strip (a
// host's picture registered with the default spec's registry: its values drawn, its readout, the SVG including it), and
// the host-drawn navigator (the viewer's `navigator(def)` model, the host's picture of its domains and windows, a drag
// through the model moving the view).

import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { CoverageView, HostSection, hostDrawnSpec, hostPicture } from '../compose.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

describe('the custom navigators’ stories, headless', () => {
  it('the coverage strip: the host’s picture drawn, its readout, the SVG including it', async () => {
    viewer = await openExample(CoverageView, {});
    const v = viewer;
    v.view.setRange('/views/main/protein', 100, 150);
    v.draw();
    const nav = v.navigatorModel('coverage')!;
    nav.draw();
    await vi.advanceTimersByTimeAsync(600);
    v.draw();
    const f = nav.draw()!;
    expect(nav.height).toBe(32);
    // the window is the view's range, placed through the model as a host would place its marks
    const [x0, x1] = nav.xRange(100, 150)!;
    expect(x0).toBeCloseTo(f.windows[0].x0, 6);
    expect(x1).toBeCloseTo(f.windows[0].x1, 6);
    // the picture's readout outside the window: the residue and its conservation
    nav.move(f.windows[0].dx1 + 200, 16);
    expect(nav.readout.value).toMatch(/^residue \d+ · conservation \d\.\d\d$/);
    expect(nav.issues).toEqual([]);
    const p = v.renderSvg({ navigators: 'include' });
    await vi.advanceTimersByTimeAsync(600);
    const r = await p;
    expect(r.svg).toContain('data-navigator="coverage"');
    expect(r.issues.filter((i) => i.code.startsWith('seq.render.navigator'))).toEqual([]);
  });

  it('the host-drawn navigator: the host’s blocks and window from the model, a drag through it pans', async () => {
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
    const v = (viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined }));
    const model = v.navigator({ id: 'host', of: HostSection, content: 'axis' });
    v.resize(720);
    void v.load(hostDrawnSpec());
    await vi.advanceTimersByTimeAsync(2000);
    v.view.setRange(HostSection, 100, 150);
    v.draw();
    model.resize(600, 36);
    model.draw();
    await vi.advanceTimersByTimeAsync(600);
    v.draw();
    model.draw();
    const pic = hostPicture(model);
    expect(pic.blocks.map((b) => b.label)).toEqual(['TAD', 'PRR', 'DBD', 'NLS', 'TET', 'CTD']);
    expect(pic.windows.length).toBe(1);
    const w = pic.windows[0];
    expect(w.x).toBeCloseTo((100 / 393) * 600, 0);
    expect(w.idle).toBe(false);
    // the host maps its pointer to the model's methods: a drag on the window pans the view
    const mid = w.x + w.width / 2;
    expect(model.down(mid, 18)).toBe(true);
    model.move(mid + (40 / 393) * 600, 18);
    model.up();
    expect(v.view.get(HostSection)!.start).toBeCloseTo(140, 0);
    // the declared rail it replaces is gone: a host's navigator binds the section
    expect(v.frame!.layout.sections[0].rail).toBeUndefined();
  });
});
