/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Navigators in the viewer (_spec/README.md "Models", NavigatorModel; "Behaviors", Navigator): two consumers, a
// navigator and the view it moves, headless on a fake clock, in the platform's event order (a click after a press that
// never moved, a double-click's second click before its `dblclick`): drag the window and the view moves, move the view
// and the window follows; a click centres (Alt+Left returns), a sweep zooms, a double-click fits or zooms in one step,
// Esc cancels; a far drag pushes one step, a nudge none; a follower's window leads; a scroll group's one window; two
// independent leaders' windows on one strip, the press focusing its slider; idle on and off; keys only; a removed
// bound section; the overview and the miniature over an alignment's rows (drags, keys, clicks, a far rows drag's
// step); the ruler rail on the canvas (dragged at its scale, grown, a click outside centres, a grip's click does not,
// a double-click fits, idle no target, touch, the grab-zone split on an alignment's ruler); the export.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, miniatureNavigator, overviewNavigator, stripNavigator } from '@molstar/seqx-core/lang';
import type { DrawLists, NavigatorRenderer } from '@molstar/seqx-core/render';
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

const Letters =
  'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD';

const P = '/views/v/protein';
const Twin = '/views/v/twin';

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

const range = (v: SeqViewer, s = P) => {
  const x = v.view.get(s)!;
  return [x.start, x.end] as const;
};

describe('a strip and the view it moves', () => {
  it('drag the window: the view moves with the hand; move the view: the window follows', async () => {
    const v = await make(spec());
    v.view.setRange(P, 100, 150);
    v.draw();
    const { nav, f } = await drawNav(v);
    const w = f.windows[0];
    const mid = (w.dx0 + w.dx1) / 2;
    expect(nav.down(mid, 10)).toBe(true);
    const px = nav.width / 393;
    nav.move(mid + 40 * px, 10);
    nav.up();
    const [s, e] = range(v);
    expect(s).toBeCloseTo(140, 6);
    expect(e).toBeCloseTo(190, 6);
    // a nudge pushes no step
    expect(v.view.canBack).toBe(false);
    // the view moves: the window follows
    v.view.setRange(P, 300, 350);
    v.draw();
    const g = nav.draw()!.windows[0];
    expect(g.x0).toBeCloseTo((300 / 393) * nav.width, 6);
  });

  it('a click outside centres there (a jump, animated); Alt+Left returns', async () => {
    const v = await make(spec());
    v.view.setRange(P, 0, 40);
    v.draw();
    const { nav } = await drawNav(v);
    const x = (300 / 393) * nav.width;
    expect(nav.hit(x, 10)?.part).toBe('outside');
    expect(nav.down(x, 10)).toBe(true);
    nav.up();
    // the platform's click after the press: a jump, its own step
    expect(nav.click(x, 10, 1)).toBe(true);
    const [s, e] = range(v);
    expect((s + e) / 2).toBeCloseTo(300, 6);
    expect(v.view.history.value.back).toBe(1);
    expect(nav.key('ArrowLeft', { altKey: true })).toBe(true);
    expect(range(v)).toEqual([0, 40]);
  });

  it('a double-click: fit, then fitted a window zoomPx wide; one step with its first click', async () => {
    const v = await make(spec());
    v.view.setRange(P, 0, 40);
    v.draw();
    const { nav } = await drawNav(v);
    const x = (300 / 393) * nav.width;
    // the platform's order: down, up, click, down, up, its second click, dblclick
    const double = () => {
      for (const detail of [1, 2]) {
        nav.down(x, 10);
        nav.up();
        nav.click(x, 10, detail);
      }
      nav.doubleClick(x, 10);
    };
    double();
    expect(range(v)).toEqual([0, 393]);
    expect(v.view.history.value.back).toBe(1);
    v.view.back();
    expect(range(v)).toEqual([0, 40]);
    v.view.forward();
    v.draw();
    nav.draw();
    double();
    const [s, e] = range(v);
    expect(e - s).toBeCloseTo((48 / nav.width) * 393, 6);
    expect((s + e) / 2).toBeCloseTo(300, 0);
    expect(v.view.history.value.back).toBe(2);
  });

  it('a sweep from outside zooms to what it swept on release (a jump); Esc during one cancels it', async () => {
    const v = await make(spec());
    v.view.setRange(P, 0, 40);
    v.draw();
    const { nav } = await drawNav(v);
    const at = (u: number) => (u / 393) * nav.width;
    nav.down(at(200), 10);
    nav.move(at(250), 10);
    nav.up();
    expect(range(v)[0]).toBeCloseTo(200, 0);
    expect(range(v)[1]).toBeCloseTo(250, 0);
    expect(v.view.canBack).toBe(true);
    v.draw();
    nav.draw();
    nav.down(at(300), 10);
    nav.move(at(350), 10);
    expect(nav.cancel()).toBe(true);
    nav.up();
    expect(range(v)[0]).toBeCloseTo(200, 0);
  });

  it('a far drag pushes the viewports at its press as one step', async () => {
    const v = await make(spec());
    v.view.setRange(P, 0, 20);
    v.draw();
    const before = range(v);
    const { nav, f } = await drawNav(v);
    const w = f.windows[0];
    const mid = (w.dx0 + w.dx1) / 2;
    nav.down(mid, 10);
    nav.move(mid + (200 / 393) * nav.width, 10);
    nav.up();
    expect(v.view.history.value.back).toBe(1);
    v.view.back();
    expect(range(v)).toEqual(before);
  });

  it('a floored window’s end dragged inward keeps the view; outward it zooms out from the drawn end', async () => {
    const v = await make(spec(), [[viewerItem('viewport.maxPxPerPosition'), 1000]]);
    v.view.setRange(P, 100, 101);
    v.draw();
    const { nav, f } = await drawNav(v);
    const w = f.windows[0];
    expect(w.dx1 - w.dx0).toBeCloseTo(16, 6);
    nav.down(w.dx1 + 1, 10);
    nav.move(w.dx1 - 4, 10);
    expect(range(v)).toEqual([100, 101]);
    nav.move(w.dx1 + 10, 10);
    const drawnEnd = (w.dx1 / nav.width) * 393;
    expect(range(v)[1]).toBeGreaterThan(drawnEnd);
    nav.up();
  });

  it('idle when the whole axis is shown at the zoom limit, live again when a zoom becomes possible', async () => {
    const short = createSeqSpec();
    const p = short.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
    short.views.view('v', (vv) => {
      const s = vv.section('protein', { of: p });
      s.track('letters').layer('letters', { type: 'letters' });
      vv.navigator('strip', stripNavigator(s));
    });
    const v = await make(short.build());
    const { nav, f } = await drawNav(v);
    expect(nav.values.value[0].disabled).toBe(true);
    expect(nav.down(400, 10)).toBe(false);
    // idle is drawn quiet, never hidden (a part that comes and goes moves the host's layout)
    expect(f.windows[0].shown).toBe(true);
    // a narrower canvas: a zoom becomes possible
    v.resize(200);
    v.draw();
    nav.resize(200, 20);
    nav.draw();
    // the sliders say it when the view rests
    frames!.timers.fire();
    frames!.flush();
    const g = nav.draw()!;
    expect(nav.values.value[0].disabled).toBe(false);
    expect(g.windows[0].shown).toBe(true);
    expect(nav.down((g.windows[0].dx0 + g.windows[0].dx1) / 2, 10)).toBe(true);
    nav.up();
  });

  it('keys only: arrows pan, Shift+arrows resize, Home and End, Enter centres on the readout', async () => {
    const v = await make(spec());
    v.view.setRange(P, 100, 140);
    v.draw();
    const { nav } = await drawNav(v);
    expect(nav.key('ArrowRight')).toBe(true);
    expect(range(v)[0]).toBeCloseTo(104, 6);
    expect(nav.key('ArrowLeft', { shiftKey: true })).toBe(true);
    expect(range(v)[1] - range(v)[0]).toBeCloseTo(40 / 1.25, 6);
    expect(nav.key('Home')).toBe(true);
    expect(range(v)[0]).toBe(0);
    expect(nav.key('End')).toBe(true);
    expect(range(v)[1]).toBe(393);
    nav.move((200 / 393) * nav.width, 10);
    expect(nav.readout.value).toBe('p:201');
    expect(nav.key('Enter')).toBe(true);
    expect((range(v)[0] + range(v)[1]) / 2).toBeCloseTo(200, 0);
    // the slider's value said on settle
    v.draw();
    nav.draw();
    expect(nav.values.value[0].text).toMatch(/^p:\d+–\d+ of p$/);
  });

  it('two independent leaders on one strip: a window each, each drag moving its own section', async () => {
    // two sections on one system that do not scroll together
    const v = await make(spec({ twin: true, same: true }), [[viewerItem('viewport.groupEqualAxes'), false]]);
    v.view.setRange(P, 0, 50);
    v.view.setRange(Twin, 200, 260);
    v.draw();
    const { nav, f } = await drawNav(v);
    expect(f.windows.map((w) => w.section)).toEqual([P, Twin]);
    const tw = f.windows[1];
    const mid = (tw.dx0 + tw.dx1) / 2;
    expect(nav.hit(mid, 10)?.section).toBe(Twin);
    nav.down(mid, 10);
    nav.move(mid + (20 / 393) * nav.width, 10);
    nav.up();
    expect(range(v, Twin)[0]).toBeCloseTo(220, 6);
    expect(range(v)).toEqual([0, 50]);
  });

  it('a follower’s window leads: dragged, its leader moves', async () => {
    const v = await make(spec({ twin: true, follow: true }));
    v.view.setRange(P, 0, 50);
    v.draw();
    const { nav, f } = await drawNav(v);
    const tw = f.windows[1];
    expect(tw.shown).toBe(true);
    const mid = (tw.dx0 + tw.dx1) / 2;
    nav.down(mid, 10);
    nav.move(mid + (30 / 393) * nav.width, 10);
    nav.up();
    expect(range(v)[0]).toBeGreaterThan(20);
  });

  it('a scroll group shares one range: one window', async () => {
    const v = await make(spec({ twin: true, group: true }));
    const { f } = await drawNav(v);
    expect(f.windows.length).toBe(1);
    expect(v.view.groupOf(P)).toContain(Twin);
  });

  it('a removed bound section: said, its window not drawn', async () => {
    const v = await make(spec());
    const { nav } = await drawNav(v);
    const root = v.spec.shown.value;
    const views = root.children!.find((c) => c.kind === 'views')!;
    const view = views.children![0];
    const edited: Node = {
      ...root,
      children: root.children!.map((c) =>
        c !== views
          ? c
          : {
              ...views,
              children: [
                {
                  ...view,
                  children: view.children!.map((n) =>
                    n.kind === 'navigators'
                      ? { ...n, children: [{ ...n.children![0], params: { of: '/views/v/gone' } }] }
                      : n,
                  ),
                },
              ],
            },
      ),
    };
    await v.load(edited);
    await v.settle();
    frames!.flush();
    v.draw();
    const f = nav.draw();
    expect(f?.windows.length ?? 0).toBe(0);
    expect(nav.issues.map((i) => i.code)).toContain('seq.render.navigator.of');
  });

  it('a host’s navigator, never saved: the same model', async () => {
    const v = await make(spec({ strip: false }));
    const m = v.navigator({ id: 'host', of: P, content: 'axis' });
    v.view.setRange(P, 10, 30);
    v.draw();
    m.draw();
    await v.settle();
    frames!.flush();
    v.draw();
    const f = m.draw()!;
    expect(f.windows[0].shown).toBe(true);
    expect(v.navigatorModel('host')).toBe(m);
    expect(JSON.stringify(v.spec.shown.value)).not.toContain('host');
  });

  it('the wheel over the navigator zooms about the window’s centre, not the pointer', async () => {
    const v = await make(spec());
    v.view.setRange(P, 100, 140);
    v.draw();
    const { nav } = await drawNav(v);
    expect(nav.wheel(700, 10, 0, -100, { mod: true })).toBe(true);
    frames!.flush();
    const [s, e] = range(v);
    expect((s + e) / 2).toBeCloseTo(120, 0);
    expect(e - s).toBeLessThan(40);
  });
});

describe('a strip’s sliders and keys', () => {
  it('a press on a window focuses its slider; PageDown and PageUp rove; the keys act on the focused window', async () => {
    const v = await make(spec({ twin: true, same: true }), [[viewerItem('viewport.groupEqualAxes'), false]]);
    v.view.setRange(P, 0, 50);
    v.view.setRange(Twin, 200, 260);
    v.draw();
    const { nav, f } = await drawNav(v);
    expect(nav.focused.value).toBe(0);
    const tw = f.windows[1];
    const mid = (tw.dx0 + tw.dx1) / 2;
    nav.down(mid, 10);
    nav.up();
    nav.click(mid, 10, 1);
    expect(nav.focused.value).toBe(1);
    // the keys move the twin, the protein stays
    expect(nav.key('ArrowRight')).toBe(true);
    expect(range(v, Twin)[0]).toBeGreaterThan(200);
    expect(range(v)).toEqual([0, 50]);
    expect(nav.key('PageDown')).toBe(true);
    expect(nav.focused.value).toBe(0);
    expect(nav.key('PageUp')).toBe(true);
    expect(nav.focused.value).toBe(1);
  });

  it('Alt+Left and Alt+Right walk the history while the window is idle', async () => {
    const short = createSeqSpec();
    const p = short.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
    const q = short.data.sequence('q', { alphabet: 'protein', letters: Letters });
    short.views.view('v', (vv) => {
      const s = vv.section('protein', { of: p });
      s.track('letters').layer('letters', { type: 'letters' });
      vv.section('other', { of: q }).track('letters').layer('letters', { type: 'letters' });
      vv.navigator('strip', stripNavigator(s));
    });
    const v = await make(short.build(), [[viewerItem('view.jumpMs'), 0]]);
    const O = '/views/v/other';
    // a step of the view's history, made on the other section
    v.view.jump(() => {
      v.view.setRange(O, 100, 140);
      return true;
    });
    v.draw();
    const { nav } = await drawNav(v);
    expect(nav.values.value[0].disabled).toBe(true);
    expect(nav.key('ArrowLeft', { altKey: true })).toBe(true);
    expect(range(v, O)).not.toEqual([100, 140]);
    expect(nav.key('ArrowRight', { altKey: true })).toBe(true);
    expect(range(v, O)).toEqual([100, 140]);
    // the window's own keys do nothing while it is idle
    expect(nav.key('ArrowRight')).toBe(false);
  });

  it('the picture shares the view’s image cache', async () => {
    const v = await make(spec());
    const { nav } = await drawNav(v);
    const r = (nav as unknown as { renderer: NavigatorRenderer }).renderer;
    expect(r.options.images).toBe(v.renderer!.images);
    expect(r.options.domains).toBe(v.renderer!.layoutState);
  });

  it('a drag across a cyclic origin: the view follows over it', async () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: Letters, cyclic: true });
    b.views.view('v', (vv) => {
      const s = vv.section('protein', { of: p });
      s.track('letters').layer('letters', { type: 'letters' });
      vv.navigator('strip', stripNavigator(s));
    });
    const v = await make(b.build());
    v.view.setRange(P, 360, 380);
    v.draw();
    const { nav, f } = await drawNav(v);
    const w = f.windows[0];
    const mid = (w.dx0 + w.dx1) / 2;
    const px = nav.width / 393;
    nav.down(mid, 10);
    nav.move(mid + 10 * px, 10);
    nav.move(mid + 25 * px, 10);
    nav.up();
    const [s, e] = range(v);
    // over the origin: its end past the axis's length
    expect(s).toBeCloseTo(385, 6);
    expect(e).toBeGreaterThan(393);
    v.draw();
    const g = nav.draw()!.windows[0];
    expect(g.wrap).toBe(true);
  });
});

/** An alignment of `n` rows 40 columns wide in a short band, an overview or a miniature over it. */
function rowsSpec(n: number, form: 'overview' | 'miniature') {
  const b = createSeqSpec({ title: 'rows' });
  const label = Array.from({ length: n }, (_, i) => `r${i}`);
  const letters = Array.from({ length: n }, (_, i) => Letters.slice(i % 40, (i % 40) + 40));
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label, letters } });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    s.track('ruler').layer('ruler', { type: 'ruler' });
    s.track('rows', { label: '', height: 60 }).layer('r', { type: 'rows', data: aln });
    v.navigator('nav', form === 'overview' ? overviewNavigator(s) : miniatureNavigator(s));
  });
  return b.build();
}

const S = '/views/v/s';
const RowsTrack = '/views/v/s/rows';
const top = (v: SeqViewer) => v.view.band(RowsTrack)?.bandTop ?? 0;

describe('an overview and a miniature over an alignment’s rows', () => {
  it('the overview: a drag moves the rows, clamped at the last; a far drag is one step; keys; a click centres', async () => {
    const v = await make(rowsSpec(60, 'overview'), [[viewerItem('view.jumpMs'), 0]]);
    v.view.setBand(RowsTrack, { bandTop: 30 });
    v.draw();
    const { nav, f } = await drawNav(v, 'nav');
    const w = f.windows[0];
    const shown = w.r1 - w.r0;
    expect(shown).toBeGreaterThan(0);
    expect(shown).toBeLessThan(30);
    // the rows fill the column
    expect(w.hy - w.ly).toBeCloseTo(nav.height, 0);
    const x = (w.dx0 + w.dx1) / 2;
    const y = (w.dy0 + w.dy1) / 2;
    expect(nav.down(x, y)).toBe(true);
    // far past the end: the rows stop at the last ones
    nav.move(x, nav.height + 400);
    expect(top(v)).toBeCloseTo(60 - shown, 6);
    nav.up();
    // farther than the rows in view: one step, back to where the press found them
    expect(v.view.history.value.back).toBe(1);
    v.view.back();
    expect(top(v)).toBe(30);
    v.draw();
    nav.draw();
    // keys: Down scrolls, clamped at the end
    expect(nav.key('ArrowDown')).toBe(true);
    expect(top(v)).toBeGreaterThan(0);
    for (let i = 0; i < 100; i++) nav.key('ArrowDown');
    expect(top(v)).toBeCloseTo(60 - shown, 6);
    // a click below the window centres the rows there (a jump), clamped
    v.view.setBand(RowsTrack, { bandTop: 0 });
    v.draw();
    const g = nav.draw()!.windows[0];
    const yy = g.hy - 2;
    nav.down(x, yy);
    nav.up();
    nav.click(x, yy, 1);
    expect(top(v)).toBeCloseTo(60 - shown, 0);
  });

  it('the miniature: its body drags both axes; a click outside centres both', async () => {
    const v = await make(rowsSpec(60, 'miniature'), [[viewerItem('view.jumpMs'), 0]]);
    v.view.setRange(S, 0, 10);
    v.draw();
    const { nav, f } = await drawNav(v, 'nav');
    const w = f.windows[0];
    expect(w.axes).toBe(3);
    const x = (w.dx0 + w.dx1) / 2;
    const y = (w.dy0 + w.dy1) / 2;
    expect(nav.hit(x, y)?.part).toBe('body');
    expect(nav.down(x, y)).toBe(true);
    const cols = nav.width / 40;
    const rows = (w.hy - w.ly) / 60;
    nav.move(x + 10 * cols, y + 10 * rows);
    nav.up();
    expect(range(v, S)[0]).toBeCloseTo(10, 0);
    expect(top(v)).toBeCloseTo(10, 0);
    v.draw();
    const g = nav.draw()!.windows[0];
    // a click outside both: the columns and the rows centre there
    const cx = g.lo + 2;
    const cy = g.hy - 2;
    nav.down(cx, cy);
    nav.up();
    nav.click(cx, cy, 1);
    expect(range(v, S)[0]).toBeCloseTo(0, 0);
    expect(top(v)).toBeGreaterThan(10);
  });

  it('on an alignment’s ruler the rail takes the lower share; a drag above it selects columns', async () => {
    const b = createSeqSpec({ title: 'aln ruler' });
    const label = Array.from({ length: 6 }, (_, i) => `r${i}`);
    const letters = label.map(() => Letters.slice(0, 393));
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label, letters } });
    b.views.view('v', (vv) => {
      const s = vv.section('s', { of: aln });
      s.track('ruler').layer('ruler', { type: 'ruler' });
      s.track('rows', { label: '', height: 60 }).layer('r', { type: 'rows', data: aln });
    });
    const v = await make(b.build());
    v.view.setRange(S, 100, 150);
    v.draw();
    const sec = v.frame!.layout.sections[0];
    const r = sec.rail!;
    const ruler = sec.tracks[0].box;
    const share = v.context.renderConfig.navigator.rail.grabShare;
    expect(r.top).toBeCloseTo(ruler.y + ruler.height * (1 - share), 6);
    const x = r.x + railFraction(r, 125);
    expect(v.hit(x, r.bottom - 1)?.navigator).toBe('rail');
    expect(v.hit(x, ruler.y + 1)?.navigator).toBeUndefined();
    // above the grab zone a drag selects columns, the view stays
    expect(v.input.down(x, ruler.y + 1)).toBe(true);
    v.input.move(x + railFraction(r, 5), ruler.y + 1);
    v.input.up();
    expect(range(v, S)).toEqual([100, 150]);
    expect(v.selection.value.value.shown.parts.length).toBeGreaterThan(0);
  });
});

describe('the ruler rail on the canvas', () => {
  async function rail() {
    const v = await make(spec({ ruler: true, strip: false }));
    v.view.setRange(P, 100, 150);
    v.draw();
    const s = v.frame!.layout.sections[0];
    return { v, s, rail: s.rail! };
  }

  it('every linear section has one by default; a drag on it pans at the rail’s scale', async () => {
    const { v, rail: r } = await rail();
    expect(r.active).toBe(true);
    const k = r.width / r.total;
    const mid = r.x + 125 * k;
    const y = r.bottom - 2;
    expect(v.hit(mid, y)?.navigator).toBe('rail');
    expect(v.input.down(mid, y)).toBe(true);
    v.input.move(mid + 30 * k, y + 80);
    expect(v.input.dragHit?.navigator).toBe('rail');
    v.input.up();
    expect(range(v)[0]).toBeCloseTo(130, 6);
    expect(v.input.cursor).toBeUndefined();
  });

  it('a click outside the window centres there (the click is the rail’s: no selection); idle no target', async () => {
    const { v, rail: r } = await rail();
    const x = r.x + railFraction(r, 300);
    v.input.down(x, r.bottom - 2);
    v.input.up();
    v.input.click(x, r.bottom - 2, { detail: 1 });
    const [s, e] = range(v);
    expect((s + e) / 2).toBeCloseTo(300, 0);
    expect(v.selection.value.value.shown.parts.length).toBe(0);
    // idle: the whole axis at the zoom limit; the ruler's own hit under the pointer, no window
    const short = createSeqSpec();
    const p = short.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
    short.views.view('v', (vv) => {
      const sec = vv.section('protein', { of: p });
      sec.track('ruler').layer('ruler', { type: 'ruler' });
      sec.track('letters').layer('letters', { type: 'letters' });
    });
    const w = await make(short.build());
    const idle = w.frame!.layout.sections[0].rail!;
    expect(idle.active).toBe(false);
    expect(w.hit(idle.x + idle.width / 2, idle.bottom - 2)?.target).not.toBe('window');
    // the press falls through to the ruler: no rail drag
    w.input.down(idle.x + idle.width / 2, idle.bottom - 2);
    expect(w.input.dragHit).toBeUndefined();
    w.input.up();
  });

  it('a click on a grip’s tab, not moved, does not centre: only a press outside the window does', async () => {
    const v = await make(spec({ ruler: true, strip: false }), [[viewerItem('viewport.maxPxPerPosition'), 1000]]);
    // a narrow window, floored: its grips' tabs outside its drawn ends
    v.view.setRange(P, 100, 102);
    v.draw();
    const r = v.frame!.layout.sections[0].rail!;
    const y = r.bottom - 2;
    const c = v.context.renderConfig.navigator.window;
    const mid = r.x + railFraction(r, 101);
    for (const x of [mid - c.minPx / 2 - 3, mid + c.minPx / 2 + 3]) {
      expect(v.hit(x, y)?.part).not.toBe('outside');
      v.input.down(x, y);
      v.input.up();
      expect(range(v)).toEqual([100, 102]);
      v.input.click(x, y, { detail: 1 });
      expect(range(v)).toEqual([100, 102]);
    }
  });

  it('a double-click outside the window: one step (down, up, click, down, up, its second click, dblclick)', async () => {
    const { v, rail: r } = await rail();
    const x = r.x + railFraction(r, 300);
    const y = r.bottom - 2;
    for (const detail of [1, 2]) {
      v.input.down(x, y);
      v.input.up();
      v.input.click(x, y, { detail });
    }
    // not fitted: the double-click fits (`navigator.doubleClick: toggle-fit`)
    expect(v.input.doubleClick(x, y)).toBe(true);
    expect(range(v)).toEqual([0, 393]);
    expect(v.view.history.value.back).toBe(1);
    v.view.back();
    expect(range(v)).toEqual([100, 150]);
  });

  it('grown while its drag runs (drawn so); by touch its grips hit wider; the cursor says grab', async () => {
    const { v, rail: r } = await rail();
    const k = r.width / r.total;
    const y = r.bottom - 2;
    const start = r.x + 100 * k;
    // a finger's hit: 12 px outside the start is its grip, a mouse's is outside
    v.input.move(start - 10, y, { pointerType: 'mouse' });
    expect(v.input.windowHit?.part).toBe('outside');
    v.input.move(start - 10, y, { pointerType: 'touch' });
    expect(v.input.windowHit?.part).toBe('start');
    expect(v.input.cursor).toBe('grab');
    const mid = r.x + 125 * k;
    v.input.down(mid, y, { pointerType: 'mouse' });
    v.input.move(mid + 30 * k, y + 80);
    expect(v.input.cursor).toBe('grabbing');
    v.draw();
    const grown = v.context.renderConfig.navigator.rail.grownPx;
    expect(readRects(v.frame!.lists, 'window').some((x) => x.h === grown)).toBe(true);
    v.input.up();
    v.draw();
    expect(readRects(v.frame!.lists, 'window').some((x) => x.h === grown)).toBe(false);
  });

  it('hidden by the host (`navigator.form: none`) or the reader (the corner button)', async () => {
    const v = await make(spec({ ruler: true, strip: false }), [[viewerItem('navigator.form'), 'none']]);
    expect(v.frame!.layout.sections[0].rail).toBeUndefined();
    const w = await make(spec({ ruler: true, strip: false }));
    expect(w.frame!.layout.sections[0].rail).toBeDefined();
    w.setRailsShown(false);
    w.draw();
    expect(w.frame!.layout.sections[0].rail).toBeUndefined();
  });

  it('a declared navigator over the section replaces its rail, and a host’s does; hidden, the rail is back', async () => {
    const v = await make(spec({ ruler: true }));
    expect(v.frame!.layout.sections[0].rail).toBeUndefined();
    const w = await make(spec({ ruler: true, strip: false }));
    expect(w.frame!.layout.sections[0].rail).toBeDefined();
    w.navigator({ id: 'host', of: P });
    frames!.flush();
    expect(w.frame!.layout.sections[0].rail).toBeUndefined();
    const root = v.spec.shown.value;
    const hidden = JSON.parse(JSON.stringify(root).replace('"axes":"positions"', '"axes":"positions","hidden":true'));
    await v.load(hidden);
    await v.settle();
    frames!.flush();
    v.draw();
    expect(v.frame!.layout.sections[0].rail).toBeDefined();
  });

  it('the viewer’s export leaves the rails out; with `navigators: include` the strip and the rails, plain', async () => {
    const v = await make(spec({ ruler: true, strip: false }));
    const plain = await v.renderSvg();
    expect(plain.svg).not.toContain('data-navigator');
    expect(plain.lists.passes.window.rects.n).toBe(0);
    const w = await make(spec({ ruler: true }));
    const none = await w.renderSvg();
    expect(none.svg).not.toContain('data-navigator');
    const p = w.renderSvg({ navigators: 'include' });
    frames!.flush();
    const withNav = await p;
    expect(withNav.svg).toContain('data-navigator="strip"');
    expect(withNav.svg).toContain('<title>Window: p:');
  });
});

/** The window pass's rects as drawn (x, height). */
function readRects(lists: DrawLists, pass: 'window') {
  const r = lists.passes[pass].rects;
  return Array.from({ length: r.n }, (_, i) => ({ x: r.x[i], h: r.h[i] }));
}

function railFraction(r: { width: number; total: number }, u: number) {
  return (u / r.total) * r.width;
}
