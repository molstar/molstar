/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Navigators in the viewer (_spec/README.md "Models", NavigatorModel; "Behaviors", Navigator): two consumers, a
// navigator and the view it moves, headless on a fake clock: drag the window and the view moves, move the view and the
// window follows; a click centres (Alt+Left returns), a sweep zooms, a double-click fits or zooms, Esc cancels; a far
// drag pushes one step, a nudge none; a follower's window leads; a scroll group's one window; two independent leaders'
// windows on one strip; idle on and off with no relayout; keys only; a removed bound section; the ruler rail on the
// canvas (dragged at its scale, grown, a click outside centres, idle no target); the export with `navigators`.

import type { Node } from '@molstar/core/repr/tree';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import { railUnitAt } from '@molstar/seqx-core/render';
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
    const [s, e] = range(v);
    expect((s + e) / 2).toBeCloseTo(300, 6);
    // the click's step waits for a double-click: its time runs out, the step is pushed
    expect(v.view.canBack).toBe(false);
    frames!.timers.fire();
    expect(v.view.canBack).toBe(true);
    expect(nav.key('ArrowLeft', { altKey: true })).toBe(true);
    expect(range(v)).toEqual([0, 40]);
  });

  it('a double-click: fit, then fitted a window zoomPx wide; one step with its first click', async () => {
    const v = await make(spec());
    v.view.setRange(P, 0, 40);
    v.draw();
    const { nav } = await drawNav(v);
    const x = (300 / 393) * nav.width;
    nav.down(x, 10);
    nav.up();
    nav.doubleClick(x, 10);
    expect(range(v)).toEqual([0, 393]);
    expect(v.view.history.value.back).toBe(1);
    v.draw();
    nav.draw();
    nav.doubleClick(x, 10);
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

  it('idle when the whole axis is shown at the zoom limit, live again with no relayout', async () => {
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
    // grown while its drag runs
    v.input.move(mid + 30 * k, y + 80);
    expect(v.input.dragHit?.navigator).toBe('rail');
    v.input.up();
    expect(range(v)[0]).toBeCloseTo(130, 6);
    expect(v.input.cursor).toBeUndefined();
  });

  it('a click outside the window centres there; idle no target', async () => {
    const { v, rail: r } = await rail();
    const x = r.x + railFraction(r, 300);
    v.input.down(x, r.bottom - 2);
    v.input.up();
    const [s, e] = range(v);
    expect((s + e) / 2).toBeCloseTo(300, 0);
    // the click after the press is the rail's: no selection
    v.input.click(x, r.bottom - 2);
    expect(v.selection.value.value.shown.parts.length).toBe(0);
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

  it('a declared navigator over the section replaces its rail', async () => {
    const v = await make(spec({ ruler: true }));
    expect(v.frame!.layout.sections[0].rail).toBeUndefined();
  });

  it('the viewer’s export draws the rail idle (the canvas’s geometry); with `navigators: include` the strip too', async () => {
    const v = await make(spec({ ruler: true }));
    const plain = await v.renderSvg();
    expect(plain.svg).not.toContain('data-navigator');
    const p = v.renderSvg({ navigators: 'include' });
    frames!.flush();
    const withNav = await p;
    expect(withNav.svg).toContain('data-navigator="strip"');
    expect(withNav.svg).toContain('<title>Window: p:');
  });
});

function railFraction(r: { width: number; total: number }, u: number) {
  return (u / r.total) * r.width;
}

// a unit at the rail's scale: what the input reads
void railUnitAt;
void SeqRuntime;
