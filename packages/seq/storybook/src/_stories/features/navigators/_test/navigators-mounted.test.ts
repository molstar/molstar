// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigators' stories in a DOM (jsdom) on a fake clock, composed as the stories compose them: a viewer mounted into a
// target with the default spec's UI (its toolbar, its canvas, the navigator parts placed around it), its sources
// answered with the story's own pacing (`registry()`, no delay given), the parts drawn when the platform's frames run,
// before the sources answer as in a browser. What each part's base canvas was painted with is recorded: chr17's bands
// are painted; the plasmid's strip draws its axis; and the parts are placed beside the view, never over it: in the
// document's flow around the canvas (jsdom lays nothing out: the boxes are read from the styles the mount sets).

import { GlobinsViews, P53Views, PlasmidViews, type ExampleView } from '@molstar/seqx-examples';
import type { NavigatorRenderer } from '@molstar/seqx-core/render/navigator';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { loadLegendEyes } from '@molstar/seqx-viewer/mount-legend';
import { afterAll, afterEach, beforeEach, describe, expect, it, vi } from 'vitest';

/** A canvas's 2D context that takes every call and records the fills of `fillRect` (their `fillStyle`). */
interface Recorded {
  readonly fills: string[];
}

function recordingContext(canvas: HTMLCanvasElement): unknown {
  const fills: string[] = [];
  const target: Record<string | symbol, unknown> = {
    canvas,
    fills,
    fillStyle: '',
    measureText: (t: string) => ({ width: 6 * t.length, actualBoundingBoxAscent: 8, actualBoundingBoxDescent: 2 }),
    getImageData: (_x: number, _y: number, w: number, h: number) => ({
      data: new Uint8ClampedArray(4 * w * h),
      width: w,
      height: h,
    }),
    createImageData: (w: number, h: number) => ({ data: new Uint8ClampedArray(4 * w * h), width: w, height: h }),
    getTransform: () => ({ a: 1, b: 0, c: 0, d: 1, e: 0, f: 0 }),
    createLinearGradient: () => ({ addColorStop() {} }),
    createRadialGradient: () => ({ addColorStop() {} }),
    createPattern: () => ({}),
    fillRect: () => fills.push(String(target.fillStyle)),
    fill: () => fills.push(String(target.fillStyle)),
  };
  return new Proxy(target, {
    get: (t, k) => (k in t ? t[k] : () => {}),
    set: (t, k, v) => {
      t[k] = v;
      return true;
    },
  });
}

const raf = DomHost.requestAnimationFrame;
let viewer: SeqViewer | undefined;
let target: HTMLElement | undefined;

beforeEach(() => {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as unknown as { __ctx?: unknown }).__ctx ??= recordingContext(this)) as never;
  });
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(456);
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  target?.remove();
  target = undefined;
  DomHost.requestAnimationFrame = raf;
  vi.useRealTimers();
  vi.restoreAllMocks();
});

/**
 * The story's viewer: the default spec's UI with the toolbar, the view's registry at its own pacing, its fixtures read a
 * while after the first frames (as a page fetches them), loaded, the clock run.
 */
async function story(view: ExampleView, gene = '') {
  target = document.createElement('div');
  document.body.append(target);
  const base = DefaultSeqViewerSpec();
  viewer = createSeqViewer(
    {
      ...base,
      ui: { ...base.ui, toolbar: mountToolbar },
      registry: [...(base.registry ?? []), ...(view.registry?.() ?? [])],
      // a fixture read as a page fetches it: a while after the first frames
      ports: {
        assets: () => ({
          reader: {
            readText: () => new Promise((done) => setTimeout(() => done(gene), 500)),
            readBytes: () => new Promise((done) => setTimeout(() => done(new TextEncoder().encode(gene)), 500)),
          },
        }),
      },
    },
    target,
  );
  void viewer.load(view.build());
  await vi.advanceTimersByTimeAsync(8000);
  return viewer;
}

/** A part's base canvas's fills, as painted. */
function painted(id: string): string[] {
  const part = target!.querySelector<HTMLElement>(`[data-navigator="${id}"]`)!;
  const canvas = part.querySelector('canvas')!;
  return ((canvas.getContext('2d') as unknown as Recorded).fills ?? []).slice();
}

// The legend's eyes are a lazy chunk the first mount in a worker starts: the file waits for it, or it lands after the file ends.
afterAll(() => loadLegendEyes());

describe('the navigators’ stories, mounted', () => {
  it('chr17’s ideogram paints its bands, at a device pixel ratio of 2 in a wide pane too', async () => {
    // 684 CSS px at 2 device pixels each: a device pixel is finer than the bands' coarse level
    vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(684);
    vi.stubGlobal('devicePixelRatio', 2);
    await story(P53Views.find((v) => v.id === 'gene-region-ideogram')!, regionSeqText());
    const fills = new Set(painted('chr17'));
    // the background and the giemsa grays (gneg to gpos100) and the centromere: many colors, not one flat one
    expect(fills.size).toBeGreaterThan(4);
    expect(viewer!.navigatorModel('chr17')!.issues.map((i) => i.code)).not.toContain('seq.render.navigator.layer');
    // and the store refused no want
    expect(viewer!.store.store.issues.value).toEqual([]);
  });

  it('the plasmid’s strip draws its axis', async () => {
    await story(PlasmidViews.find((v) => v.id === 'detail-strip')!);
    const m = viewer!.navigatorModel('strip')!;
    const picture = (m as unknown as { renderer?: NavigatorRenderer }).renderer?.pictureSection;
    expect(picture?.tracks.length ?? 0).toBeGreaterThan(0);
    expect(new Set(painted('strip')).size).toBeGreaterThan(1);
  });

  it('the parts beside the view, never over it: a strip above the canvas, below the toolbar; a side part at its width', {
    tags: ['medium'],
  }, async () => {
    await story(GlobinsViews.find((v) => v.id === 'scale-navigators')!);
    const toolbar = target!.querySelector<HTMLElement>('.seqx-toolbar')!;
    const canvas = target!.querySelector<HTMLElement>('[role="application"]')!;
    const overview = target!.querySelector<HTMLElement>('[data-navigator="overview"]')!;
    const mini = target!.querySelector<HTMLElement>('[data-navigator="miniature"]')!;
    const follows = (a: Node, b: Node) => (a.compareDocumentPosition(b) & Node.DOCUMENT_POSITION_FOLLOWING) !== 0;
    // in the document's flow: the toolbar, the miniature, then the canvas and the overview beside it in a row
    expect(follows(toolbar, mini)).toBe(true);
    expect(follows(mini, canvas)).toBe(true);
    expect(follows(canvas, overview)).toBe(true);
    for (const p of [mini, overview]) {
      for (let e: HTMLElement | null = p; e && e !== target; e = e.parentElement)
        expect(getComputedStyle(e).position, e.outerHTML.slice(0, 60)).not.toMatch(/absolute|fixed/);
    }
    // each at its form's box (a miniature 160 × 96, an overview 96 wide), never the host's width
    expect([mini.style.width, mini.style.height]).toEqual(['160px', '96px']);
    expect(overview.style.width).toBe('96px');
    // the canvas's column shrinks beside the overview: the viewer is laid out at the column's width, not the row's
    const middle = canvas.parentElement!;
    expect(middle.style.minWidth).toBe('0px');
    expect(viewer!.navigatorModel('overview')!.width).toBe(96);
    expect(viewer!.navigatorModel('miniature')!.width).toBe(160);
  });
});
