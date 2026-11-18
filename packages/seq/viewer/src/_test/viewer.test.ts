/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqContext } from '@molstar/seqx-core/context';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { CanvasBackend } from '../canvas/backend.js';
import { FakeCanvas } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { HoverProvider } from '../behaviors/providers.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { viewerItem } from '../config.js';
import { Gene, Protein, proteinOnGene } from './fixture.js';
import { installFrames, type Frames } from './frames.js';

// The viewer composed from the default spec, headless, with a fake frame host and a recording canvas.

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

// the fixture's four residues fill the width: zooming in needs a larger most pixels a position
const TinyConfig = [[viewerItem('viewport.maxPxPerPosition'), 1000]] as const;

/** The protein section's letters cell `i` of 4 at width 600: the data box starts at x = 128, 116 px a residue. */
const cellX = (i: number) => 128 + 116 * i + 58;
const ProteinY = 18;

async function make(options: { canvas?: boolean; width?: number } = {}) {
  frames = installFrames();
  const spec = { ...DefaultSeqViewerSpec(), ui: undefined, config: TinyConfig };
  const viewer = createSeqViewer(spec);
  viewers.push(viewer);
  const base = new FakeCanvas();
  const overlay = new FakeCanvas();
  const labels = new FakeCanvas();
  const backend = new CanvasBackend(base, overlay, labels);
  if (options.canvas) viewer.useBackend(backend);
  await viewer.load(proteinOnGene());
  viewer.resize(options.width ?? 600);
  await viewer.settle();
  frames.flush();
  return { viewer, base, overlay, labels, backend, frames };
}

describe('createSeqViewer', () => {
  it('composes the models, the bindings and the default behaviors from the spec, and shows a spec', async () => {
    const { viewer } = await make();
    expect([...viewer.behaviors.keys()]).toEqual(['hover', 'navigate', 'addressing']);
    expect(viewer.bindings.issues).toEqual([]);
    expect(viewer.commands.map((c) => c.id)).toContain('view.fit');
    expect(viewer.store.store.coverage(Protein)?.systems[0].length).toBe(4);
    expect(viewer.store.store.coverage(Gene)?.systems[0].length).toBe(24);
    const f = viewer.frame!;
    expect(f.layout.sections.map((s) => s.section.address)).toEqual(['/views/main/protein', '/views/main/gene']);
    expect(viewer.layoutIssues.value).toEqual([]);
    expect(viewer.height).toBeGreaterThan(0);
  });

  it('moves layoutIssues only when the issues change, not with each layout', async () => {
    const { viewer, frames } = await make();
    const before = viewer.layoutIssues.value;
    const layouts = viewer.counts.layouts;
    let told = 0;
    const off = viewer.layoutIssues.listen(() => told++);
    viewer.resize(640);
    frames.flush();
    expect(viewer.counts.layouts).toBeGreaterThan(layouts);
    expect(viewer.layoutIssues.value).toBe(before);
    expect(told).toBe(0);
    off();
  });

  it('lists a behavior once and an empty spec makes a viewer with nothing to do', () => {
    const spec = DefaultSeqViewerSpec();
    const twice = { ...spec, behaviors: [...spec.behaviors, spec.behaviors[0]] };
    expect(() => createSeqViewer(twice)).toThrow(/twice/);
    const empty = createSeqViewer({ behaviors: [] });
    expect(empty.behaviors.size).toBe(0);
    empty.retire();
  });

  it('takes config items from the spec', () => {
    const spec = { ...DefaultSeqViewerSpec(), config: [[viewerItem('viewport.maxPxPerPosition'), 99]] as const };
    const v = createSeqViewer(spec);
    expect(v.config.viewport.maxPxPerPosition).toBe(99);
    expect(v.config.hover.projectLinked).toBe(true);
    v.retire();
  });

  it('a hit under the pointer is the layout’s: a letters cell with its position', async () => {
    const { viewer } = await make();
    const h = viewer.hit(cellX(2), ProteinY);
    expect(h).toMatchObject({ system: Protein, position: 2 });
    expect(viewer.hit(2, 2)).toBeUndefined();
  });

  it('the SVG of the viewer equals the context’s of the same view, size and window', async () => {
    const { viewer } = await make();
    const mine = await viewer.renderSvg();
    const ctx = createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
    const theirs = await ctx.renderSvg(proteinOnGene(), {
      size: { width: viewer.width, height: viewer.height },
      ...(viewer.viewport ? { viewport: viewer.viewport } : {}),
    });
    expect(mine.settled).toBe(true);
    expect(mine.svg).toBe(theirs.svg);
    expect(mine.svg.length).toBeGreaterThan(500);
    ctx.retire();
  });
});

describe('the Canvas2D drawing', () => {
  it('draws the picture on the stacked canvases and sizes them for the device ratio', async () => {
    const { viewer, base, overlay, labels } = await make({ canvas: true });
    expect(labels.width).toBe(base.width);
    expect(base.width).toBe(viewer.width * viewer.dpr);
    expect(base.height).toBe(viewer.height * viewer.dpr);
    expect(overlay.width).toBe(base.width);
    expect(base.context.calls('fillRect').length).toBeGreaterThan(3);
    // the letters are data, on the base canvas; the track labels are above everything
    expect(base.context.calls('fillText').length).toBeGreaterThan(3);
    expect(labels.context.calls('fillText').length).toBeGreaterThan(0);
  });

  it('a hover redraws the overlay only', async () => {
    const { viewer, base, overlay, labels, frames } = await make({ canvas: true });
    base.context.take();
    labels.context.take();
    overlay.context.take();
    const layouts = viewer.counts.layouts;
    viewer.input.move(cellX(1), ProteinY);
    frames.flush();
    expect(viewer.hover.active).toBe(true);
    expect(base.context.take()).toEqual([]);
    expect(labels.context.take()).toEqual([]);
    expect(overlay.context.calls('clearRect').length).toBe(1);
    expect(overlay.context.calls('strokeRect').length).toBeGreaterThan(0);
    // no layout was made for it: only the emphasis passes
    expect(viewer.counts.layouts).toBe(layouts);
    // it moves to the next cell: the overlay again, still no base
    overlay.context.take();
    viewer.input.move(cellX(2), ProteinY);
    frames.flush();
    expect(base.context.take()).toEqual([]);
    expect(overlay.context.calls('strokeRect').length).toBeGreaterThan(0);
    // it leaves: the overlay is cleared
    overlay.context.take();
    viewer.input.leave();
    frames.flush();
    expect(viewer.hover.active).toBe(false);
    expect(overlay.context.calls('clearRect').length).toBe(1);
    expect(overlay.context.calls('strokeRect').length).toBe(0);
    expect(base.context.take()).toEqual([]);
    expect(viewer.counts.layouts).toBe(layouts);
  });
});

describe('input through the bindings', () => {
  it('a wheel with the zoom modifier zooms about the pointer, a plain wheel and Shift pan', async () => {
    const { viewer, frames } = await make();
    viewer.input.mac = false;
    const v0 = viewer.view.get('/views/main/protein')!;
    expect([v0.start, v0.end]).toEqual([0, 4]);
    // Ctrl+wheel up: in
    expect(viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: -100, ctrlKey: true, deltaMode: 1 })).toBe(true);
    frames.flush();
    const v1 = viewer.view.get('/views/main/protein')!;
    expect(v1.end - v1.start).toBeLessThan(4);
    // a plain wheel pans (the view is zoomed in, so there is room)
    const before = viewer.view.get('/views/main/protein')!.start;
    viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: 40 });
    frames.flush();
    expect(viewer.view.get('/views/main/protein')!.start).toBeGreaterThan(before);
    // Shift with a vertical wheel pans sideways
    const b2 = viewer.view.get('/views/main/protein')!.start;
    viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: -20, shiftKey: true });
    frames.flush();
    expect(viewer.view.get('/views/main/protein')!.start).toBeLessThan(b2);
  });

  it('a drag on empty space pans, one on data does not (the default preset)', async () => {
    const { viewer, frames } = await make();
    viewer.input.mac = false;
    viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: -300, ctrlKey: true, deltaMode: 1 });
    frames.flush();
    const start = viewer.view.get('/views/main/protein')!.start;
    // on a letter: no pan binding for a drag on a cell
    expect(viewer.input.down(cellX(1), ProteinY)).toBe(false);
    // the label column is empty space of the section
    expect(viewer.hit(40, ProteinY)).toBeUndefined();
    expect(viewer.input.down(40, ProteinY)).toBe(true);
    viewer.input.move(80, ProteinY);
    viewer.input.up();
    frames.flush();
    expect(viewer.view.get('/views/main/protein')!.start).toBeLessThan(start);
    // after the button went up a move does not pan
    const after = viewer.view.get('/views/main/protein')!.start;
    viewer.input.move(120, ProteinY);
    expect(viewer.view.get('/views/main/protein')!.start).toBe(after);
  });

  it('a trackpad pinch (a Ctrl+wheel with small deltas) is the pinch gesture: it zooms about the pointer', async () => {
    const { viewer, frames } = await make();
    viewer.input.mac = true;
    const v0 = viewer.view.get('/views/main/protein')!;
    expect(viewer.input.wheel(cellX(2), ProteinY, { deltaX: 0, deltaY: -10, ctrlKey: true })).toBe(true);
    frames.flush();
    const v1 = viewer.view.get('/views/main/protein')!;
    // exp(0.1): the window shrinks by that factor
    expect((v0.end - v0.start) / (v1.end - v1.start)).toBeCloseTo(Math.exp(0.1), 3);
    // out
    viewer.input.wheel(cellX(2), ProteinY, { deltaX: 0, deltaY: 10, ctrlKey: true });
    frames.flush();
    const v2 = viewer.view.get('/views/main/protein')!;
    expect(v2.end - v2.start).toBeCloseTo(v0.end - v0.start, 3);
  });

  it('a Ctrl+wheel with a mouse’s delta is the wheel with Ctrl, not a pinch; with the setting off, a pinch is too', async () => {
    const { viewer, frames } = await make();
    viewer.input.mac = true;
    const start = () => viewer.view.get('/views/main/protein')!;
    // elsewhere a Ctrl+wheel of a mouse's size is the wheel with Ctrl, which is the zoom modifier there
    viewer.input.mac = false;
    viewer.input.wheel(cellX(2), ProteinY, { deltaX: 0, deltaY: -100, ctrlKey: true });
    frames.flush();
    expect(4 / (start().end - start().start)).toBeCloseTo(1.25, 5);
    // on macOS Ctrl is no zoom modifier: a Ctrl+wheel is always a pinch, whatever its delta (exp(1) for 100 px)
    viewer.input.mac = true;
    viewer.view.fit('/views/main/protein');
    viewer.input.wheel(cellX(2), ProteinY, { deltaX: 0, deltaY: -100, ctrlKey: true });
    frames.flush();
    expect(4 / (start().end - start().start)).toBeCloseTo(Math.exp(1), 5);

    const off = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [...TinyConfig, [viewerItem('pointer.pinchFromCtrlWheel'), false]],
    });
    viewers.push(off);
    await off.load(proteinOnGene());
    off.resize(600);
    await off.settle();
    off.input.mac = true;
    off.input.wheel(cellX(2), ProteinY, { deltaX: 0, deltaY: -10, ctrlKey: true });
    frames.flush();
    expect(off.view.get('/views/main/protein')!.end).toBe(4);
  });

  it('keys run the commands with the shortcut on the section under the pointer', async () => {
    const { viewer, frames } = await make();
    viewer.input.mac = false;
    viewer.input.move(cellX(1), ProteinY);
    expect(viewer.input.key({ key: '=' })).toBe(true);
    frames.flush();
    const v = viewer.view.get('/views/main/protein')!;
    expect(v.end - v.start).toBeLessThan(4);
    expect(viewer.input.key({ key: '0' })).toBe(true);
    expect(viewer.input.key({ key: 'q' })).toBe(false);
    expect(viewer.view.get('/views/main/protein')!.end).toBe(4);
  });

  it('a double-click on a cell focuses it, on empty space fits', async () => {
    const { viewer } = await make();
    expect(viewer.input.doubleClick(cellX(1), ProteinY)).toBe(true);
    const v = viewer.view.get('/views/main/protein')!;
    expect([v.start, v.end]).toEqual([1, 2]);
    expect(viewer.input.doubleClick(40, ProteinY)).toBe(true);
    expect(viewer.view.get('/views/main/protein')!.end).toBe(4);
  });
});

describe('a resting pointer while the view moves (EC-143)', () => {
  it('is hit again once a frame: the hover follows the data under it', async () => {
    const { viewer, frames } = await make();
    viewer.input.move(cellX(0), ProteinY);
    frames.flush();
    expect(viewer.hover.current.position).toBe(0);
    // zoom in on the right half with the pointer resting: the cell under it changes with no pointer event
    viewer.view.setRange('/views/main/protein', 2, 4);
    frames.flush();
    expect(viewer.hover.current.position).toBe(2);
    viewer.view.pan('/views/main/protein', -1);
    frames.flush();
    expect(viewer.hover.current.position).toBe(1);
  });

  it('once a frame, not once a change: several moves of the view in one frame hit once', async () => {
    const { viewer, frames } = await make();
    viewer.input.move(cellX(0), ProteinY);
    frames.flush();
    const hits = vi.spyOn(viewer, 'hit');
    viewer.view.setRange('/views/main/protein', 1, 4);
    viewer.view.setRange('/views/main/protein', 2, 4);
    viewer.view.setRange('/views/main/protein', 3, 4);
    frames.flush();
    expect(hits).toHaveBeenCalledTimes(1);
    expect(viewer.hover.current.position).toBe(3);
  });

  it('the setting turns it off: the hover waits for the next pointer move', async () => {
    frames = installFrames();
    const v = createSeqViewer({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      config: [...TinyConfig, [viewerItem('pointer.rehit'), false]],
    });
    viewers.push(v);
    await v.load(proteinOnGene());
    v.resize(600);
    await v.settle();
    v.input.move(cellX(0), ProteinY);
    frames.flush();
    v.view.setRange('/views/main/protein', 2, 4);
    frames.flush();
    expect(v.hover.current.position).toBe(0);
    v.input.move(cellX(0) + 1, ProteinY);
    expect(v.hover.current.position).toBe(2);
  });

  it('a pointer that left is not hit again', async () => {
    const { viewer, frames } = await make();
    viewer.input.move(cellX(0), ProteinY);
    viewer.input.leave();
    frames.flush();
    viewer.view.setRange('/views/main/protein', 2, 4);
    frames.flush();
    expect(viewer.hover.active).toBe(false);
  });
});

describe('the picture follows a pan under a resting pointer (EC-143)', () => {
  it('draws the base canvas again: the layout moved under the pointer, whether or not the hover changed with it', async () => {
    const { viewer, backend, base, labels, frames } = await make({ canvas: true });
    viewer.input.move(cellX(0), ProteinY);
    frames.flush();
    const before = { ...backend.counts };
    base.context.take();
    labels.context.take();
    const frameCount = viewer.counts.frames;
    viewer.view.setRange('/views/main/protein', 2, 4);
    frames.flush();
    // the hover moved with the data (cell 0 → 2), and still the base was drawn for the new window
    expect(viewer.hover.current.position).toBe(2);
    expect(backend.counts.base).toBe(before.base + 1);
    expect(backend.counts.labels).toBe(before.labels + 1);
    expect(base.context.calls('fillText').length).toBeGreaterThan(0);
    // one frame did it: the re-hit asks for no frame of its own
    expect(viewer.counts.frames).toBe(frameCount + 1);
  });

  it('draws the base canvas when the pan leaves the hover as it was', async () => {
    const { viewer, backend, frames } = await make({ canvas: true });
    // over the label column: nothing to hit before or after
    viewer.input.move(40, ProteinY);
    frames.flush();
    const before = backend.counts.base;
    viewer.view.setRange('/views/main/protein', 2, 4);
    frames.flush();
    expect(backend.counts.base).toBe(before + 1);
  });
});

describe('behaviors the spec does not list', () => {
  it('a gesture that needs one is an issue, once, and does nothing', async () => {
    frames = installFrames();
    const spec = DefaultSeqViewerSpec();
    const v = createSeqViewer({ ...spec, ui: undefined, behaviors: [], config: TinyConfig });
    viewers.push(v);
    await v.load(proteinOnGene());
    v.resize(600);
    await v.settle();
    v.input.move(cellX(1), ProteinY);
    v.input.move(cellX(2), ProteinY);
    expect(v.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: 10 })).toBe(false);
    expect(v.hover.active).toBe(false);
    expect(v.behaviorIssues.value.map((i) => [i.code, i.address])).toEqual([
      ['seq.viewer.behavior.missing', 'hover'],
      ['seq.viewer.behavior.missing', 'navigate'],
    ]);
  });

  it('finds a behavior by its provider, whatever the spec calls it', async () => {
    const { viewer } = await make();
    expect(viewer.behavior(HoverProvider)).toBeDefined();
    expect(viewer.behaviorIssues.value).toEqual([]);
  });
});

describe('wheel details', () => {
  it('a line is pointer.lineHeight px, a page the viewport’s own axis', async () => {
    const { viewer, frames } = await make();
    viewer.input.mac = false;
    const span = () => {
      const v = viewer.view.get('/views/main/protein')!;
      return v.end - v.start;
    };
    viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: -3, ctrlKey: true, deltaMode: 1 });
    frames.flush();
    // 3 lines of 16 px: the wheel step's power of 48 / 100
    expect(4 / span()).toBeCloseTo(1.25 ** 0.48, 5);
    // a page: scaled by the height for a vertical delta, by the width for a horizontal one
    viewer.view.fit('/views/main/protein');
    const before = viewer.view.get('/views/main/protein')!.start;
    viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: -0.5, ctrlKey: true, deltaMode: 2 });
    frames.flush();
    expect(4 / span()).toBeCloseTo(1.25 ** ((viewer.height * 0.5) / 100), 5);
    expect(before).toBe(0);
  });

  it('over the label column the zoom is about the section’s first visible position', async () => {
    const { viewer, frames } = await make();
    viewer.input.mac = false;
    viewer.view.setRange('/views/main/protein', 1, 3);
    frames.flush();
    viewer.input.wheel(10, ProteinY, { deltaX: 0, deltaY: -100, ctrlKey: true });
    frames.flush();
    const v = viewer.view.get('/views/main/protein')!;
    expect(v.start).toBeCloseTo(1, 6);
    expect(v.end).toBeLessThan(3);
  });
});

describe('rested', () => {
  // data arriving while the view is still is driven by the Gene by region story's test (a delayed source)
  it('moves once each time the picture rests, after the settled frame is drawn: a move of the view, a resize', async () => {
    const { viewer, frames } = await make();
    /** The settle timers run, then the frames they ask for. */
    const rest = () => {
      frames.timers.fire();
      frames.flush();
      frames.timers.fire();
      frames.flush();
    };
    rest();
    let rested = 0;
    const off = viewer.rested.listen(() => rested++);

    viewer.view.setRange('/views/main/protein', 1, 3);
    frames.flush();
    expect(rested).toBe(0);
    rest();
    expect(rested).toBe(1);

    viewer.resize(720);
    frames.flush();
    rest();
    expect(rested).toBe(2);

    // nothing changed: nothing rests again
    rest();
    expect(rested).toBe(2);

    // a move before the settled frame is drawn: no rest for the first, one for the second
    viewer.view.setRange('/views/main/protein', 0, 2);
    frames.flush();
    frames.timers.fire();
    viewer.view.setRange('/views/main/protein', 1, 3);
    frames.flush();
    expect(rested).toBe(2);
    rest();
    expect(rested).toBe(3);
    off();
  });
});

describe('load', () => {
  it('opens the sources once', async () => {
    frames = installFrames();
    const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(v);
    const sync = vi.spyOn(v.store, 'sync');
    await v.load(proteinOnGene());
    expect(sync).toHaveBeenCalledTimes(1);
  });
});

describe('retire', () => {
  it('leaves no frame, timer or listener behind: a frame asked for before does nothing', async () => {
    const { viewer, base, frames } = await make({ canvas: true });
    const timers = frames.timers;
    viewer.input.move(cellX(1), ProteinY);
    viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: -100, ctrlKey: true, deltaMode: 1 });
    // a frame and the settle timer are pending (and the rest timer of the load's last layout)
    expect(frames.callbacks).toBeGreaterThan(0);
    expect(timers.live).toBeGreaterThanOrEqual(1);
    const drawn = viewer.counts.frames;
    base.context.take();
    viewer.retire();
    // retire cancels the frame it asked for, and the wheel's own: none is pending, no timer is live
    expect(frames.callbacks).toBe(0);
    expect(timers.live).toBe(0);
    frames.flush();
    expect(viewer.counts.frames).toBe(drawn);
    expect(base.context.take()).toEqual([]);
    // events after retire change nothing
    viewer.input.wheel(cellX(1), ProteinY, { deltaX: 0, deltaY: -100, ctrlKey: true, deltaMode: 1 });
    frames.flush();
    expect(viewer.counts.frames).toBe(drawn);
    expect(frames.callbacks).toBe(0);
    expect(frames.pending).toBe(0);
    expect(timers.live).toBe(0);
  });

  it('retired first: the UI’s teardown and the behaviors’ asks no frame', async () => {
    const { viewer, frames } = await make({ canvas: true });
    frames.flush();
    expect(frames.callbacks).toBe(0);
    let askedWhileRetiring = 0;
    const ask = vi.spyOn(viewer, 'schedule');
    viewer.useBackend(undefined);
    frames.flush();
    ask.mockClear();
    // what the disposers do: the backend is dropped, the models retire; none asks for a frame
    viewer.retire();
    askedWhileRetiring = frames.callbacks;
    expect(askedWhileRetiring).toBe(0);
    expect(viewer.retired).toBe(true);
  });

  it('disposes the ports it made and leaves a host’s', async () => {
    const owned = vi.fn();
    const kept = vi.fn();
    const spec = (disposeWithViewer: boolean, dispose: () => void) => ({
      ...DefaultSeqViewerSpec(),
      ui: undefined,
      ports: {
        assets: () => ({ manager: { dispose } as never, disposeWithViewer }),
      },
    });
    createSeqViewer(spec(true, owned)).retire();
    createSeqViewer(spec(false, kept)).retire();
    expect(owned).toHaveBeenCalledTimes(1);
    expect(kept).not.toHaveBeenCalled();
  });
});

describe('DomHost', () => {
  it('is the host’s clock and frames: the viewer asks no other', async () => {
    const raf = vi.spyOn(DomHost, 'requestAnimationFrame');
    const { viewer } = await make();
    viewer.schedule();
    expect(raf).toHaveBeenCalled();
  });
});
