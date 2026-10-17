// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EmptyLoci } from '@molstar/model/model/loci';
import { StructureElement } from '@molstar/model/model/structure';
import { afterAll, beforeAll, describe, expect, it, vi } from 'vitest';
import { ensemble } from '../../api/_test/trajectory-fixtures.js';
import { cbs, expectOk } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createDefaultPlugin } from '../../sources/create-plugin.js';
import { createViewerCore, type ViewerCore } from '../../viewer.js';
import { type MountedShell, mountShell } from '../shell.js';

// The idle and interaction budget (`test:perf`; prism AGENTS.md: no Prism code keeps a loop running while idle, and
// frame work joins the one frame tick, rendering.md §2.3). In the real shell with structures loaded: nothing changing,
// Prism asks for no animation frame and queues no microtask over 2 s; while input comes at two events a frame (a hover
// sweep over the structure, the scene tree scrolled, a splitter dragged), Prism asks for at most one frame a browser
// frame. Prism's are the calls whose stack runs through packages/prism or extensions/ (Mol*'s own are not counted; with
// no canvas mounted, it runs no loop here). The positive controls (a frame loop left running; a frame asked for per
// event) are caught.

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

const nativeRaf = globalThis.requestAnimationFrame.bind(globalThis);
const nativeMicrotask = globalThis.queueMicrotask.bind(globalThis);
const counts = { on: false, raf: 0, micro: 0, sites: new Map<string, number>() };
const isPrism = (stack: string) =>
  stack
    .split('\n')
    .slice(2)
    .find((l) => /\/packages\/prism\/|\/extensions\//.test(l));
globalThis.requestAnimationFrame = ((fn: FrameRequestCallback) => {
  if (counts.on) {
    const stack = new Error().stack ?? '';
    const at = isPrism(stack);
    if (at) {
      counts.raf++;
      // Named by the first frame outside the frame lane's own (who asked for it).
      const who =
        stack
          .split('\n')
          .slice(2)
          .find((l) => /\/packages\/prism\/|\/extensions\//.test(l) && !/reactive\/frame-state/.test(l)) ?? at;
      counts.sites.set(who.trim(), (counts.sites.get(who.trim()) ?? 0) + 1);
    }
  }
  return nativeRaf(fn);
}) as typeof requestAnimationFrame;
globalThis.queueMicrotask = ((fn: VoidFunction) => {
  if (counts.on && isPrism(new Error().stack ?? '')) counts.micro++;
  return nativeMicrotask(fn);
}) as typeof queueMicrotask;

const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
const nextFrame = () => new Promise<number>((r) => nativeRaf(r));

/** Counts over `run`: Prism's frame requests and microtasks, and the most frame requests between two browser frames. */
async function measure(run: () => Promise<void>) {
  counts.raf = 0;
  counts.micro = 0;
  counts.sites.clear();
  let most = 0;
  let last = 0;
  let frames = 0;
  let on = true;
  const tick = () => {
    if (!on) return;
    frames++;
    most = Math.max(most, counts.raf - last);
    last = counts.raf;
    nativeRaf(tick);
  };
  nativeRaf(tick);
  counts.on = true;
  try {
    await run();
    await nextFrame();
    await nextFrame();
  } finally {
    counts.on = false;
    on = false;
  }
  most = Math.max(most, counts.raf - last);
  return { raf: counts.raf, micro: counts.micro, most, frames, sites: [...counts.sites] };
}

let core: ViewerCore;
let mounted: MountedShell;
let host: HTMLElement;

beforeAll(async () => {
  vi.stubGlobal('ResizeObserver', NoResize);
  // jsdom has no pointer capture (a splitter's drag takes it).
  HTMLElement.prototype.setPointerCapture ??= () => {};
  HTMLElement.prototype.releasePointerCapture ??= () => {};
  HTMLElement.prototype.hasPointerCapture ??= () => false;
  core = await createViewerCore({
    spec: DefaultPrismViewerSpec(),
    createPlugin: createDefaultPlugin,
    modules: ['files', 'page'],
  });
  vi.spyOn(core, 'mountCanvas').mockResolvedValue(false);
  host = document.createElement('div');
  document.body.append(host);
  mounted = await mountShell(core, host);
  expectOk(await core.api.load.files({ files: [cbs()] }));
  expectOk(await core.api.load.files({ files: [ensemble(5)] }));
  mounted.shell.show('structure');
  // Settled: the loads' late work done (the sequence's first draw comes a while after, its chunk and its document in),
  // until 1.5 s pass with no Prism frame asked for (at most 15 s).
  for (let k = 0, quiet = 0; k < 30 && quiet < 3; k++) quiet = (await measure(() => wait(500))).raf > 0 ? 0 : quiet + 1;
}, 60_000);

afterAll(async () => {
  await wait(300);
  mounted.dispose();
  await core.dispose();
  host.remove();
  vi.unstubAllGlobals();
});

/** Two events a frame for `frames` frames. */
async function twiceAFrame(frames: number, event: (i: number) => void) {
  for (let i = 0; i < frames; i++) {
    event(2 * i);
    await wait(4);
    event(2 * i + 1);
    await nextFrame();
  }
}

const hoverSweep = () => {
  const structure = core.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  if (!structure) throw new Error('no structure');
  const all = StructureElement.Loci.all(structure);
  const first = StructureElement.Loci.firstResidue(all);
  return (i: number) => {
    // The label Mol* sets for what is under the pointer (its pick, on a canvas this test has not): the hover label's.
    core.plugin.behaviors.labels.highlight.next({ labels: i % 3 === 2 ? [] : [`Residue ${i}`] } as never);
    core.plugin.behaviors.interaction.hover.next({
      current: { loci: i % 3 === 2 ? EmptyLoci : i % 3 === 0 ? first : all },
      buttons: 0,
      button: 0,
      modifiers: { shift: false, control: false, alt: false, meta: false },
    } as never);
  };
};

describe('frames and microtasks', () => {
  it('idle: Prism asks for no frame and queues no microtask over 2 s', async () => {
    const r = await measure(() => wait(2000));
    expect(r.raf, `frames asked for: ${JSON.stringify(r.sites)}`).toBe(0);
    expect(r.micro).toBe(0);
  });

  it('a frame loop left running is caught (positive control)', async () => {
    let on = true;
    const loop = () => {
      if (on) requestAnimationFrame(loop);
    };
    const r = await measure(async () => {
      requestAnimationFrame(loop);
      await wait(500);
      on = false;
    });
    expect(r.raf).toBeGreaterThan(0);
  });

  it('a hover sweep at two events a frame: at most one Prism frame a browser frame', async () => {
    const hover = hoverSweep();
    const r = await measure(() => twiceAFrame(40, hover));
    expect(r.raf).toBeGreaterThan(0);
    expect(r.most, `the most in one frame: ${JSON.stringify(r.sites)}`).toBeLessThanOrEqual(1);
  });

  it('the scene tree scrolled and a splitter dragged, two events a frame: at most one Prism frame a browser frame', async () => {
    const viewport = mounted.root
      .querySelector<HTMLElement>('[role=tree]')
      ?.closest('[data-axis]')
      ?.querySelector<HTMLElement>(':scope > div');
    const splitter = mounted.root.querySelector<HTMLElement>('[role=separator]');
    let top = 0;
    if (viewport)
      Object.defineProperty(viewport, 'scrollTop', {
        get: () => top,
        set: (v: number) => (top = v),
        configurable: true,
      });
    splitter?.dispatchEvent(
      new PointerEvent('pointerdown', { bubbles: true, clientX: 300, clientY: 300, pointerId: 1 }),
    );
    const r = await measure(() =>
      twiceAFrame(40, (i) => {
        if (viewport) {
          top = (i % 20) * 26;
          viewport.dispatchEvent(new Event('scroll'));
        }
        splitter?.dispatchEvent(
          new PointerEvent('pointermove', { bubbles: true, clientX: 300 + i, clientY: 300, pointerId: 1 }),
        );
      }),
    );
    splitter?.dispatchEvent(new PointerEvent('pointerup', { bubbles: true, clientX: 380, clientY: 300, pointerId: 1 }));
    expect(viewport || splitter, 'something to scroll or drag').toBeTruthy();
    expect(r.most, `the most in one frame: ${JSON.stringify(r.sites)}`).toBeLessThanOrEqual(1);
  });

  it('a frame asked for per event is caught (positive control)', async () => {
    const hover = hoverSweep();
    const r = await measure(() =>
      twiceAFrame(10, (i) => {
        hover(i);
        requestAnimationFrame(() => {});
      }),
    );
    expect(r.most).toBeGreaterThan(1);
  });
});
