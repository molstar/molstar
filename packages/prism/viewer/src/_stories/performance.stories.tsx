/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { For } from '@molstar/prism-ui/dom/control-flow';
import { render } from '@molstar/prism-ui/dom/render';
import { flush, flushFrame } from '@molstar/prism-ui/dom/scope';
import { ReactiveCounts } from '@molstar/prism-ui/reactive/state';
import type { StateTreeModel, StateTreeNode } from '../model/state-tree.js';
import { type ObservableState, observableState } from '@molstar/prism-ui/reactive/observable-state';
import { flattenTree, RecycledListStats, VirtualList } from '@molstar/prism-ui';
import { StateTreeTab } from '../views/state-tree/state-tree-tab.js';
import { SceneTree } from '../views/structure/scene-tree.js';
import { makeScene } from './performance-fixtures.js';
import { story } from '@molstar/prism-ui/testing/index';

import { derive } from '@molstar/prism-ui/reactive/state';
const meta: Meta = { title: 'Internals/Performance/Trees' };
export default meta;

// The renderer check: the state tree at 5,000 cells under bursts of updates, in the `dom` runtime.
// Every number is wall time in this browser (ms), layout forced, measured over a frame.

// Two frames when the page shows; a hidden page gets no frames and throttled timers, so there a microtask (script,
// style and layout counted, layout forced; paint not).
const frame = () =>
  document.hidden
    ? Promise.resolve()
    : new Promise<void>((r) => requestAnimationFrame(() => requestAnimationFrame(() => r())));

/** 25 structures × 4 models × 49 cells under each (+ the models and structures): 5,025 cells. */
function makeTree(structures = 25): StateTreeNode[] {
  const roots: StateTreeNode[] = [];
  for (let s = 0; s < structures; s++) {
    const models: StateTreeNode[] = [];
    for (let m = 0; m < 4; m++) {
      const cells: StateTreeNode[] = [];
      for (let c = 0; c < 49; c++)
        cells.push({
          ref: `s${s}m${m}c${c}`,
          label: `Component ${c}`,
          description: `${(c + 1) * 13} elements`,
          icon: 'cartoon',
          transform: 'Component',
          visible: c % 7 !== 0,
        });
      models.push({ ref: `s${s}m${m}`, label: `Model ${m + 1}`, icon: 'layers', transform: 'Model', children: cells });
    }
    roots.push({
      ref: `s${s}`,
      label: `Structure ${s + 1}`,
      icon: 'structure',
      transform: 'Structure',
      children: models,
    });
  }
  return roots;
}

function count(nodes: readonly StateTreeNode[]): number {
  let n = 0;
  for (const x of nodes) n += 1 + count(x.children ?? []);
  return n;
}

/** Path-copying update of one cell's label (how an immutable model hands out a change). */
function relabel(nodes: readonly StateTreeNode[], ref: string, label: string): StateTreeNode[] {
  return nodes.map((n) => {
    if (n.ref === ref) return { ...n, label };
    if (n.children && ref.startsWith(n.ref)) return { ...n, children: relabel(n.children, ref, label) };
    return n;
  });
}

function treeModel(roots: ObservableState<readonly StateTreeNode[]>): StateTreeModel {
  const selected = observableState<string | undefined>(undefined);
  return {
    roots: roots.asReadonly,
    selected: selected.asReadonly,
    select: (ref) => selected.set(ref),
    inspect: () => undefined,
    toggleVisible: () => {},
    remove: () => {},
  };
}

interface Result {
  case: string;
  /** Until the result is all there (a sliced list renders over several frames). */
  ms: number;
  /** The longest task in that time (how long the page could not respond); 0 when none passed 50 ms. */
  longestTask: number;
  note?: string;
}

/** Results as they come (window.__mspPerf), so a slow case does not hide the others. */
const log: Result[] = [];

async function time(
  label: string,
  run: () => void | Promise<void>,
  host: HTMLElement,
  note?: string,
  settled: () => boolean = () => true,
): Promise<Result> {
  await frame();
  let longest = 0;
  const observer =
    typeof PerformanceObserver !== 'undefined' && PerformanceObserver.supportedEntryTypes?.includes('longtask')
      ? new PerformanceObserver((list) => {
          for (const e of list.getEntries()) longest = Math.max(longest, e.duration);
        })
      : undefined;
  observer?.observe({ type: 'longtask' });
  const t0 = performance.now();
  await run();
  const tRun = performance.now();
  void host.offsetHeight;
  const tLayout = performance.now();
  await frame();
  (window as any).__mspSplit = {
    label,
    run: Math.round(tRun - t0),
    layout: Math.round(tLayout - tRun),
    frame: Math.round(performance.now() - tLayout),
  };
  // A sliced list: wait until all of it is there.
  for (let n = 0; n < 600 && !settled(); n++) await frame();
  const ms = Math.round(performance.now() - t0);
  await frame();
  observer?.disconnect();
  const r = { case: label, ms, longestTask: Math.round(longest), note };
  log.push(r);
  (window as any).__mspPerf = { hidden: document.hidden, results: log.slice() };
  return r;
}

/** How many relabels a burst makes (the immutable case is slow enough that 200 can block the page for a minute). */
const Burst = 20;

/** The real State tree tab (its rows, CSS, subscriptions) at 5,025 cells. */
async function realTab(host: HTMLElement, out: Result[]) {
  const roots = observableState<readonly StateTreeNode[]>([]);
  const model = treeModel(roots);
  host.replaceChildren();
  const box = document.createElement('div');
  box.style.cssText = 'height: 600px; width: 320px; display: flex; flex-direction: column;';
  host.append(box);
  const dispose = render(() => <StateTreeTab model={model} />, box);
  const tree = makeTree();
  const cells = count(tree);
  // Virtualized (view-models.md §3): the window is all there is to render.
  const all = () => box.querySelectorAll('[role=treeitem]').length > 10;
  out.push(
    await time(
      `State tree tab: first render, ${cells} cells`,
      () => roots.set(tree),
      host,
      'only the window is rendered',
      all,
    ),
  );
  out.push(await time('State tree tab: select a cell (every row re-checks)', () => model.select('s12m2c30'), host));
  // A burst as path-copying records: 200 label changes, one per tick of the source.
  out.push(
    await time(
      `State tree tab: ${Burst} relabels as immutable records (each copies its path)`,
      async () => {
        let t = roots.value;
        for (let i = 0; i < Burst; i++) {
          t = relabel(t, `s${i % 25}m${i % 4}c${i % 49}`, `Renamed ${i}`);
          roots.set(t);
        }
      },
      host,
      'a changed cell and its ancestors render again, if on screen; nothing else',
      all,
    ),
  );
  // The structures' rows (level 1) as they come into the window: closed one by one, then opened again.
  const toggleAll = (open: boolean) => {
    const at = `[aria-level="1"][aria-expanded="${open}"] [data-tree-toggle]`;
    for (let b = box.querySelector<HTMLButtonElement>(at), n = 0; b && n < 100; b = box.querySelector(at), n++)
      b.click();
  };
  out.push(
    await time(
      'State tree tab: collapse and expand the 25 structures',
      () => {
        toggleAll(true);
        toggleAll(false);
      },
      host,
    ),
  );
  dispose();
}

/** The real Scene tree (its rows, CSS, subscriptions) at 1,000 components, each with 3 representations. */
async function sceneTree(host: HTMLElement, out: Result[]) {
  const fixture = makeScene(1000, 3);
  host.replaceChildren();
  const box = document.createElement('div');
  box.style.cssText = 'height: 600px; width: 320px; display: flex; flex-direction: column;';
  host.append(box);
  const dispose = render(() => <SceneTree scene={fixture.scene} />, box);
  const shown = () => box.querySelectorAll('[role=treeitem]').length > 10;
  // The structure is open from the start: 1,001 rows (+ Contents), the representations closed under their components.
  out.push(
    await time(
      `Scene tree: first render, ${fixture.components.length} components`,
      () => {},
      host,
      'only the window is rendered',
      shown,
    ),
  );
  out.push(await time('Scene tree: pick a component', () => fixture.scene.pick(fixture.components[500].node), host));
  out.push(
    await time(
      `Scene tree: ${Burst} relabels`,
      () => {
        for (let i = 0; i < Burst; i++) fixture.components[i * 7].label.set(`Renamed ${i}`);
      },
      host,
    ),
  );
  out.push(
    await time(
      'Scene tree: open and close 100 components',
      () => {
        for (const c of fixture.components.slice(0, 100)) c.expanded.set(true);
        for (const c of fixture.components.slice(0, 100)) c.expanded.set(false);
      },
      host,
      'each a flattening of the shown rows',
    ),
  );
  out.push(
    await time(
      `Scene tree: open all ${fixture.components.length} components (${fixture.representations} representation rows)`,
      () => {
        for (const c of fixture.components) c.expanded.set(true);
      },
      host,
      undefined,
      shown,
    ),
  );
  const viewport = box.querySelector<HTMLElement>('[data-axis] > div')!;
  out.push(
    await time(
      'Scene tree: scroll to the bottom',
      () => {
        viewport.scrollTop = viewport.scrollHeight;
      },
      host,
    ),
  );
  out.push(
    await time(
      'Scene tree: close all components',
      () => {
        for (const c of fixture.components) c.expanded.set(false);
      },
      host,
    ),
  );
  dispose();
}

/** The same rows as stable entries that follow their record by id: the shape sources would hand out. */
async function stableRows(host: HTMLElement, out: Result[]) {
  interface Row {
    ref: string;
    depth: number;
  }
  const rows: Row[] = [];
  const records = new Map<string, ObservableState<StateTreeNode>>();
  const walk = (nodes: readonly StateTreeNode[], depth: number) => {
    for (const n of nodes) {
      rows.push({ ref: n.ref, depth });
      records.set(n.ref, observableState(n));
      walk(n.children ?? [], depth + 1);
    }
  };
  walk(makeTree(), 0);
  const list = observableState<readonly Row[]>([]);
  host.replaceChildren();
  const box = document.createElement('div');
  box.style.cssText = 'height: 600px; width: 320px; overflow: auto;';
  host.append(box);
  const Line = (props: { row: Row; record: Observable<StateTreeNode> }) => (
    <div
      style={`display: flex; gap: var(--msp-space-3); height: 26px; align-items: center; padding-left: ${props.row.depth * 12}px`}
    >
      <span>
        {props.record.pipe(
          map((n) => n.label),
          distinctUntilChanged(),
        )}
      </span>
      <span style="color: var(--msp-color-text-subtle)">
        {props.record.pipe(
          map((n) => n.description ?? ''),
          distinctUntilChanged(),
        )}
      </span>
    </div>
  );
  const dispose = render(
    () => (
      <For each={list} key={(r) => r.ref}>
        {(row) => <Line row={row} record={records.get(row.ref)!.asReadonly} />}
      </For>
    ),
    box,
  );
  out.push(
    await time(
      `Stable rows: first render, ${rows.length} rows`,
      () => list.set(rows),
      host,
      undefined,
      () => box.childElementCount >= rows.length,
    ),
  );
  out.push(
    await time(
      `Stable rows: ${Burst} relabels, each a record of its own`,
      () => {
        for (let i = 0; i < Burst; i++) {
          const ref = `s${i % 25}m${i % 4}c${i % 49}`;
          const r = records.get(ref)!;
          r.set({ ...r.value, label: `Renamed ${i}` });
        }
      },
      host,
      'only the changed rows update; no entry is re-created',
    ),
  );
  dispose();
}

/**
 * The tree virtualized: its shown nodes flattened into rows, only the window rendered (`VirtualList`), each row following
 * its node's own record (per-node records, as sources hand them out), expansion a set of refs.
 */
async function virtualTree(host: HTMLElement, out: Result[]) {
  const records = new Map<string, ObservableState<StateTreeNode>>();
  const index = (nodes: readonly StateTreeNode[]) => {
    for (const n of nodes) {
      records.set(n.ref, observableState(n));
      index(n.children ?? []);
    }
  };
  const roots = makeTree();
  index(roots);
  const collapsed = observableState<ReadonlySet<string>>(new Set());
  const selected = observableState<string | undefined>(undefined);
  const rows = observableState<readonly { node: StateTreeNode; depth: number }[]>([]);
  const reflow = () =>
    rows.set(
      flattenTree(
        roots,
        (n) => n.children,
        (n) => !collapsed.value.has(n.ref),
      ),
    );
  host.replaceChildren();
  const box = document.createElement('div');
  box.style.cssText = 'height: 600px; width: 320px; display: flex; flex-direction: column;';
  host.append(box);
  const Row = (props: { ref: string; depth: number }) => {
    const record = records.get(props.ref)!.asReadonly;
    return (
      <div
        role="treeitem"
        aria-selected={derive(selected.asReadonly, (s) => s === props.ref)}
        style={`display: flex; gap: var(--msp-space-3); height: 26px; align-items: center; padding-left: ${props.depth * 14 + 4}px`}
      >
        <span>{derive(record, (n) => n.label)}</span>
        <span style="color: var(--msp-color-text-subtle)">{derive(record, (n) => n.description ?? '')}</span>
      </div>
    );
  };
  const dispose = render(
    () => (
      <VirtualList items={rows.asReadonly} key={(r) => r.node.ref} rowHeight={26} role="tree" aria-label="Tree">
        {(r) => <Row ref={r.node.ref} depth={r.depth} />}
      </VirtualList>
    ),
    box,
  );
  const shown = () => box.querySelectorAll('[role=treeitem]').length > 10;
  out.push(
    await time(`Virtual tree: first render, ${records.size} cells`, reflow, host, 'only the window is rendered', shown),
  );
  out.push(await time('Virtual tree: select a cell', () => selected.set('s0m0c3'), host));
  out.push(
    await time(
      `Virtual tree: ${Burst} relabels, per-node records`,
      () => {
        for (let i = 0; i < Burst; i++) {
          const r = records.get(`s0m${i % 4}c${i % 49}`)!;
          r.set({ ...r.value, label: `Renamed ${i}` });
        }
      },
      host,
    ),
  );
  out.push(
    await time(
      'Virtual tree: collapse and expand the 25 structures',
      () => {
        collapsed.set(new Set(roots.map((r) => r.ref)));
        reflow();
        collapsed.set(new Set());
        reflow();
      },
      host,
    ),
  );
  const viewport = box.querySelector<HTMLElement>('[data-axis] > div')!;
  out.push(
    await time(
      'Virtual tree: scroll to the bottom',
      () => {
        viewport.scrollTop = viewport.scrollHeight;
      },
      host,
    ),
  );
  dispose();
}

/** A button that runs `cases` in the fixed host and the table of their results. */
function benchmark(cases: (host: HTMLElement, out: Result[]) => Promise<void>) {
  return story(() => {
    const results = observableState<readonly Result[]>([]);
    const running = observableState(false);
    // A fixed, contained area like the dock's panel, outside Storybook's root (whose addons watch every change there),
    // themed like the story.
    const host = document.createElement('div');
    host.className = 'msp-root';
    host.style.cssText = 'position: fixed; right: 16px; top: 16px; width: 320px; height: 600px; contain: strict;';
    const run = async () => {
      document.body.append(host);
      running.set(true);
      const out: Result[] = [];
      log.length = 0;
      await cases(host, out);
      host.replaceChildren();
      host.remove();
      results.set(out);
      running.set(false);
      (window as any).__mspPerf = { hidden: document.hidden, results: out, done: true };
    };
    return (
      <div style="display: grid; gap: var(--msp-space-4); padding: var(--msp-space-5); font-size: var(--msp-font-size-sm); color: var(--msp-color-text);">
        <button
          type="button"
          style="justify-self: start; height: 28px; padding: 0 var(--msp-space-4); border-radius: 6px; background: var(--msp-color-control);"
          disabled={running.asReadonly}
          onClick={run}
        >
          {derive(running.asReadonly, (r) => (r ? 'Running…' : 'Run'))}
        </button>
        <table style="border-collapse: collapse;">
          <For each={results.asReadonly}>
            {(r) => (
              <tr>
                <td style="padding: var(--msp-space-1) var(--msp-space-4) var(--msp-space-1) 0">{r.case}</td>
                <td style="padding: var(--msp-space-1) var(--msp-space-4); text-align: right; font-variant-numeric: tabular-nums">
                  {r.ms} ms
                </td>
                <td style="padding: var(--msp-space-1) var(--msp-space-4); text-align: right; font-variant-numeric: tabular-nums">
                  {r.longestTask ? `longest task ${r.longestTask} ms` : 'no long task'}
                </td>
                <td style="padding: var(--msp-space-1) 0; color: var(--msp-color-text-subtle)">{r.note ?? ''}</td>
              </tr>
            )}
          </For>
        </table>
      </div>
    );
  });
}

// Storybook's accessibility checks run on every change in the story: at 5,000 rows they, not the UI, would be timed.
const noA11y = { a11y: { disable: true, test: 'off' } };

/**
 * The renderer check: the state tree at 5,025 cells in the `dom` runtime. Run it, then read the table (wall time in
 * this browser, layout forced).
 */
export const StateTree: StoryObj = {
  name: 'State tree',
  parameters: noA11y,
  render: () =>
    benchmark(async (host, out) => {
      await realTab(host, out);
      await stableRows(host, out);
      await virtualTree(host, out);
    }),
};

/** The same check for the scene tree: 1,000 components (3,000 representation rows when open), bursts of changes. */
export const SceneTreeBenchmark: StoryObj = {
  name: 'Scene tree',
  parameters: noA11y,
  render: () => benchmark(sceneTree),
};

/**
 * Both real trees at their largest, to use by hand: scroll, and the keyboard (Up/Down, Home/End, Left/Right to close
 * and open, Enter to select). The readout counts the rows in the DOM against the rows of the tree. Storybook's
 * accessibility checks stay on: a window of rows is all they see.
 */
export const Interactive: StoryObj = {
  name: 'Virtualized trees (interactive)',
  render: () =>
    story(() => {
      const fixture = makeScene(1000, 3);
      const roots = observableState<readonly StateTreeNode[]>(makeTree());
      const model = treeModel(roots);
      const cells = count(roots.value);
      const panel = (title: string, body: () => HTMLElement) => {
        const readout = document.createElement('div');
        readout.style.cssText = 'font-variant-numeric: tabular-nums; color: var(--msp-color-text-subtle);';
        const box = document.createElement('div');
        box.style.cssText =
          'height: 600px; width: 320px; display: flex; flex-direction: column; border: 1px solid var(--msp-color-border, #8884); --msp-surface: var(--msp-color-bg);';
        const wrap = document.createElement('div');
        wrap.style.cssText = 'display: grid; gap: var(--msp-space-2); align-content: start;';
        const heading = document.createElement('strong');
        heading.textContent = title;
        wrap.append(heading, readout, box);
        const el = body();
        box.append(el);
        // The readout follows the window's rows (a mutation of the list's own window, not a poll of the page: the poll
        // and a body-wide observer cost 5–10 ms a call, a hitch while scrolling).
        let watched: MutationObserver | undefined;
        let pending = false;
        const update = () => {
          pending = false;
          const tree = box.querySelector<HTMLElement>('[role=tree]');
          if (!tree) return;
          const rendered = tree.querySelectorAll('[role=treeitem]').length;
          const total = Math.round((parseFloat(tree.style.height) || 0) / 26);
          readout.textContent = `${rendered} rows rendered of ${total}`;
          if (!watched) {
            watched = new MutationObserver(() => {
              if (!wrap.isConnected) return watched?.disconnect();
              if (pending) return;
              pending = true;
              requestAnimationFrame(update);
            });
            watched.observe(tree, { childList: true, subtree: true, attributes: true, attributeFilter: ['style'] });
          }
        };
        // The tree mounts a moment later; until it does, look once a frame (for a few seconds at most).
        let tries = 300;
        const first = () => {
          if (box.querySelector('[role=tree]')) update();
          else if (--tries > 0) requestAnimationFrame(first);
        };
        requestAnimationFrame(first);
        return wrap;
      };
      const mount = (view: () => Node) => () => {
        const el = document.createElement('div');
        el.style.cssText = 'display: contents;';
        // Rendered into a node of its own, the tree fills the box.
        render(view, el);
        return el;
      };
      return (
        <div
          class="msp-root"
          style="display: flex; gap: var(--msp-space-5); padding: var(--msp-space-5); font-size: var(--msp-font-size-sm); color: var(--msp-color-text);"
        >
          {panel(
            `Scene tree: ${fixture.components.length} components, ${fixture.representations} representations`,
            mount(() => <SceneTree scene={fixture.scene} />),
          )}
          {panel(
            `State tree: ${cells} cells`,
            mount(() => <StateTreeTab model={model} />),
          )}
        </div>
      );
    }),
};

/**
 * Probes for profiling from the console: render the real tab at `structures` × 201 cells, optionally with the fades
 * (scroll-driven animations) switched off; returns ms (script, style, layout).
 */
(window as any).__mspProbe = async (structures: number, fades = true, frames = false) => {
  const host = document.createElement('div');
  document.body.append(host);
  const off = document.createElement('style');
  off.textContent = '* { animation: none !important; mask-image: none !important; }';
  if (!fades) document.head.append(off);
  const roots = observableState<readonly StateTreeNode[]>([]);
  const box = document.createElement('div');
  box.style.cssText = 'height: 600px; width: 320px; display: flex; flex-direction: column;';
  host.append(box);
  const dispose = render(() => <StateTreeTab model={treeModel(roots)} />, box);
  if (frames) await frame();
  const t0 = performance.now();
  roots.set(makeTree(structures));
  const tScript = performance.now() - t0;
  void box.offsetHeight;
  const tLayout = performance.now() - t0;
  if (frames) await frame();
  const total = performance.now() - t0;
  dispose();
  host.remove();
  off.remove();
  return {
    cells: structures * 201,
    script: Math.round(tScript),
    layout: Math.round(tLayout - tScript),
    frames: Math.round(total - tLayout),
  };
};

// --- Rows made -----------------------------------------------------------------------------------

/**
 * What a row costs to make, for the trees that make rows as they scroll: both real trees fully open (the scene tree's
 * 1,000 components and 3,000 representations, the state tree's 5,025 cells), mounted outside Storybook's root (its
 * highlight addon re-scans the page's styles on every row added or removed there, ~13 ms) with a11y off. Probes for
 * a driver (`__mspRows`; Playwright adds garbage collection and the heap, CDP): `jump` sets the scroll far away and
 * dispatches its event, so the rows made are timed in the same task (script, then layout); `steps` scrolls a row at a
 * time; `mount` leaves a tall tree mounted (a window of ~900 rows) for the heap.
 */
type TreeKind = 'scene' | 'state';

/** The scene tree's node rows: recycled slots (`RecycledList`, the default), or keyed rows with `?recycle=0`. */
const recycle = !(typeof location !== 'undefined' && /[?&]recycle=0/.test(location.search));

function mountTree(kind: TreeKind, height = 600) {
  const host = document.createElement('div');
  host.className = 'msp-root';
  host.style.cssText = `position: fixed; left: 16px; top: 16px; width: 320px; height: ${height}px; contain: strict;`;
  document.body.append(host);
  const box = document.createElement('div');
  box.style.cssText = `height: ${height}px; width: 320px; display: flex; flex-direction: column;`;
  host.append(box);
  let dispose: () => void;
  if (kind === 'scene') {
    const fixture = makeScene(1000, 3);
    for (const c of fixture.components) c.expanded.set(true);
    const fx = fixture;
    dispose = render(() => <SceneTree scene={fx.scene} recycle={recycle} />, box);
  } else {
    const roots = observableState<readonly StateTreeNode[]>(makeTree());
    dispose = render(() => <StateTreeTab model={treeModel(roots)} />, box);
  }
  const viewport = () => box.querySelector<HTMLElement>('[data-axis] > div')!;
  const rows = () => box.querySelectorAll('[role=treeitem]');
  return {
    host,
    viewport,
    rows,
    dispose: () => {
      dispose();
      host.remove();
    },
  };
}

const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
const round = (x: number, d = 3) => +x.toFixed(d);
const pct = (xs: number[], p: number) => {
  const s = xs.slice().sort((a, b) => a - b);
  return s[Math.min(s.length - 1, Math.floor(s.length * p))] ?? 0;
};

async function settle(t: ReturnType<typeof mountTree>) {
  for (let n = 0; n < 300 && t.rows().length < 10; n++) await frame();
  // The window grows once the scrolling rests (VirtualList, 120 ms).
  await wait(200);
  await frame();
}

/** Far jumps (a scrollbar's drag): each makes a window of new rows in the scroll event's task. */
async function jumps(kind: TreeKind, count = 40) {
  const t = mountTree(kind);
  await settle(t);
  const vp = t.viewport();
  const span = vp.scrollHeight - vp.clientHeight;
  const script: number[] = [];
  const layout: number[] = [];
  const total: number[] = [];
  let made = 0;
  for (let k = 0; k < count; k++) {
    const before = new Set(t.rows());
    // Spread over the list, never near the last one.
    const y = Math.round((((k * 37) % 97) / 97) * span);
    const t0 = performance.now();
    vp.scrollTop = y;
    vp.dispatchEvent(new Event('scroll'));
    flushFrame();
    flush();
    const t1 = performance.now();
    // A row's place: the window's layout (the host is contained, its own size needs none).
    const all = t.rows();
    void (all[all.length - 1] as HTMLElement | undefined)?.offsetTop;
    const t2 = performance.now();
    let n = 0;
    for (const r of t.rows()) if (!before.has(r)) n++;
    made += n;
    if (k > 1) {
      script.push(t1 - t0);
      layout.push(t2 - t1);
      total.push(t2 - t0);
    }
    await frame();
  }
  t.dispose();
  const sum = (xs: number[]) => xs.reduce((a, b) => a + b, 0);
  const rowsPerJump = made / count;
  return {
    kind,
    jumps: count,
    rowsPerJump: round(rowsPerJump, 1),
    scriptPerJump: round(sum(script) / script.length, 2),
    layoutPerJump: round(sum(layout) / layout.length, 2),
    scriptPerRow: round(sum(script) / script.length / rowsPerJump),
    p95Script: round(pct(script, 0.95), 2),
    medianJumpMs: round(pct(total, 0.5), 2),
    p90JumpMs: round(pct(total, 0.9), 2),
  };
}

/** Small steps (a wheel's notch, the keyboard): a row at a time, `count` of them. */
async function steps(kind: TreeKind, count = 300, by = 26) {
  const t = mountTree(kind);
  await settle(t);
  const vp = t.viewport();
  vp.scrollTop = Math.round((vp.scrollHeight - vp.clientHeight) / 3);
  vp.dispatchEvent(new Event('scroll'));
  await wait(200);
  await frame();
  const times: number[] = [];
  for (let k = 0; k < count; k++) {
    const t0 = performance.now();
    vp.scrollTop += by;
    vp.dispatchEvent(new Event('scroll'));
    flushFrame();
    flush();
    void (vp.querySelector('[role=treeitem]') as HTMLElement | null)?.offsetTop;
    times.push(performance.now() - t0);
    if (k % 10 === 9) await frame();
  }
  t.dispose();
  return {
    kind,
    steps: count,
    meanMs: round(times.reduce((a, b) => a + b, 0) / count, 3),
    p95Ms: round(pct(times, 0.95), 2),
    maxMs: round(Math.max(...times), 2),
  };
}

let mounted: ReturnType<typeof mountTree> | undefined;

// --- The update path: a mounted row, from a state's change to the DOM showing it -----------------

/**
 * The scene tree with its 1,000 components closed (each with a swatch), `height` px tall (600: about 48 rows in the
 * window), for `updateCase`: one row's own states (its label, hidden, swatch), the picked node moving between two
 * rows (every row reads it), the tree's places map replaced by a reflow far below (the shown rows' places stay), and
 * bursts of label sets in one task and across tasks. Each case reports per update: script until `set` returns
 * (`syncMs`) and until the DOM is right, after `flush()` (`domMs`); whether a read right after `set` saw the old DOM
 * (`stale`); and the runtime's counts (`ReactiveCounts`: listener calls, derived recomputes, `watch` renders).
 */
let upd: { fx: ReturnType<typeof makeScene>; box: HTMLElement; dispose(): void } | undefined;

const lineOf = (id: string) => upd!.box.querySelector<HTMLElement>(`[data-node-id="${id}"]`)!;
const labelOf = (id: string) => lineOf(id).querySelector('[data-tip-part]')!.textContent;

type UpdateCase = 'label' | 'hidden' | 'swatch' | 'pick' | 'reflow' | 'burst' | 'burstTasks';

/**
 * One update of a case (`k` its count): what it sets, and whether the DOM shows it. Bursts are 100 label sets over the
 * first 10 rows (in the window at every size).
 */
function updateOf(name: Exclude<UpdateCase, 'burstTasks'>, k: number): { mutate: () => void; shows: () => boolean } {
  const { fx } = upd!;
  const c = fx.components;
  const tick = `${name} ${k}`;
  switch (name) {
    case 'label':
      return { mutate: () => c[3].label.set(tick), shows: () => labelOf('c3') === tick };
    case 'hidden': {
      const to = k % 2 === 0;
      return {
        mutate: () => c[3].hidden.set(to),
        shows: () => /hidden/.test(lineOf('c3').className) === to,
      };
    }
    case 'swatch': {
      const to = k % 2 ? 0xff0000 : 0x00ff00;
      const rgb = to === 0xff0000 ? '255, 0, 0' : '0, 255, 0';
      return {
        mutate: () => c[3].swatch!.set(to),
        shows: () =>
          !!lineOf('c3').querySelector<HTMLElement>('[aria-hidden="true"][style]')?.style.background.includes(rgb),
      };
    }
    case 'pick': {
      const [to, from] = k % 2 ? ['c3', 'c5'] : ['c5', 'c3'];
      return {
        mutate: () => fx.scene.pick(c[to === 'c3' ? 3 : 5].node),
        shows: () => /picked/.test(lineOf(to).className) && !/picked/.test(lineOf(from).className),
      };
    }
    case 'reflow':
      return { mutate: () => c[900].expanded.set(k % 2 === 0), shows: () => true };
    case 'burst':
      return {
        mutate: () => {
          for (let i = 0; i < 100; i++) c[i % 10].label.set(`${tick} ${i}`);
        },
        shows: () => labelOf('c9') === `${tick} 99`,
      };
  }
}

/**
 * `n` updates of a case. Timed as loops (the page's clock is coarse: 0.1 ms): `domMs`, each update set and flushed;
 * `syncMs`, the sets alone (flushed once after the loop, untimed). `stale`: a read right after a set (before the
 * microtask) saw the old DOM. Counts are per update. `burstTasks`: 100 sets, each in a task of its own, timed one by
 * one (set and flush), per burst.
 */
async function updateCase(name: UpdateCase, n: number) {
  const zero = () => {
    ReactiveCounts.listener = ReactiveCounts.derive = ReactiveCounts.render = 0;
  };
  const per = (x: number, d = n) => +(x / d).toFixed(4);
  let k0 = Math.floor(Math.random() * 1e6) * 2;
  if (name === 'burstTasks') {
    const c = upd!.fx.components;
    let dom = 0;
    let sync = 0;
    let stale = false;
    zero();
    for (let b = 0; b < n; b++) {
      for (let i = 0; i < 100; i++) {
        const text = `tasks ${k0++}`;
        const t0 = performance.now();
        c[i % 10].label.set(text);
        const t1 = performance.now();
        flush();
        const t2 = performance.now();
        sync += t1 - t0;
        dom += t2 - t0;
        await new Promise((r) => setTimeout(r, 0));
      }
    }
    const counts = { ...ReactiveCounts };
    for (let i = 0; i < 3; i++) {
      const text = `tasks ${k0++}`;
      c[i].label.set(text);
      if (labelOf(`c${i}`) !== text) stale = true;
      flush();
    }
    return {
      syncMs: per(sync),
      domMs: per(dom),
      stale,
      listener: per(counts.listener),
      derive: per(counts.derive),
      render: per(counts.render),
    };
  }
  // A read right after a set, then the DOM after the flush (untimed).
  let stale = false;
  for (let i = 0; i < 3; i++) {
    const u = updateOf(name, k0++);
    u.mutate();
    if (!u.shows()) stale = true;
    flush();
    if (!u.shows()) throw new Error(`${name}: the DOM is wrong after flush`);
  }
  const updates = Array.from({ length: n }, () => updateOf(name, k0++));
  zero();
  let t0 = performance.now();
  for (const u of updates) {
    u.mutate();
    flush();
  }
  const dom = performance.now() - t0;
  const counts = { ...ReactiveCounts };
  if (!updates[n - 1].shows()) throw new Error(`${name}: the DOM is wrong`);
  const again = Array.from({ length: n }, () => updateOf(name, k0++));
  t0 = performance.now();
  for (const u of again) u.mutate();
  const sync = performance.now() - t0;
  flush();
  return {
    syncMs: per(sync),
    domMs: per(dom),
    stale,
    listener: per(counts.listener),
    derive: per(counts.derive),
    render: per(counts.render),
  };
}

/**
 * `n` updates of a case, set and flushed, with nothing else allocated (no closures, no checks): for a driver that reads
 * the heap before and after (allocated per update). Bursts across tasks are not run here.
 */
function updateHeap(name: Exclude<UpdateCase, 'burstTasks'>, n: number, labels: readonly string[]) {
  const { fx } = upd!;
  const c = fx.components;
  for (let k = 0; k < n; k++) {
    switch (name) {
      case 'label':
        c[3].label.set(labels[k]);
        break;
      case 'hidden':
        c[3].hidden.set(k % 2 === 0);
        break;
      case 'swatch':
        c[3].swatch!.set(k % 2 ? 0xff0000 : 0x00ff00);
        break;
      case 'pick':
        fx.scene.pick(c[k % 2 ? 3 : 5].node);
        break;
      case 'reflow':
        c[900].expanded.set(k % 2 === 0);
        break;
      case 'burst':
        for (let i = 0; i < 100; i++) c[i % 10].label.set(labels[(k * 100 + i) % labels.length]);
        break;
    }
    flush();
  }
}

/** Labels made ahead, so a heap reading counts the update, not the strings. */
let heapLabels: string[] = [];

/** The rest's reorder pass (RecycledList): after each of `count` far jumps, rested, its ms and the nodes it moved. */
async function settleCost(count = 20) {
  const t = mountTree('scene');
  await settle(t);
  const vp = t.viewport();
  const span = vp.scrollHeight - vp.clientHeight;
  const ms: number[] = [];
  const moved: number[] = [];
  for (let k = 0; k < count; k++) {
    const before = { ...RecycledListStats };
    vp.scrollTop = Math.round((((k * 37) % 97) / 97) * span);
    vp.dispatchEvent(new Event('scroll'));
    await wait(200);
    if (RecycledListStats.settles > before.settles) {
      ms.push(RecycledListStats.lastMs);
      moved.push(RecycledListStats.moved - before.moved);
    }
  }
  t.dispose();
  return {
    passes: ms.length,
    medianMs: round(pct(ms, 0.5), 3),
    maxMs: round(Math.max(...ms), 3),
    medianMoved: pct(moved, 0.5),
    rows: t.rows().length,
  };
}

(window as any).__mspRows = {
  settleCost,
  /** Makes `n` distinct labels ahead of `updateHeap`. */
  heapLabels(n: number) {
    const base = Math.floor(Math.random() * 1e9);
    heapLabels = Array.from({ length: n }, (_, i) => `heap ${base + i}`);
  },
  updateHeap: (name: Exclude<UpdateCase, 'burstTasks'>, n: number) => updateHeap(name, n, heapLabels),
  /** Mounts the closed tree with swatches for `updateCase`; returns the rows in the window. */
  async updates(height = 600) {
    upd?.dispose();
    mounted?.dispose();
    const host = document.createElement('div');
    host.className = 'msp-root';
    host.style.cssText = `position: fixed; left: 16px; top: 16px; width: 320px; height: ${height}px; contain: strict;`;
    document.body.append(host);
    const box = document.createElement('div');
    box.style.cssText = `height: ${height}px; width: 320px; display: flex; flex-direction: column;`;
    host.append(box);
    const fx = makeScene(1000, 3, true);
    const dispose = render(() => <SceneTree scene={fx.scene} recycle={recycle} />, box);
    upd = {
      fx,
      box,
      dispose: () => {
        dispose();
        host.remove();
      },
    };
    for (let k = 0; k < 300 && box.querySelectorAll('[role=treeitem]').length < 4; k++) await frame();
    await wait(200);
    await frame();
    return box.querySelectorAll('[role=treeitem]').length;
  },
  updateCase,
  jumps,
  steps,
  /** Mounts a tall tree (`height` px) and leaves it; returns the rows in the DOM once settled. */
  async mount(kind: TreeKind, height = 8000) {
    mounted?.dispose();
    mounted = mountTree(kind, height);
    await settle(mounted);
    return mounted.rows().length;
  },
  unmount() {
    mounted?.dispose();
    mounted = undefined;
    upd?.dispose();
    upd = undefined;
  },
  /** Mounts the normal tree for a wheel test by a driver: the viewport's center in page coordinates. */
  async wheelTarget(kind: TreeKind) {
    mounted?.dispose();
    mounted = mountTree(kind);
    await settle(mounted);
    const r = mounted.viewport().getBoundingClientRect();
    return { x: r.left + r.width / 2, y: r.top + r.height / 2 };
  },
  /** Starts recording frames (rAF intervals and long animation frames); `stop` returns them. */
  record() {
    const intervals: number[] = [];
    const loaf: number[] = [];
    let uiScript = 0;
    let otherScript = 0;
    let last = performance.now();
    let on = true;
    const tick = (now: number) => {
      intervals.push(now - last);
      last = now;
      if (on) requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
    const observer = PerformanceObserver.supportedEntryTypes?.includes('long-animation-frame')
      ? new PerformanceObserver((list) => {
          for (const e of list.getEntries() as any[]) {
            loaf.push(e.duration);
            for (const sc of e.scripts ?? []) {
              const url = sc.sourceURL ?? '';
              if (/\/packages\/(prism|sequences)\//.test(url) && !/_stories/.test(url)) uiScript += sc.duration;
              else otherScript += sc.duration;
            }
          }
        })
      : undefined;
    observer?.observe({ type: 'long-animation-frame' });
    return () => {
      on = false;
      observer?.disconnect();
      return {
        frames: intervals.length,
        over20ms: intervals.filter((x) => x > 20).length,
        maxFrame: round(Math.max(...intervals), 1),
        longFrames: loaf.length,
        longFramesMs: round(
          loaf.reduce((a, b) => a + b, 0),
          1,
        ),
        longFramesUiScriptMs: round(uiScript, 1),
        longFramesOtherScriptMs: round(otherScript, 1),
        rowsInDom: mounted?.rows().length,
      };
    };
  },
};

/** The probes above (`window.__mspRows`); the story itself is a note. Run them from a driver or the console. */
export const RowsMade: StoryObj = {
  name: 'Rows made',
  parameters: noA11y,
  render: () =>
    story(() => (
      <p style="padding: var(--msp-space-5); color: var(--msp-color-text);">
        Probes: <code>await __mspRows.jumps('scene')</code>, <code>__mspRows.steps('state')</code>,{' '}
        <code>__mspRows.mount('scene')</code>. They mount outside Storybook's root.
      </p>
    )),
};
