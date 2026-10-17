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
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import type { StateTreeModel, StateTreeNode } from '../model/state-tree.js';
import { type ObservableState, observableState } from '@molstar/prism-ui/reactive/observable-state';
import { flattenTree, VirtualList } from '@molstar/prism-ui';
import { StateTreeTab } from '../views/state-tree/state-tree-tab.js';
import { SceneTree } from '../views/structure/scene-tree.js';
import { makeScene } from './performance-fixtures.js';
import { story } from '@molstar/prism-ui/testing/index';

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
        aria-selected={selected.asReadonly.pipe(
          map((s) => s === props.ref),
          distinctUntilChanged(),
        )}
        style={`display: flex; gap: var(--msp-space-3); height: 26px; align-items: center; padding-left: ${props.depth * 14 + 4}px`}
      >
        <span>
          {record.pipe(
            map((n) => n.label),
            distinctUntilChanged(),
          )}
        </span>
        <span style="color: var(--msp-color-text-subtle)">
          {record.pipe(
            map((n) => n.description ?? ''),
            distinctUntilChanged(),
          )}
        </span>
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
          {running.asReadonly.pipe(map((r) => (r ? 'Running…' : 'Run')))}
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
 * and open, Enter to select). The readout counts the rows in the DOM against the rows of the tree. Mounted outside
 * Storybook's root and with its accessibility checks off: its addons watch every DOM change under the root (axe runs
 * again, the highlight addon re-scans styles, ~13 ms a time), so a scrollbar drag, which changes rows every frame, would
 * time them, and more the longer it goes.
 */
export const Interactive: StoryObj = {
  name: 'Virtualized trees (interactive)',
  parameters: noA11y,
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
      const content = () => interactivePanels(fixture, model, cells, panel, mount);
      const host = document.createElement('div');
      host.className = 'msp-root';
      host.style.cssText = 'position: fixed; inset: 0; overflow: auto; background: var(--msp-color-bg);';
      document.body.append(host);
      const dispose = render(() => content(), host);
      // A/B for the skeleton (virtual-list.module.css): off, the lists' tall content boxes are not painted, so a slow
      // scroll that is the raster of those boxes (not script) shows as the difference between the two.
      const off = document.createElement('style');
      off.textContent = '.msp-virtual-list-content { background-image: none !important; }';
      const toggle = document.createElement('label');
      toggle.style.cssText =
        'position: fixed; right: 16px; top: 12px; display: flex; gap: 6px; align-items: center; font: 13px system-ui; color: var(--msp-color-text);';
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.checked = true;
      box.addEventListener('change', () => (box.checked ? off.remove() : document.head.append(off)));
      toggle.append(box, 'Skeleton behind the rows');
      host.append(toggle);
      onCleanup(() => {
        dispose();
        off.remove();
        host.remove();
      });
      return <p>The trees are mounted outside Storybook's root (its addons watch every change there).</p>;
    }),
};

function interactivePanels(
  fixture: ReturnType<typeof makeScene>,
  model: ReturnType<typeof treeModel>,
  cells: number,
  panel: (title: string, body: () => HTMLElement) => HTMLElement,
  mount: (view: () => Node) => () => HTMLElement,
) {
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
}

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
