/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { distinctUntilChanged, map, type Observable } from 'rxjs';
import { For } from '../dom/control-flow.js';
import { render } from '../dom/render.js';
import type { StateTreeModel, StateTreeNode } from '../model/state-tree.js';
import { type ObservableState, observableState } from '../reactive/observable-state.js';
import { flattenTree, VirtualList } from '../components/virtual-list/virtual-list.js';
import { StateTreeTab } from '../views/state-tree/state-tree-tab.js';
import { story } from './story-root.js';

const meta: Meta = { title: 'Internals/Performance' };
export default meta;

// The renderer check (plan.md Phase 2): the state tree at 5,000 cells under bursts of updates, in the `dom` runtime.
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
  const all = () => box.querySelectorAll('[role=treeitem]').length >= cells;
  out.push(await time(`State tree tab: first render, ${cells} cells`, () => roots.set(tree), host, undefined, all));
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
      'For re-creates an entry when its record changes: every touched root subtree is rebuilt',
      all,
    ),
  );
  const collapseAll = box.querySelectorAll<HTMLButtonElement>('[aria-label=Expand]');
  out.push(
    await time(
      'State tree tab: collapse and expand the 25 structures',
      () => {
        for (let i = 0; i < collapseAll.length && i < 5000; i++) {
          const b = collapseAll[i];
          if (b.closest('[role=treeitem]')?.parentElement?.getAttribute('role') === 'tree') b.click();
        }
        for (const b of box.querySelectorAll<HTMLButtonElement>('[aria-label=Expand]')) {
          if (b.closest('[role=treeitem]')?.parentElement?.getAttribute('role') === 'tree') b.click();
        }
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
    <div style={`display: flex; gap: 8px; height: 26px; align-items: center; padding-left: ${props.row.depth * 12}px`}>
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
        style={`display: flex; gap: 8px; height: 26px; align-items: center; padding-left: ${props.depth * 14 + 4}px`}
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

/**
 * The renderer check: the state tree at 5,025 cells in the `dom` runtime. Run it, then read the table (wall time in
 * this browser, layout forced).
 */
export const StateTree: StoryObj = {
  // Storybook's accessibility checks run on every change in the story: at 5,000 rows they, not the UI, would be timed.
  parameters: { a11y: { disable: true, test: 'off' } },
  render: () =>
    story(() => {
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
        await realTab(host, out);
        await stableRows(host, out);
        await virtualTree(host, out);
        host.replaceChildren();
        host.remove();
        results.set(out);
        running.set(false);
        (window as any).__mspPerf = { hidden: document.hidden, results: out, done: true };
      };
      return (
        <div style="display: grid; gap: 12px; padding: 16px; font-size: var(--msp-font-size-sm); color: var(--msp-color-text);">
          <button
            type="button"
            style="justify-self: start; height: 28px; padding: 0 12px; border-radius: 6px; background: var(--msp-color-control);"
            disabled={running.asReadonly}
            onClick={run}
          >
            {running.asReadonly.pipe(map((r) => (r ? 'Running…' : 'Run')))}
          </button>
          <table style="border-collapse: collapse;">
            <For each={results.asReadonly}>
              {(r) => (
                <tr>
                  <td style="padding: 2px 12px 2px 0">{r.case}</td>
                  <td style="padding: 2px 12px; text-align: right; font-variant-numeric: tabular-nums">{r.ms} ms</td>
                  <td style="padding: 2px 12px; text-align: right; font-variant-numeric: tabular-nums">
                    {r.longestTask ? `longest task ${r.longestTask} ms` : 'no long task'}
                  </td>
                  <td style="padding: 2px 0; color: var(--msp-color-text-subtle)">{r.note ?? ''}</td>
                </tr>
              )}
            </For>
          </table>
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
