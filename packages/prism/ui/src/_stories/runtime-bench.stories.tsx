/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { map, type Observable } from 'rxjs';
import { For, timeSlicing } from '../dom/control-flow.js';
import { render } from '../dom/render.js';
import { keyedSelector } from '../reactive/keyed-selector.js';
import { batch, type ObservableState, observableState } from '../reactive/observable-state.js';
import { story } from '../testing/story-root.js';

const meta: Meta = { title: 'Internals/Performance' };
export default meta;

// The renderer bench (spike.md, Next): the js-framework-benchmark operations on the `dom` runtime, rows written the way
// the app writes them (a record per row, keyed `For`, the selection a reactive class). Script is the operation alone;
// total adds forced layout and a frame. Medians of `Runs` after a warm-up; time slicing off (raw throughput).

const Runs = 5;

const adjectives = ['pretty', 'large', 'big', 'small', 'tall', 'short', 'long', 'handsome', 'plain', 'quaint', 'clean'];
const colours = ['red', 'yellow', 'blue', 'green', 'pink', 'brown', 'purple', 'orange', 'white', 'black', 'grey'];
const nouns = ['table', 'chair', 'house', 'bbq', 'desk', 'car', 'pony', 'cookie', 'sandwich', 'burger', 'pizza'];

interface RowData {
  id: number;
  label: ObservableState<string>;
}

let nextId = 1;
let seed = 1;
// A fixed sequence (no Math.random: runs compare).
const pick = (list: string[]) => {
  seed = (seed * 16807) % 2147483647;
  return list[seed % list.length];
};
function buildRows(count: number): RowData[] {
  const rows: RowData[] = new Array(count);
  for (let i = 0; i < count; i++)
    rows[i] = { id: nextId++, label: observableState(`${pick(adjectives)} ${pick(colours)} ${pick(nouns)}`) };
  return rows;
}

/** The same row without bindings (a floor: element creation and the list alone). */
function StaticRow(props: { row: RowData }) {
  return (
    <tr>
      <td class="col-md-1">{props.row.id}</td>
      <td class="col-md-4">
        <a>{props.row.label.value}</a>
      </td>
      <td class="col-md-1">
        <a>×</a>
      </td>
      <td class="col-md-6" />
    </tr>
  );
}

function Row(props: {
  row: RowData;
  selected: { is(id: number): Observable<boolean> };
  select(id: number): void;
  remove(id: number): void;
}) {
  const { row } = props;
  return (
    <tr classList={{ danger: props.selected.is(row.id) }}>
      <td class="col-md-1">{row.id}</td>
      <td class="col-md-4">
        <a onClick={() => props.select(row.id)}>{row.label.asReadonly}</a>
      </td>
      <td class="col-md-1">
        <a onClick={() => props.remove(row.id)}>×</a>
      </td>
      <td class="col-md-6" />
    </tr>
  );
}

interface Bench {
  rows: ObservableState<readonly RowData[]>;
  selected: ObservableState<number | undefined>;
  table: HTMLElement;
  dispose(): void;
}

function mountBench(host: HTMLElement): Bench {
  const rows = observableState<readonly RowData[]>([]);
  const selected = observableState<number | undefined>(undefined);
  // The selected row as the app should do it: only the two rows that change hear of it.
  const isSelected = keyedSelector(selected.asReadonly);
  const box = document.createElement('div');
  box.style.cssText = 'height: 100%; overflow: auto; font-size: 12px;';
  host.replaceChildren(box);
  const dispose = render(
    () => (
      <table>
        <tbody>
          <For each={rows.asReadonly} key={(r) => r.id}>
            {(row) =>
              (window as any).__mspStaticRows ? (
                <StaticRow row={row} />
              ) : (
                <Row
                  row={row}
                  selected={isSelected}
                  select={(id) => selected.set(id)}
                  remove={(id) => rows.set(rows.value.filter((r) => r.id !== id))}
                />
              )
            }
          </For>
        </tbody>
      </table>
    ),
    box,
  );
  return { rows, selected, table: box, dispose };
}

// Two frames when the page shows; a hidden page gets none (script is still measured, total is then not meaningful).
const frame = () =>
  document.hidden
    ? Promise.resolve()
    : new Promise<void>((r) => requestAnimationFrame(() => requestAnimationFrame(() => r())));
const median = (xs: number[]) => {
  const s = xs.slice().sort((a, b) => a - b);
  return s[Math.floor(s.length / 2)];
};

interface Op {
  name: string;
  /** Untimed: brings the bench to the state the operation starts from. */
  setup(b: Bench): void;
  run(b: Bench): void;
}

const Ops: Op[] = [
  { name: 'create 1,000 rows', setup: (b) => b.rows.set([]), run: (b) => b.rows.set(buildRows(1000)) },
  {
    name: 'replace all 1,000 rows',
    setup: (b) => b.rows.set(buildRows(1000)),
    run: (b) => b.rows.set(buildRows(1000)),
  },
  {
    name: 'update every 10th row (of 1,000)',
    setup: (b) => b.rows.set(buildRows(1000)),
    run: (b) => {
      const rows = b.rows.value;
      for (let i = 0; i < rows.length; i += 10) rows[i].label.set(`${rows[i].label.value} !!!`);
    },
  },
  {
    name: 'select a row (of 1,000)',
    setup: (b) => {
      b.rows.set(buildRows(1000));
      b.selected.set(undefined);
    },
    run: (b) => b.selected.set(b.rows.value[500].id),
  },
  {
    name: 'swap rows 2 and 999',
    setup: (b) => b.rows.set(buildRows(1000)),
    run: (b) => {
      const rows = b.rows.value.slice();
      const t = rows[1];
      rows[1] = rows[998];
      rows[998] = t;
      b.rows.set(rows);
    },
  },
  {
    name: 'remove a row (of 1,000)',
    setup: (b) => b.rows.set(buildRows(1000)),
    run: (b) => b.rows.set(b.rows.value.filter((_, i) => i !== 500)),
  },
  { name: 'create 10,000 rows', setup: (b) => b.rows.set([]), run: (b) => b.rows.set(buildRows(10000)) },
  {
    name: 'append 1,000 rows (to 10,000)',
    setup: (b) => b.rows.set(buildRows(10000)),
    run: (b) => b.rows.set([...b.rows.value, ...buildRows(1000)]),
  },
  { name: 'clear 10,000 rows', setup: (b) => b.rows.set(buildRows(10000)), run: (b) => b.rows.set([]) },
  // A burst: one change made of many steps (a model relabelling as it goes), without and with `batch`.
  {
    name: 'relabel all 1,000 rows 10 times',
    setup: (b) => b.rows.set(buildRows(1000)),
    run: (b) => relabelTenTimes(b.rows.value),
  },
  {
    name: 'relabel all 1,000 rows 10 times, batched',
    setup: (b) => b.rows.set(buildRows(1000)),
    run: (b) => batch(() => relabelTenTimes(b.rows.value)),
  },
];

function relabelTenTimes(rows: readonly RowData[]) {
  for (let k = 0; k < 10; k++) for (const r of rows) r.label.set(`${r.id} step ${k}`);
}

interface Result {
  op: string;
  script: number;
  total: number;
}

async function runBench(host: HTMLElement, only?: string): Promise<Result[]> {
  timeSlicing(false);
  const out: Result[] = [];
  try {
    for (const op of Ops) {
      if (only && !op.name.includes(only)) continue;
      const scripts: number[] = [];
      const totals: number[] = [];
      for (let i = 0; i <= Runs; i++) {
        const b = mountBench(host);
        op.setup(b);
        void host.offsetHeight;
        await frame();
        const t0 = performance.now();
        op.run(b);
        const t1 = performance.now();
        void host.offsetHeight;
        await frame();
        const t2 = performance.now();
        // The first run warms up.
        if (i > 0) {
          scripts.push(t1 - t0);
          totals.push(t2 - t0);
        }
        b.dispose();
      }
      out.push({ op: op.name, script: +median(scripts).toFixed(1), total: +median(totals).toFixed(1) });
      (window as any).__mspBench = { results: out.slice() };
    }
  } finally {
    timeSlicing('auto');
  }
  return out;
}

/**
 * The renderer bench: the js-framework-benchmark operations on the `dom` runtime (Run, then read the table). Script is
 * the operation's own time; total adds layout and a frame (two frames are ~33 ms of it).
 */
export const Runtime: StoryObj = {
  parameters: { a11y: { disable: true, test: 'off' } },
  render: () =>
    story(() => {
      const results = observableState<readonly Result[]>([]);
      const running = observableState(false);
      // Outside Storybook's root (its addons watch every change there), a fixed, contained area.
      const host = document.createElement('div');
      host.className = 'msp-root';
      host.style.cssText = 'position: fixed; right: 16px; top: 16px; width: 480px; height: 600px; contain: strict;';
      const run = async () => {
        running.set(true);
        document.body.append(host);
        const out = await runBench(host);
        host.remove();
        results.set(out);
        running.set(false);
        (window as any).__mspBench = { results: out, done: true };
      };
      (window as any).__mspBenchRun = async (only?: string) => {
        document.body.append(host);
        const out = await runBench(host, only);
        host.remove();
        return out;
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
                  <td style="padding: 2px 12px 2px 0">{r.op}</td>
                  <td style="padding: 2px 12px; text-align: right; font-variant-numeric: tabular-nums">
                    script {r.script} ms
                  </td>
                  <td style="padding: 2px 0; text-align: right; font-variant-numeric: tabular-nums">
                    total {r.total} ms
                  </td>
                </tr>
              )}
            </For>
          </table>
        </div>
      );
    }),
};
