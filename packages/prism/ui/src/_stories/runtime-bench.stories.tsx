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
import { flush, watch } from '../dom/scope.js';
import { setClass, setText } from '../dom/writes.js';
import { keyedSelector } from '../reactive/keyed-selector.js';
import { batch, type ObservableState, observableState } from '../reactive/observable-state.js';
import { story } from '../testing/story-root.js';
import { CommandManager } from '../commands/manager.js';
import type { CommandDefinition } from '../commands/types.js';
import { CommandMenu } from '../views/command-menu/command-menu.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsForm } from '../components/params/params-form.js';
import { ParamsFormModel } from '../forms/params-form.js';

const meta: Meta = { title: 'Internals/Runtime performance' };
export default meta;

// The renderer bench (rendering.md §1): the js-framework-benchmark operations on the `dom` runtime, rows written the
// way the app writes them (a record per row, keyed `For`, the selection a reactive class). Script is the operation
// alone; total adds forced layout and a frame. Medians of `Runs` after a warm-up; time slicing off (raw throughput).

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

/**
 * The row as view = f(state) (rendering.md §2): its DOM made once, one `watch` render that writes what differs. The
 * selection is read whole (every row renders when it changes, and writes only if its own flag flipped).
 */
function WatchRow(props: {
  row: RowData;
  selected: ObservableState<number | undefined>;
  select(id: number): void;
  remove(id: number): void;
}) {
  const { row, selected } = props;
  const label = document.createTextNode('');
  const tr = (
    <tr>
      <td class="col-md-1">{row.id}</td>
      <td class="col-md-4">
        <a onClick={() => props.select(row.id)}>{label}</a>
      </td>
      <td class="col-md-1">
        <a onClick={() => props.remove(row.id)}>×</a>
      </td>
      <td class="col-md-6" />
    </tr>
  ) as HTMLElement;
  watch([row.label, selected], () => {
    setText(label, row.label.value);
    setClass(tr, 'danger', selected.value === row.id);
  });
  return tr;
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
              ) : (window as any).__mspWatchRows ? (
                <WatchRow
                  row={row}
                  selected={selected}
                  select={(id) => selected.set(id)}
                  remove={(id) => rows.set(rows.value.filter((r) => r.id !== id))}
                />
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
        flush();
        void host.offsetHeight;
        await frame();
        const t0 = performance.now();
        op.run(b);
        // Watched rows render in a microtask: in the operation's time.
        flush();
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
                  <td style="padding: var(--msp-space-1) var(--msp-space-4) var(--msp-space-1) 0">{r.op}</td>
                  <td style="padding: var(--msp-space-1) var(--msp-space-4); text-align: right; font-variant-numeric: tabular-nums">
                    script {r.script} ms
                  </td>
                  <td style="padding: var(--msp-space-1) 0; text-align: right; font-variant-numeric: tabular-nums">
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

// --- The command menu's rows ---------------------------------------------------------------------

/**
 * Probes for a driver (`window.__mspMenu`; Playwright adds garbage collection and the heap): the command menu with
 * `n` commands in 20 categories, mounted outside Storybook's root (its highlight addon re-scans styles on every change
 * there). `open` times the first render; `grow` scrolls it to the end, rendering every row; `type` times keystrokes
 * (each a new list of rows); `moves` times the highlight moving down.
 */
let menu: { host: HTMLElement; dispose(): void } | undefined;
(window as any).__mspMenu = {
  open(n = 1000) {
    menu?.dispose();
    const manager = new CommandManager<object>({
      context: observableState<object>({}).asReadonly,
      formats: { get: () => undefined },
    });
    const commands: CommandDefinition<object>[] = [];
    for (let i = 0; i < n; i++)
      commands.push({
        id: `cmd.${i}`,
        title: `${pick(adjectives)} ${pick(colours)} ${pick(nouns)} ${i}`,
        category: `Category ${i % 20}`,
        icon: 'cartoon',
        run: () => ({ message: '' }),
      });
    manager.register(commands);
    const host = document.createElement('div');
    host.className = 'msp-root';
    host.style.cssText = 'position: fixed; right: 16px; top: 16px; width: 480px; height: 600px; contain: strict;';
    document.body.append(host);
    const t0 = performance.now();
    const dispose = render(() => <CommandMenu manager={manager} onClose={() => {}} />, host);
    flush();
    const ms = performance.now() - t0;
    menu = {
      host,
      dispose: () => {
        dispose();
        host.remove();
      },
    };
    return { ms: +ms.toFixed(2), rows: host.querySelectorAll('[role=option]').length };
  },
  /** Scrolls the list to its end until every row is rendered. */
  async grow() {
    const vp = menu!.host.querySelector<HTMLElement>('[role=listbox]')!.closest('[data-axis]')!
      .firstElementChild as HTMLElement;
    for (let k = 0; k < 400; k++) {
      const before = menu!.host.querySelectorAll('[role=option]').length;
      vp.scrollTop = vp.scrollHeight;
      vp.dispatchEvent(new Event('scroll'));
      flush();
      await frame();
      if (menu!.host.querySelectorAll('[role=option]').length === before) break;
    }
    return menu!.host.querySelectorAll('[role=option]').length;
  },
  type(text: string) {
    const input = menu!.host.querySelector('input')!;
    const times: number[] = [];
    for (let i = 1; i <= text.length; i++) {
      const t0 = performance.now();
      input.value = text.slice(0, i);
      input.dispatchEvent(new Event('input', { bubbles: true }));
      flush();
      times.push(performance.now() - t0);
    }
    return {
      msPerKey: +(times.reduce((a, b) => a + b, 0) / times.length).toFixed(2),
      rows: menu!.host.querySelectorAll('[role=option]').length,
    };
  },
  moves(k = 50) {
    const input = menu!.host.querySelector('input')!;
    const t0 = performance.now();
    for (let i = 0; i < k; i++) {
      input.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowDown', bubbles: true }));
      flush();
    }
    return { msPerMove: +((performance.now() - t0) / k).toFixed(3) };
  },
  close() {
    menu?.dispose();
    menu = undefined;
  },
};

// --- The params form's fields --------------------------------------------------------------------

/**
 * Probes for a driver (`window.__mspForm`): a flat form of `n` fields (sliders, numbers, switches, texts, colors in
 * turn), outside Storybook's root. `open` times the first render; `setAll` sets every field in one batch (a preset);
 * `setOne` one field, `count` times (a slider dragged).
 */
let form: { model: ParamsFormModel<PD.Params>; host: HTMLElement; dispose(): void } | undefined;
(window as any).__mspForm = {
  open(n = 200) {
    form?.dispose();
    const params: PD.Params = {};
    for (let i = 0; i < n; i++) {
      const k = i % 5;
      params[`f${i}`] =
        k === 0
          ? PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 })
          : k === 1
            ? PD.Numeric(3)
            : k === 2
              ? PD.Boolean(true)
              : k === 3
                ? PD.Text(`text ${i}`)
                : PD.Color(0x2e86c1 as any);
    }
    const model = new ParamsFormModel(params);
    const host = document.createElement('div');
    host.className = 'msp-root';
    host.style.cssText = 'position: fixed; right: 16px; top: 16px; width: 360px; height: 600px; overflow: auto;';
    document.body.append(host);
    const t0 = performance.now();
    const dispose = render(() => <ParamsForm form={model} />, host);
    flush();
    const ms = performance.now() - t0;
    form = {
      model,
      host,
      dispose: () => {
        dispose();
        host.remove();
      },
    };
    return { ms: +ms.toFixed(2), fields: n };
  },
  setAll(round = 1) {
    const m = form!.model;
    const values: Record<string, unknown> = {};
    for (const [k, p] of Object.entries(m.params)) {
      values[k] =
        p.type === 'number'
          ? (p.max ?? 10) * ((round % 9) / 10)
          : p.type === 'boolean'
            ? round % 2 === 0
            : p.type === 'text'
              ? `text ${k} ${round}`
              : p.type === 'color'
                ? 0x100000 * (round % 15)
                : m.values.value[k];
    }
    const t0 = performance.now();
    m.load(values as any);
    flush();
    return +(performance.now() - t0).toFixed(2);
  },
  setOne(count = 100) {
    const m = form!.model;
    const t0 = performance.now();
    for (let i = 0; i < count; i++) {
      m.set('f0', (i % 100) / 100);
      flush();
    }
    return +((performance.now() - t0) / count).toFixed(3);
  },
  close() {
    form?.dispose();
    form = undefined;
  },
};
