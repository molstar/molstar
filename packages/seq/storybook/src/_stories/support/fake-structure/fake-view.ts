/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Cell, Span } from './chain-table.js';
import { frameBox } from './frame-box.js';
import type { FakeStructureViewer } from './fake-viewer.js';

// The fake structure viewer's view (a thin one: the model is `fake-viewer.ts`): a chain instance a row, a residue a
// cell (unobserved ones dashed gaps, a ligand one chip), the state's colors, the pointer's hover as an outline, what
// the host lit as a ring, the selection as a bar under the cells, and the camera as a frame that moves over the time
// it was given. Its grammar is the model's: hover; click selects; Cmd/Ctrl+click toggles; Esc clears; double-click
// focuses; the two buttons make selections of their own.

const CellWidth = 4;
const RowHeight = 18;
const LabelWidth = 76;

const Css = `
.fv { display: grid; gap: var(--msui-space-3); font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-sans); color: var(--msui-text); outline: none; }
.fv-tools { display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap; }
.fv-status { color: var(--msui-text-muted); font-family: var(--msui-font-mono); }
.fv-rows { position: relative; display: grid; gap: 4px; overflow-x: auto; padding: 6px; background: var(--msui-surface); border: 1px solid var(--msui-line); border-radius: var(--msui-radius-md); }
.fv-row { display: flex; align-items: center; height: ${RowHeight}px; }
.fv-label { flex: none; width: ${LabelWidth}px; color: var(--msui-text-muted); font-family: var(--msui-font-mono); }
.fv-cells { display: flex; position: relative; }
.fv-cell { box-sizing: border-box; width: ${CellWidth}px; height: ${RowHeight - 4}px; flex: none; }
.fv-cell[data-gap] { background: transparent !important; border: 1px dashed var(--msui-line); border-width: 1px 0; }
.fv-cell[data-chip] { width: auto; min-width: 30px; padding: 0 6px; border-radius: 9px; font: var(--msui-text-xs) / ${RowHeight - 4}px var(--msui-font-mono); color: var(--msui-bg); text-align: center; }
.fv-cell[data-hover] { box-shadow: 0 0 0 2px var(--msui-text); position: relative; z-index: 2; }
.fv-cell[data-lit] { box-shadow: 0 0 0 2px var(--msui-accent, #d6336c); position: relative; z-index: 1; }
.fv-cell[data-selected] { border-bottom: 3px solid var(--msui-text); }
.fv-frame { position: absolute; pointer-events: none; box-sizing: border-box; border: 2px dashed var(--msui-accent, #d6336c); border-radius: 4px; z-index: 3; }
`;

/** The cell a pointer event is on. */
function cellOf(e: Event): Cell | null {
  const el = (e.target as HTMLElement | null)?.closest<HTMLElement>('.fv-cell');
  return el ? { row: Number(el.dataset.row), index: Number(el.dataset.index) } : null;
}

const within = (spans: readonly Span[], row: number, index: number) =>
  spans.some((s) => s.row === row && index >= s.start && index < s.end);

class FakeStructureElement extends UiElement<FakeStructureViewer> {
  static readonly tag = 'seqx-fake-structure';

  protected render(fake: FakeStructureViewer): HTMLElement {
    const root = document.createElement('div');
    root.className = 'fv';
    root.tabIndex = 0;
    const style = document.createElement('style');
    style.textContent = Css;

    const tools = document.createElement('div');
    tools.className = 'fv-tools';
    const status = document.createElement('span');
    status.className = 'fv-status';
    tools.append(
      Button({ children: 'Select chain', onClick: () => fake.selectChain() }) as HTMLElement,
      Button({ children: 'Within 5 Å', onClick: () => fake.within() }) as HTMLElement,
      Button({ children: 'Clear', onClick: () => fake.clear() }) as HTMLElement,
      status,
    );

    const rowsEl = document.createElement('div');
    rowsEl.className = 'fv-rows';
    const cells: HTMLElement[][] = [];
    for (const row of fake.rows) {
      const r = document.createElement('div');
      r.className = 'fv-row';
      const label = document.createElement('span');
      label.className = 'fv-label';
      label.textContent = `${row.asym} ${row.operator}`;
      label.title = `label_asym_id ${row.asym}, auth_asym_id ${row.auth}, entity ${row.entity}, operator ${row.operator}`;
      const line = document.createElement('div');
      line.className = 'fv-cells';
      const els: HTMLElement[] = [];
      row.residues.forEach((res, index) => {
        const c = document.createElement('span');
        c.className = 'fv-cell';
        c.dataset.row = String(row.index);
        c.dataset.index = String(index);
        c.title = `${row.asym} ${res.auth}${res.ins} ${res.comp}${res.observed ? '' : ' (not modeled)'}`;
        if (!res.observed) c.dataset.gap = '';
        if (!row.polymer) {
          c.dataset.chip = '';
          c.textContent = res.comp;
        }
        els.push(c);
        line.append(c);
      });
      cells.push(els);
      r.append(label, line);
      rowsEl.append(r);
    }
    const frame = document.createElement('div');
    frame.className = 'fv-frame';
    rowsEl.append(frame);
    root.append(style, tools, rowsEl);

    const draw = () => {
      const hover = fake.hover.value;
      const lit = fake.highlight.value;
      const sel = fake.selection.value;
      fake.rows.forEach((row) => {
        cells[row.index].forEach((c, i) => {
          if (!row.residues[i].observed) return;
          c.style.background = fake.colorOf(row.index, i);
          c.toggleAttribute('data-hover', within(hover, row.index, i));
          c.toggleAttribute('data-lit', within(lit, row.index, i));
          c.toggleAttribute('data-selected', within(sel, row.index, i));
        });
      });
      status.textContent = `${fake.entry}: ${fake.rows.length} rows, ${sel.reduce((n, s) => n + s.end - s.start, 0)} selected`;
    };
    const place = () => {
      const cam = fake.camera.value;
      if (!cam) {
        frame.style.display = 'none';
        return;
      }
      const rect = (el: HTMLElement) => el.getBoundingClientRect();
      const boxes = cam.spans.flatMap((s) => [rect(cells[s.row][s.start]), rect(cells[s.row][s.end - 1])]);
      const { left, top, width, height } = frameBox(
        boxes,
        rect(rowsEl),
        { left: rowsEl.scrollLeft, top: rowsEl.scrollTop },
        6,
      );
      frame.style.transition = `left ${cam.durationMs}ms, top ${cam.durationMs}ms, width ${cam.durationMs}ms, height ${cam.durationMs}ms`;
      frame.style.display = 'block';
      frame.style.left = `${left}px`;
      frame.style.top = `${top}px`;
      frame.style.width = `${width}px`;
      frame.style.height = `${height}px`;
    };
    draw();
    place();
    for (const s of [fake.hover, fake.highlight, fake.selection, fake.mvs]) {
      const off = s.listen(draw);
      onCleanup(() => off());
    }
    const offCamera = fake.camera.listen(place);
    onCleanup(() => offCamera());

    rowsEl.addEventListener('pointerover', (e) => fake.pointerOver(cellOf(e)));
    rowsEl.addEventListener('pointerleave', () => fake.pointerOver(null));
    rowsEl.addEventListener('click', (e) => {
      const c = cellOf(e);
      if (c) fake.click(c, e.metaKey || e.ctrlKey);
    });
    rowsEl.addEventListener('dblclick', (e) => {
      const c = cellOf(e);
      if (c) fake.doubleClick(c);
    });
    root.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') fake.clear();
    });
    return root;
  }
}

defineElements([FakeStructureElement]);

/** The element of a fake structure viewer (the story places it; the viewer is the story's to dispose). */
export function fakeStructureElement(fake: FakeStructureViewer): HTMLElement {
  const el = document.createElement(FakeStructureElement.tag) as FakeStructureElement;
  el.model = fake;
  return el;
}
