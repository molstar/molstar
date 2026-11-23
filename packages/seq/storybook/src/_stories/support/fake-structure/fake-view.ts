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
// focuses; the buttons make selections of their own, and "Drop chain" edits the table (the chain of the selection or the
// pointer, else the last polymer's, with its ligands).

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
.fv-card { position: absolute; display: none; pointer-events: none; z-index: 4; padding: 2px 6px; white-space: nowrap; font-family: var(--msui-font-mono); box-shadow: 0 0 0 1px var(--msui-line), var(--msui-shadow-float); }
.fv-note { margin: 4px 0; width: fit-content; padding: 2px 6px; font-family: var(--msui-font-mono); box-shadow: 0 0 0 1px var(--msui-line); }
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
      Button({
        children: 'Drop chain',
        onClick: () => {
          const auth = fake.chainToDrop();
          if (auth !== undefined) fake.dropChain(auth);
        },
      }) as HTMLElement,
      status,
    );

    const rowsEl = document.createElement('div');
    rowsEl.className = 'fv-rows';
    const frame = document.createElement('div');
    frame.className = 'fv-frame';
    const card = document.createElement('div');
    card.className = 'fv-card msui-surface';
    card.setAttribute('aria-hidden', 'true');
    // the card of "Within 5 Å": how many residues, of which chain's neighbors
    const note = document.createElement('div');
    note.className = 'fv-note msui-surface';
    note.setAttribute('role', 'status');
    note.style.display = 'none';
    let cells: HTMLElement[][] = [];
    // the rows as the table has them now (an edit builds them again)
    const buildRows = () => {
      cells = [];
      const children: HTMLElement[] = [];
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
        children.push(r);
      }
      rowsEl.replaceChildren(...children, frame, card);
    };
    buildRows();
    root.append(style, tools, note, rowsEl);

    const draw = () => {
      const hover = fake.hover.value;
      const lit = fake.highlight.value;
      const sel = fake.selection.value;
      fake.rows.forEach((row) => {
        // rows an edit made are built a moment later (the edit's own listener): until then, nothing to draw on
        if (cells[row.index]?.length !== row.residues.length) return;
        cells[row.index].forEach((c, i) => {
          if (!row.residues[i].observed) return;
          c.style.background = fake.colorOf(row.index, i);
          c.toggleAttribute('data-hover', within(hover, row.index, i));
          c.toggleAttribute('data-lit', within(lit, row.index, i));
          c.toggleAttribute('data-selected', within(sel, row.index, i));
        });
      });
      // the card names the hovered residue, whoever hovered it (the pointer here, the host's: a first user looked for names)
      const name = fake.hoverName();
      const at = fake.hover.value[0];
      const cell = at ? cells[at.row]?.[at.start] : undefined;
      if (name && cell) {
        const box = cell.getBoundingClientRect();
        const host = rowsEl.getBoundingClientRect();
        card.textContent = name;
        card.style.display = 'block';
        card.style.left = `${box.left - host.left + rowsEl.scrollLeft}px`;
        card.style.top = `${box.bottom - host.top + rowsEl.scrollTop + 4}px`;
      } else card.style.display = 'none';
      status.textContent = `${fake.entry}: ${fake.rows.length} rows, ${sel.reduce((n, s) => n + s.end - s.start, 0)} selected`;
    };
    const drawNote = () => {
      const text = fake.note.value;
      note.textContent = text ?? '';
      note.style.display = text ? 'block' : 'none';
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
    // an edit of the table: the rows built again, then drawn
    const offEdits = fake.edits.listen(() => {
      buildRows();
      draw();
      place();
    });
    onCleanup(() => offEdits());
    for (const s of [fake.hover, fake.highlight, fake.selection, fake.mvs]) {
      const off = s.listen(draw);
      onCleanup(() => off());
    }
    drawNote();
    const offNote = fake.note.listen(drawNote);
    onCleanup(() => offNote());
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
