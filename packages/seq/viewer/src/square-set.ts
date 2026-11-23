/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SquareGripAxis, SquareNavigator } from '@molstar/seqx-core/render/navigator-square';

// The square over a set (_spec/README.md "Composition", Navigators, The square over a set): one section on both axes
// with `navigator.switch: include` (a ring's square), its window the blocks the section shows. The run's box is dragged
// along the diagonal (a move) or by an edge or corner (a resize), snapped to whole blocks, drawn as a preview while it
// runs and written once on release (one patch, one undo step); a click on a cell shows its pair, on a label its block,
// with Cmd or Ctrl adding or removing them; the keys move, resize and pick cells. What it writes is the bound section's
// `include` (every block: none), through `write`.

/** What a set's change writes: the lattice's names in its order, or undefined to show every block. */
export type SquareSetWrite = (include: readonly string[] | undefined) => boolean;

/** The cursor over the set's parts: a move, a resize along x, y or the diagonal. */
export type SquareSetCursor = 'move' | 'ew-resize' | 'ns-resize' | 'nwse-resize';

/** A grip's cursor: the body moves, an edge resizes along its axis, a corner along the diagonal. */
export function setCursor(grip: 'body' | 'start' | 'end', axis: SquareGripAxis): SquareSetCursor {
  if (grip === 'body') return 'move';
  return axis === 'x' ? 'ew-resize' : axis === 'y' ? 'ns-resize' : 'nwse-resize';
}

/**
 * The include for a set of lattice blocks (`shown(k)`): their names in the lattice's order, undefined for every block;
 * an empty set is none (`null`: nothing to write, a ring of nothing is refused).
 */
export function includeOf(square: SquareNavigator, shown: (k: number) => boolean): string[] | undefined | null {
  const names: string[] = [];
  const n = square.blockCount;
  for (let k = 0; k < n; k++) if (shown(k)) names.push(square.blockName(k));
  if (names.length === 0) return null;
  return names.length === n ? undefined : names;
}

/** A readout's parts (the navigator's `ReadoutParts`). */
interface Readout {
  readonly name: string;
  readonly value: string;
  readonly rest: string;
}

/** What a cell says: its two blocks (`chr4 × chr7`; one on the diagonal), and what the picture adds. */
export function setCellReadout(square: SquareNavigator, i: number, j: number): Readout | undefined {
  if (i < 0 || j < 0 || i >= square.blockCount || j >= square.blockCount) return undefined;
  const pair = i === j ? square.blockName(i) : `${square.blockName(i)} × ${square.blockName(j)}`;
  const more = square.cellText(i, j);
  return { name: '', value: pair, rest: more ? ` · ${more}` : '' };
}

/** What a point on the square over a set says: a cell's pair, a label's block, else nothing. */
export function setReadout(square: SquareNavigator, x: number, y: number): Readout | undefined {
  const label = square.labelAt(x, y);
  if (label >= 0) return { name: '', value: square.blockName(label), rest: '' };
  const cell = square.cellAt(x, y);
  return cell ? setCellReadout(square, cell.x, cell.y) : undefined;
}

/** The set's gesture on the run's box, and the keyboard's cell. Its drag allocates nothing. */
export class SquareSetGesture {
  /** The grip pressed (`none`: no drag of the box), its axis, where, and the run then. */
  private grip: 'none' | 'body' | 'start' | 'end' = 'none';
  private axis: SquareGripAxis = 'both';
  private b0 = 0;
  private i0 = -1;
  private i1 = -1;
  /** The run the drag shows now. */
  private n0 = -1;
  private n1 = -1;
  /** The keyboard's cell (lattice indices, -1: none yet): what Enter shows. */
  keyX = -1;
  keyY = -1;
  /** The square (its set) the keyboard's cell was placed on. */
  private keyOf: SquareNavigator | undefined;

  constructor(private readonly write: SquareSetWrite) {}

  /** Whether a drag of the box runs. */
  get dragging(): boolean {
    return this.grip !== 'none';
  }

  /** The grip pressed, for the cursor while it runs. */
  get cursor(): SquareSetCursor {
    return this.grip === 'none' ? 'move' : setCursor(this.grip, this.axis);
  }

  /** A press: on the run's box a drag begins (true); elsewhere none. */
  press(square: SquareNavigator, x: number, y: number, touch = false): boolean {
    const g = square.gripAt(x, y, touch);
    const r = square.run;
    if (!g || r.i0 < 0) {
      this.grip = 'none';
      return false;
    }
    this.grip = g.grip;
    this.axis = g.axis;
    this.i0 = this.n0 = r.i0;
    this.i1 = this.n1 = r.i1;
    this.b0 = square.blockAtPoint(x, y, 'both');
    return true;
  }

  /** The pointer moved in a drag: the run it would write, snapped to whole blocks, drawn as the preview. */
  drag(square: SquareNavigator): void {
    // the point is the square's `pointer` (set by the caller: nothing boxed a move)
    if (this.grip === 'none') return;
    const n = square.blockCount;
    const len = this.i1 - this.i0;
    if (this.grip === 'body') {
      const shift = square.blockAtPointer('both') - this.b0;
      this.n0 = Math.max(0, Math.min(n - len, this.i0 + shift));
      this.n1 = this.n0 + len;
    } else if (this.grip === 'start') {
      this.n0 = Math.max(0, Math.min(this.i1 - 1, square.edgeAtPointer(this.axis)));
      this.n1 = this.i1;
    } else {
      this.n0 = this.i0;
      this.n1 = Math.max(this.i0 + 1, Math.min(n, square.edgeAtPointer(this.axis)));
    }
    square.setPreview(this.n0, this.n1);
  }

  /** The drag ends: its run written when it changed (one patch); true when it was. */
  release(square: SquareNavigator): boolean {
    const was = this.grip;
    this.grip = 'none';
    square.setPreview(-1, -1);
    if (was === 'none' || (this.n0 === this.i0 && this.n1 === this.i1)) return false;
    return this.writeRun(square, this.n0, this.n1);
  }

  /** The drag is cancelled: nothing written, the shown run drawn again. */
  cancel(square: SquareNavigator): void {
    this.grip = 'none';
    square.setPreview(-1, -1);
  }

  /** Writes a run `[i0, i1)` of lattice blocks (every block: no include). */
  writeRun(square: SquareNavigator, i0: number, i1: number): boolean {
    const include = includeOf(square, (k) => k >= i0 && k < i1);
    return include !== null && this.write(include);
  }

  /** A click on a cell: its pair alone; with `toggle` (Cmd or Ctrl) its two blocks added, or removed when both shown. */
  cell(square: SquareNavigator, i: number, j: number, toggle: boolean): boolean {
    if (!toggle) return this.writeChange(square, (k) => k === i || k === j);
    const both = square.shows(i) && square.shows(j);
    return this.writeChange(square, (k) => (k === i || k === j ? !both : square.shows(k)));
  }

  /** A click on a label: its block alone; with `toggle` added or removed. */
  label(square: SquareNavigator, k: number, toggle: boolean): boolean {
    if (!toggle) return this.writeChange(square, (b) => b === k);
    return this.writeChange(square, (b) => (b === k ? !square.shows(b) : square.shows(b)));
  }

  /** Writes a set when it differs from the shown one (none shown: refused). */
  private writeChange(square: SquareNavigator, next: (k: number) => boolean): boolean {
    let same = true;
    for (let k = 0; k < square.blockCount && same; k++) if (next(k) !== square.shows(k)) same = false;
    if (same) return false;
    const include = includeOf(square, next);
    return include !== null && this.write(include);
  }

  /**
   * A key on the focused square (true when taken): arrows move the run a block (Left and Up back), Shift+arrows grow or
   * shrink its end; Cmd or Ctrl+arrows move the keyboard's cell, as plain arrows do over a scattered set; Enter shows
   * the cell's pair, with Cmd or Ctrl adds or removes it.
   */
  key(square: SquareNavigator, key: string, shift: boolean, toggle: boolean): boolean {
    const n = square.blockCount;
    const r = square.run;
    const back = key === 'ArrowLeft' || key === 'ArrowUp';
    const arrow = back || key === 'ArrowRight' || key === 'ArrowDown';
    this.ensureKeyCell(square);
    if (key === 'Enter') return this.cell(square, this.keyX, this.keyY, toggle);
    if (!arrow || n === 0) return false;
    if (toggle || r.i0 < 0) {
      if (shift) return false;
      const step = back ? -1 : 1;
      if (key === 'ArrowLeft' || key === 'ArrowRight') this.keyX = Math.max(0, Math.min(n - 1, this.keyX + step));
      else this.keyY = Math.max(0, Math.min(n - 1, this.keyY + step));
      return true;
    }
    if (shift) {
      const i1 = Math.max(r.i0 + 1, Math.min(n, r.i1 + (back ? -1 : 1)));
      return i1 !== r.i1 && this.writeRun(square, r.i0, i1);
    }
    const len = r.i1 - r.i0;
    const i0 = Math.max(0, Math.min(n - len, r.i0 + (back ? -1 : 1)));
    if (i0 === r.i0) return false;
    // the keyboard's cell moves with the run
    this.keyX = Math.max(0, Math.min(n - 1, this.keyX + i0 - r.i0));
    this.keyY = Math.max(0, Math.min(n - 1, this.keyY + i0 - r.i0));
    return this.writeRun(square, i0, i0 + len);
  }

  /**
   * The keyboard's cell, first the run's first diagonal cell (or the first shown block's); picked again when the set
   * shown changed under it (`square` a new set's: a cell of blocks no longer shown is none of the window's).
   */
  ensureKeyCell(square: SquareNavigator): void {
    const n = square.blockCount;
    if (square !== this.keyOf) {
      this.keyOf = square;
      if (!square.shows(this.keyX) || !square.shows(this.keyY)) this.keyX = this.keyY = -1;
    }
    if (this.keyX >= 0 && this.keyX < n && this.keyY >= 0 && this.keyY < n) return;
    let k = Math.max(0, square.run.i0);
    if (square.run.i0 < 0) for (k = 0; k < n && !square.shows(k); k++);
    this.keyX = this.keyY = Math.min(k, Math.max(0, n - 1));
  }
}
