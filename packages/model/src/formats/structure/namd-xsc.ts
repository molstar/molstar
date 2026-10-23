/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PeriodicBox } from '@molstar/core/math/geometry/periodic-box';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { XscEntry, XscFile } from '@molstar/io/reader/xsc/parser';
import { Coordinates, type Frame } from '@molstar/model/model/structure/coordinates';

// A NAMD extended system file's boxes given to coordinates: a DCD keeps a box's lengths and
// angles but not where it is (NAMD's `cellOrigin`, its center), nor its vectors' orientation. Each frame takes the
// entry of its step: the last at or before it (a box is written every `xstFreq` steps, frames every `dcdFreq`).

/** The box of an entry: its vectors, its corner the center less half of a + b + c. */
export function xscPeriodicBox(entry: XscEntry): { box: PeriodicBox; cell: Cell } | undefined {
  const a = Vec3.create(...entry.a);
  const b = Vec3.create(...entry.b);
  const c = Vec3.create(...entry.c);
  const corner = Vec3.create(...entry.origin);
  Vec3.scaleAndAdd(corner, corner, a, -0.5);
  Vec3.scaleAndAdd(corner, corner, b, -0.5);
  Vec3.scaleAndAdd(corner, corner, c, -0.5);
  const box = PeriodicBox.create(a, b, c, corner);
  if (!PeriodicBox.isValid(box)) return;
  return { box, cell: Cell.fromBasis(a, b, c) };
}

/** The entry for a frame: by its step when it has one, else by index when the counts agree, else the last box. */
function entryFor(xsc: XscFile, frame: Frame, index: number, frameCount: number): XscEntry {
  const { entries } = xsc;
  if (frame.step !== undefined) {
    let found = entries[0];
    for (const e of entries) {
      if (e.step > frame.step) break;
      found = e;
    }
    return found;
  }
  if (entries.length === frameCount) return entries[index];
  return entries[entries.length - 1];
}

/** The coordinates with each frame's box (and cell) from the extended system file; the frames' arrays are shared. */
export function coordinatesWithXsc(coordinates: Coordinates, xsc: XscFile): Coordinates {
  const n = coordinates.frames.length;
  const frames = coordinates.frames.map((f, i) => {
    const periodic = xscPeriodicBox(entryFor(xsc, f, i, n));
    // the same `xyzOrdering` object: a frame's arrays are reordered in place once, whichever frame object asks
    return periodic ? { ...f, cell: periodic.cell, box: periodic.box } : f;
  });
  return Coordinates.create(frames, coordinates.deltaTime, coordinates.timeOffset);
}
