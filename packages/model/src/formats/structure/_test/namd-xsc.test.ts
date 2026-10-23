/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { parseDcd } from '@molstar/io/reader/dcd/parser';
import { parseXsc } from '@molstar/io/reader/xsc/parser';
import { coordinatesFromDcd } from '@molstar/model/formats/structure/dcd';
import { coordinatesWithXsc } from '@molstar/model/formats/structure/namd-xsc';
import { Coordinates, type Frame, Time } from '@molstar/model/model/structure/coordinates';

// NAMD extended system files (.xsc, .xst) give a DCD's frames their boxes with origins.

// A trajectory file of two steps: a 60 Å cube centered at 0, then a 60.5 Å one centered at (0.1, 0, 0)
const xst = `# NAMD extended system trajectory file
#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w
0 60 0 0 0 60 0 0 0 60 0 0 0 0 0 0 0 0 0
500 60.5 0 0 0 60.5 0 0 0 60.5 0.1 0 0 0 0 0 0 0 0
`;

async function xsc(text: string) {
  const r = await parseXsc(text).run();
  if (r.isError) throw new Error(r.message);
  return r.result;
}

const frames = (steps: (number | undefined)[]): Frame[] =>
  steps.map((step, i) => ({
    elementCount: 1,
    time: Time(i, 'step'),
    x: [0],
    y: [0],
    z: [0],
    xyzOrdering: { isIdentity: true },
    ...(step !== undefined ? { step } : {}),
  }));

describe('NAMD extended system files', () => {
  it('reads the steps, box vectors and centers by their labels', async () => {
    const file = await xsc(xst);
    expect(file.entries).toHaveLength(2);
    expect(file.entries[1]).toEqual({
      step: 500,
      a: [60.5, 0, 0],
      b: [0, 60.5, 0],
      c: [0, 0, 60.5],
      origin: [0.1, 0, 0],
    });
    // without labels: NAMD's order, and an origin at 0 when not given
    const plain = await xsc('0 10 0 0 0 10 0 0 0 10\n');
    expect(plain.entries[0].origin).toEqual([0, 0, 0]);
    expect((await parseXsc('# nothing\n').run()).isError).toBe(true);
  });

  it("gives each frame its step's box, the corner half a box from the center", async () => {
    const file = await xsc(xst);
    const original = frames([0, 250, 500, 750]);
    const c = coordinatesWithXsc(Coordinates.create(original, Time(1, 'step'), Time(0, 'step')), file);
    // step 250 takes step 0's box; 750 the last one
    expect(Array.from(c.frames[1].box!.origin)).toEqual([-30, -30, -30]);
    expect(c.frames[2].box!.basis[0]).toBeCloseTo(60.5);
    expect(c.frames[3].box!.origin[0]).toBeCloseTo(0.1 - 30.25);
    expect(c.frames[3].cell!.size[0]).toBeCloseTo(60.5);
    // the frames' arrays and their ordering state are shared: reordered in place once, whichever frame asks
    expect(c.frames[0].x).toBe(original[0].x);
    expect(c.frames[0].xyzOrdering).toBe(original[0].xyzOrdering);
  });

  it('matches frames without steps by index when the counts agree, else takes the last box', async () => {
    const file = await xsc(xst);
    const two = coordinatesWithXsc(
      Coordinates.create(frames([undefined, undefined]), Time(1, 'step'), Time(0, 'step')),
      file,
    );
    expect(two.frames[0].box!.basis[0]).toBeCloseTo(60);
    expect(two.frames[1].box!.basis[0]).toBeCloseTo(60.5);
    const three = coordinatesWithXsc(
      Coordinates.create(frames([undefined, undefined, undefined]), Time(1, 'step'), Time(0, 'step')),
      file,
    );
    for (const f of three.frames) expect(f.box!.basis[0]).toBeCloseTo(60.5);
  });

  it("gives a DCD's frames their steps", async () => {
    const data = fs.readFileSync(
      path.resolve(import.meta.dirname, '../../../../../../data/examples/trajectory/protein.dcd'),
    );
    const parsed = await parseDcd(new Uint8Array(data)).run();
    if (parsed.isError) throw new Error(parsed.message);
    const { ISTART, NSAVC } = parsed.result.header;
    const c = await coordinatesFromDcd(parsed.result).run();
    expect(NSAVC).toBeGreaterThan(0);
    expect(c.frames.map((f) => f.step)).toEqual(c.frames.map((_, i) => ISTART + i * NSAVC));
  });
});
