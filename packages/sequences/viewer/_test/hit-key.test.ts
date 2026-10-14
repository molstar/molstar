/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover element's key (seqspec:interaction#input-events): `HIT_ELEMENT` and the straight-line `HitKey` cannot drift
// apart.

import { describe, expect, it } from "vitest";
import { HIT_ELEMENT, HitKey } from "../src/hit.js";
import { Hit } from "../src/index.js";

describe("HitKey", () => {
  it("writes the named fields in order, then the system's kind, id and row", () => {
    const h = Hit.create();
    const k = new Int32Array(HitKey.SIZE);
    HIT_ELEMENT.forEach((name, i) => {
      h[name] = 100 + i;
    });
    [h.system.kind, h.system.id, h.system.row] = [200, 201, 202];
    HitKey.write(h, k);
    expect(Array.from(k)).toEqual([...HIT_ELEMENT.map((_, i) => 100 + i), 200, 201, 202]);
    expect(HitKey.SIZE).toBe(HIT_ELEMENT.length + 3);
    expect(HitKey.same(h, k)).toBe(true);
  });

  it("changing any one field of the list, or any slot of the system, is a new element", () => {
    const k = new Int32Array(HitKey.SIZE);
    const h = Hit.create();
    HitKey.write(h, k);
    expect(HitKey.same(h, k)).toBe(true);
    for (const name of HIT_ELEMENT) {
      const was = h[name];
      h[name] = was + 1;
      expect([name, HitKey.same(h, k)]).toEqual([name, false]);
      h[name] = was;
      expect([name, HitKey.same(h, k)]).toEqual([name, true]);
    }
    for (const slot of ["kind", "id", "row"] as const) {
      h.system[slot] += 1;
      expect([slot, HitKey.same(h, k)]).toEqual([slot, false]);
      h.system[slot] -= 1;
    }
    // Fields outside the key (the point, the model index, a location) are not an element's identity.
    h.x += 5;
    h.index += 1;
    h.seq.index = 7 as never;
    expect(HitKey.same(h, k)).toBe(true);
  });

  it("-2 matches no hit: after a leave or a cancel, the next move is a new element", () => {
    const k = new Int32Array(HitKey.SIZE).fill(-2);
    expect(HitKey.same(Hit.create(), k)).toBe(false);
  });
});
