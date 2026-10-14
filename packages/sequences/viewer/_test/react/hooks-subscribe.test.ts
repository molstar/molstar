/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The React hooks' store subscription (`useSeqState`, `useDerived`): a state's replay on subscribe is not a change, so
// it calls nothing back; each change calls back once.

import { describe, expect, it } from "vitest";
import { subscriber } from "../../src/react/hooks.js";
import { live, state } from "../../src/state.js";

describe("the hooks' store subscription", () => {
  it("skips each state's replay on subscribe and calls back once per change", () => {
    const a = state(1);
    const l = live({ x: 0 });
    let calls = 0;
    const off = subscriber([a, l.version])(() => calls++);
    expect(calls).toBe(0);
    a.set(2);
    l.bump();
    expect(calls).toBe(2);
    off();
    a.set(3);
    expect(calls).toBe(2);
  });
});
