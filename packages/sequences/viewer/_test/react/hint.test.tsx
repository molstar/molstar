/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The React UI's default `Hint` never blocks the pointer (seqspec:interaction#accessibility): the cases are
// dom/hint-cases.ts, shared with the DOM UI's.

import { act, useLayoutEffect, useRef } from "react";
import { createRoot } from "react-dom/client";
import { describe } from "vitest";
import { Hint } from "../../src/react/index.ts";
import { hintCases } from "../dom/hint-cases.ts";

(globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

/** The test's own button element as the trigger, in a span of the hint's. */
function Adopt({ el }: { el: HTMLElement }) {
  const at = useRef<HTMLSpanElement>(null);
  useLayoutEffect(() => {
    at.current?.append(el);
  }, [el]);
  return <span ref={at} />;
}

describe("the React Hint", () => {
  hintCases(async (host, buttons) => {
    const root = createRoot(host);
    await act(async () =>
      root.render(
        buttons.map((b, i) => (
          <Hint key={b.textContent} label={`Hint ${i}`}>
            <Adopt el={b} />
          </Hint>
        )),
      ),
    );
    return () => act(async () => root.unmount());
  });
});
