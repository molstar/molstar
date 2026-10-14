/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Nested sources are asked per depth, siblings together: the genes' manifest keys 5 sources, each of whose manifests
// keys 2 more. With every answer after a latency, the 5 are asked at once and the 10 once those answered, so the whole
// tree is in after two latencies, not after one per source (the Storybook story "Templates and nested sources"). The
// clock is fake (`vi.useFakeTimers`): the latencies pass when the test advances them, and the times are exact.

import type { DataPatch, Loaders, Manifest, StreamNode } from "@molstar/seq-core/format";
import { afterEach, describe, expect, it, vi } from "vitest";
import { opened } from "./layout-fixtures.ts";
import { nested, P, range, SECTION, templateDoc } from "./template-fixtures.ts";

const LATENCY = 60;

describe("nested sources, asked by depth", () => {
  afterEach(() => {
    vi.useRealTimers();
  });

  it("the siblings of a depth are asked together, each depth after the one above answered", async () => {
    vi.useFakeTimers({ toFake: ["setTimeout", "clearTimeout", "performance"] });
    const manifests: Record<string, StreamNode[]> = { genes: [], notes: [] };
    for (let g = 1; g <= 5; g++) {
      manifests.genes.push(range(`G${g}`, g), nested(`G${g}/tx`, `tx-G${g}`));
      manifests[`tx-G${g}`] = [1, 2].flatMap((k) => {
        manifests[`dom-G${g}.${k}`] = [range("d1", g + k)];
        return [range(`G${g}.${k}`, g + k), nested(`G${g}.${k}/dom`, `dom-G${g}.${k}`)];
      });
    }
    const t0 = performance.now();
    const asked: [string, number][] = [];
    const loaders: Loaders = {
      load: (r) => Promise.reject(new Error(`no file ${r.uri}`)),
      stream: (s) =>
        s.stream && s.stream in manifests
          ? {
              fetch(r) {
                asked.push([s.stream as string, performance.now() - t0]);
                const m: Manifest = {
                  systems: [{ key: "p", extent: [0, P.length], numbering: null, label: null }],
                  levels: [],
                };
                const at = { source: r.source, system: null, start: 0, stop: 0, row_start: 0, row_stop: 0, level: -1 };
                const patch = { ...at, kind: "manifest", manifest: { ...m, nodes: manifests[s.stream as string] } };
                const wait = s.stream === "genes" || s.stream === "notes" ? 0 : LATENCY;
                return new Promise((ok) => setTimeout(ok, wait)).then(() => patch as DataPatch);
              },
            }
          : null,
    };
    const doc = templateDoc({ for_each: { kind: "annotation", source: "genes", depth: null } });
    const opening = opened(doc, { loaders, size: { [SECTION]: [1000, 600] } });
    await vi.runAllTimersAsync();
    const t = await opening;
    const total = performance.now() - t0;
    const at = (prefix: string) => asked.filter(([id]) => id.startsWith(prefix)).map(([, ms]) => ms);
    const spread = (xs: number[]) => Math.max(...xs) - Math.min(...xs);
    expect([at("tx-").length, at("dom-").length]).toEqual([5, 10]);
    // Together: a serial walk would put them LATENCY apart each (4 × LATENCY across the five).
    expect(spread(at("tx-"))).toBeLessThan(LATENCY);
    expect(spread(at("dom-"))).toBeLessThan(LATENCY);
    expect(Math.min(...at("dom-"))).toBeGreaterThanOrEqual(Math.max(...at("tx-")) + LATENCY * 0.8);
    expect(total).toBeLessThan(LATENCY * 3.5);
    t.ctx.dispose();
  });
});
