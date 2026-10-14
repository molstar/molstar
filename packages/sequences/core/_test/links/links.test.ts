/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Invocation } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { createLinkRunner, DefaultLinkKeys, formatLink, type LinkKey, parseLink } from "../../src/links/index.js";

/** A recording `run`; the highlight channels it leaves lit are `lit()`. */
function recorder(fail?: "reject" | "throw") {
  const calls: Invocation[] = [];
  const run = (inv: Invocation) => {
    calls.push(inv);
    if (fail === "throw") throw new Error("boom");
    return fail === "reject" ? Promise.reject(new Error("boom")) : undefined;
  };
  const lit = () => {
    const on = new Set<string>();
    for (const c of calls) {
      if (c.command !== "seqv:highlight") continue;
      const p = c.params as { channel: string; target: unknown };
      if (p.target) on.add(p.channel);
      else on.delete(p.channel);
    }
    return [...on];
  };
  return { calls, run, lit };
}

describe("codec", () => {
  it("round trips keys, values and args through URI escapes", () => {
    const href = formatLink("seq-focus", "a b,c&d", { start: 3, stop: 9, section: "x=y/ü" });
    expect(href.startsWith("!seq-focus=a%20b,c%26d&")).toBe(true);
    const [p] = parseLink(href);
    expect(p.key).toBe("seq-focus");
    expect(p.value).toBe("a b,c&d");
    expect(p.args).toEqual({ "seq-focus": "a b,c&d", start: "3", stop: "9", section: "x=y/ü" });
  });
  it("gives one entry per key, Mol*'s included, and none for anything else", () => {
    const l = parseLink("!seq-select=r1&focus-refs=a&seq-highlight=r2", "L");
    expect(l.map((x) => [x.key, x.value, x.link])).toEqual([
      ["seq-select", "r1", "L"],
      ["focus-refs", "a", "L"],
      ["seq-highlight", "r2", "L"],
    ]);
    expect(parseLink("https://example.com/?a=b")).toEqual([]);
    expect(parseLink("#anchor")).toEqual([]);
  });
  it("skips a malformed escape", () => {
    expect(parseLink("!seq-view=%E0%A4%A&seq-select=r").map((x) => x.key)).toEqual(["seq-select"]);
  });
});

describe("default keys", () => {
  it("enter highlights with args converted from strings", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    const href = formatLink("seq-highlight", "u1,u2", {
      start: 3,
      stop: "9",
      items: "i1,i2",
      beg: "10",
      numbering: "pdb",
    });
    expect(runner.run(href, "enter", "L1")).toBe(true);
    expect(r.calls).toEqual([
      {
        command: "seqv:highlight",
        params: {
          channel: "links:L1",
          closure: "containing",
          target: { refs: ["u1", "u2"], args: { start: 3, stop: 9, beg: "10", numbering: "pdb", items: ["i1", "i2"] } },
        },
      },
    ]);
  });
  it("leaves absent and non-integer args out", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    runner.run(formatLink("seq-highlight", "u1", { start: "1.5", stop: "", items: "" }), "enter", "L");
    expect((r.calls[0].params as { target: unknown }).target).toEqual({ refs: ["u1"], args: {} });
  });
  it("seq-select clicks as a set, with the closure", () => {
    const r = recorder();
    createLinkRunner(r.run).run(formatLink("seq-select", "u1"), "click", "L");
    expect(r.calls).toEqual([
      { command: "seqv:select", params: { target: { refs: ["u1"], args: {} }, mode: "set", closure: "containing" } },
    ]);
  });
  it("seq-focus clicks with its section, and only it", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    runner.run(formatLink("seq-focus", "u1", { section: "main" }), "click", "L");
    runner.run(formatLink("seq-focus", "u1"), "click", "L");
    runner.run(formatLink("seq-select", "u1", { section: "main" }), "click", "L");
    expect(r.calls.map((c) => (c.params as { section?: string }).section)).toEqual(["main", undefined, undefined]);
  });
  it("snapshot and view clicks take their ref", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    runner.run("!seq-snapshot=next", "click", "L");
    runner.run("!seq-snapshot=", "click", "L");
    runner.run("!seq-view=v2", "click", "L");
    expect(runner.run("!seq-view=", "click", "L")).toBe(true);
    expect(r.calls.map((c) => [c.command, c.params])).toEqual([
      ["seqv:setSnapshot", { ref: "next" }],
      ["seqv:setSnapshot", { ref: null }],
      ["seqv:setView", { ref: "v2" }],
    ]);
  });
  it("runs nothing for a ref-less target, and an event a key has no entry for", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    expect(runner.run("!seq-select=", "enter", "L")).toBe(true);
    expect(runner.run("!seq-snapshot=next", "enter", "L")).toBe(true);
    expect(r.calls).toEqual([]);
  });
  it("ignores other viewers' keys and other links", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    expect(runner.run("!focus-refs=a", "click", "L")).toBe(false);
    expect(runner.run("!nope=a&focus-refs=b", "enter", "L")).toBe(false);
    expect(runner.run("https://example.com", "click", "L")).toBe(false);
    expect(r.calls).toEqual([]);
  });
  it("runs a key beside Mol*'s in one link", () => {
    const r = recorder();
    expect(createLinkRunner(r.run).run("!focus-refs=a&seq-view=v1", "click", "L")).toBe(true);
    expect(r.calls).toHaveLength(1);
  });
});

describe("the runner", () => {
  const href = formatLink("seq-focus", "u1");
  it("leave clears only its link's channel", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    runner.run(href, "enter", "A");
    runner.run(href, "enter", "B");
    runner.run(href, "leave", "A");
    expect(r.lit()).toEqual(["links:B"]);
    runner.release("A");
    expect(r.calls).toHaveLength(3);
  });
  it("release runs the pending leave of a link removed while hovered, once", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    runner.run(href, "enter", "A");
    runner.release("A");
    runner.release("A");
    expect(r.lit()).toEqual([]);
    expect(r.calls).toHaveLength(2);
  });
  it("dispose leaves every entered link", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run);
    for (const l of ["A", "B", "C"]) runner.run(href, "enter", l);
    runner.run(href, "leave", "B");
    expect(r.lit()).toEqual(["links:A", "links:C"]);
    runner.dispose();
    expect(r.lit()).toEqual([]);
    runner.dispose();
    expect(r.calls).toHaveLength(3 + 1 + 2);
  });
  it("catches a rejected or throwing run, and keeps running", async () => {
    for (const fail of ["reject", "throw"] as const) {
      const r = recorder(fail);
      const unhandled: unknown[] = [];
      const on = (e: unknown) => unhandled.push(e);
      process.on("unhandledRejection", on);
      const runner = createLinkRunner(r.run);
      expect(runner.run(href, "enter", "A")).toBe(true);
      expect(() => runner.run(href, "click", "A")).not.toThrow();
      await new Promise((res) => setTimeout(res, 5));
      process.off("unhandledRejection", on);
      expect(unhandled).toEqual([]);
      expect(r.calls).toHaveLength(2);
    }
  });
  it("takes an application's keys beside the defaults", () => {
    const r = recorder();
    const pulse: LinkKey = {
      key: "lab:pulse",
      value: "refs",
      click: (l) => ({ command: "lab:pulse", params: { refs: l.value.split(",") } }),
    };
    const runner = createLinkRunner(r.run, { keys: [...DefaultLinkKeys, pulse] });
    expect(runner.run("!lab:pulse=a,b", "click", "L")).toBe(true);
    expect(r.calls).toEqual([{ command: "lab:pulse", params: { refs: ["a", "b"] } }]);
  });
  it("keeps a bounded cache of parsed hrefs", () => {
    const r = recorder();
    const runner = createLinkRunner(r.run, { cacheLinks: 2 });
    for (const v of ["a", "b", "c", "a"]) runner.run(formatLink("seq-select", v), "click", "L");
    expect(r.calls.map((c) => (c.params as { target: { refs: string[] } }).target.refs[0])).toEqual([
      "a",
      "b",
      "c",
      "a",
    ]);
  });
});
