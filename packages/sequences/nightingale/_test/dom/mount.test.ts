/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `mount` against Nightingale's own elements (NG-14), in Chromium: what a manager shows after its children's `change`
// events (a hover's highlight, a zoom's display range), a section's `highlights` and `focus` kept as Nightingale
// moves (NG-01), a scroll group synced both ways, a host's highlight, a feature's hover event and extent, and a group's
// toggle.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import {
  compile,
  createRegistries,
  KernelRegistry,
  MvsRegistry,
  register,
  resolveView,
  type SeqViewSpecDoc,
} from "@molstar/seq-core/format";
import { KeyKind, type KeyTarget, type Loci, OrderedSet, type SeqSet, type UnitId } from "@molstar/seq-core/model";
import { afterEach, describe, expect, it } from "vitest";
import { materialize, NightingaleSupported } from "../../src/index.js";
import { type Mounted, mount, type NightingaleEvent } from "../../src/mount.js";

const reg = createRegistries();
register(reg, [KernelRegistry, MvsRegistry]);
const SEQ = "MKVLAAGIVGLLLAQWERTYIPASDFGHKLCVNMQWERTYIPASDFGHKLCVNMKVLAAG";

type Builder = ReturnType<typeof createSeqViewSpec>;
let mounted: Mounted | null = null;
let container: HTMLElement | null = null;
afterEach(() => {
  mounted?.dispose();
  container?.remove();
  [mounted, container] = [null, null];
});

async function open(build: (b: Builder, s: string) => void, onEvent?: (e: NightingaleEvent) => void) {
  const b = createSeqViewSpec({ title: "mount" });
  b.data.sequence({ ref: "s", value: SEQ, label: "S" });
  build(b, "s");
  const doc: SeqViewSpecDoc = b.getDocument();
  const { set } = await compile(doc, reg, { load: () => Promise.reject(new Error("no files")) });
  const { view } = resolveView(doc, set, "v", reg, NightingaleSupported);
  const { root } = materialize(view, set);
  container = document.createElement("div");
  container.style.width = "800px";
  document.body.appendChild(container);
  mounted = mount(root, container, { view, set }, onEvent);
  await mounted.ready;
  const managers = [...container.querySelectorAll("nightingale-manager")];
  /** A manager's drawn elements (Nightingale's, the manager aside). */
  const elements = (m: Element) => [...m.querySelectorAll("*")].filter((x) => x.localName.startsWith("nightingale-"));
  return { set, managers, elements };
}
/** What a child sends its manager: Nightingale's `change`, bubbling. */
const change = (from: Element, detail: Record<string, unknown>) =>
  from.dispatchEvent(new CustomEvent("change", { detail, bubbles: true }));
/** Loci over `s`'s unit, indices [a, b). */
function lociOf(set: SeqSet, a: number, b: number): Loci {
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  set.byName("s", t);
  let unit = t.id;
  if (t.kind === KeyKind.molecule)
    for (let u = 0; u < set.units.count; u++) if (set.units.molecule[u] === t.id) unit = u;
  return {
    kind: "seq-loci",
    set,
    elements: [{ unit: unit as UnitId, indices: OrderedSet.ofRange(a, b) as never }],
  } as Loci;
}
const shown = (els: Element[], attr: string) => [...new Set(els.map((x) => x.getAttribute(attr)))];

const sites = (b: Builder, s: string) =>
  b.data.annotation({
    ref: "sites",
    geometry: "ranges",
    target: s,
    value: { start: [2, 9], stop: [4, 12], id: ["a", "b"] },
  });

describe("a section's highlights (NG-01)", () => {
  it("stay lit under a component's hover and after its leave; a host's highlight joins them", async () => {
    const { set, managers, elements } = await open((b, s) => {
      sites(b, s);
      const highlights = [{ data: "sites", items: ["a", "b"] }];
      const v = b.view({ ref: "v", label: "V" });
      v.section({ ref: "one", highlights } as never)
        .block({ sequence: s })
        .track({ label: "R" })
        .letters();
    });
    const els = elements(managers[0]);
    expect(shown(els, "highlight")).toEqual(["3:4,10:12"]);
    change(els[0], { highlight: "20:20" });
    expect(shown(els, "highlight")).toEqual(["3:4,10:12,20:20"]);
    change(els[0], { highlight: undefined });
    expect(shown(els, "highlight")).toEqual(["3:4,10:12"]);
    mounted?.highlight(lociOf(set, 29, 31));
    expect(shown(els, "highlight")).toEqual(["3:4,10:12,30:31"]);
    mounted?.highlight(null);
    expect(shown(els, "highlight")).toEqual(["3:4,10:12"]);
  });
});

describe("a section's focus (NG-01)", () => {
  it("opens the section; a zoom moves it, and a later hover keeps the zoom", async () => {
    const { managers, elements } = await open((b, s) => {
      const v = b.view({ ref: "v", label: "V" });
      const one = v.section({ ref: "one" });
      one.focus({ target: s, start: 9, stop: 40 });
      one.block({ sequence: s }).track({ label: "R" }).letters();
    });
    const els = elements(managers[0]);
    expect([shown(els, "display-start"), shown(els, "display-end")]).toEqual([["10"], ["40"]]);
    change(els[0], { "display-start": 30, "display-end": 50 });
    expect([shown(els, "display-start"), shown(els, "display-end")]).toEqual([["30"], ["50"]]);
    change(els[0], { highlight: "35:35" });
    expect([shown(els, "display-start"), shown(els, "display-end")]).toEqual([["30"], ["50"]]);
  });
});

describe("a scroll group", () => {
  it("copies a zoom to the group's other sections, both ways, and each can still zoom after", async () => {
    const { managers, elements } = await open((b, s) => {
      const v = b.view({ ref: "v", label: "V" });
      v.section({ ref: "one", scroll_group: "x" } as never)
        .block({ sequence: s })
        .track({ label: "A" })
        .letters();
      v.section({ ref: "two", scroll_group: "x" } as never)
        .block({ sequence: s })
        .track({ label: "B" })
        .letters();
    });
    const [a, b] = [elements(managers[0]), elements(managers[1])];
    change(a[0], { "display-start": 5, "display-end": 25 });
    expect([shown(b, "display-start"), shown(b, "display-end")]).toEqual([["5"], ["25"]]);
    change(b[0], { "display-start": 12, "display-end": 33 });
    expect([shown(a, "display-start"), shown(a, "display-end")]).toEqual([["12"], ["33"]]);
    change(a[0], { "display-start": 2, "display-end": 20 });
    expect([shown(a, "display-start"), shown(b, "display-start")]).toEqual([["2"], ["2"]]);
  });
});

describe("events", () => {
  it("a feature's hover names its block, layer, item and every fragment; leave clears", async () => {
    const got: NightingaleEvent[] = [];
    const { managers } = await open(
      (b, s) => {
        sites(b, s);
        b.view({ ref: "v", label: "V" }).block({ sequence: s }).track({ label: "A" }).ranges({ data: "sites" });
      },
      (e) => got.push(e),
    );
    const track = managers[0].querySelector("nightingale-track, nightingale-track-canvas") as Element;
    const feature = (track as unknown as { data: { svs: unknown; start: number; end: number }[] }).data[1];
    change(track, { eventType: "mouseover", feature, highlight: `${feature.start}:${feature.end}` });
    change(track, { eventType: "mouseout", feature, highlight: undefined });
    expect(got.map((e) => [e.kind, e.block, e.annotation >= 0, e.item, e.extent])).toEqual([
      ["hover", 0, true, 1, [{ block: 0, start: 9, stop: 12 }]],
      ["leave", -1, true, 1, []],
    ]);
  });

  it("a group's header collapses and expands it in place", async () => {
    const { managers } = await open((b, s) => {
      const v = b.view({ ref: "v", label: "V" });
      const sec = v.block({ sequence: s });
      const g = sec.group({ label: "G" } as never);
      g.track({ label: "A" }).letters();
      g.track({ label: "B" }).letters();
    });
    const group = managers[0].querySelector(".svs-group") as Element;
    const button = group.querySelector(".svs-toggle") as HTMLElement;
    button.click();
    expect([group.hasAttribute("data-collapsed"), button.getAttribute("aria-expanded")]).toEqual([true, "false"]);
    button.click();
    expect([group.hasAttribute("data-collapsed"), button.getAttribute("aria-expanded")]).toEqual([false, "true"]);
  });
});
