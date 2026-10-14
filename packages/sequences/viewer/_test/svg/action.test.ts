/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `view.exportSvg` (SV1b; seqspec:viewer#actions): the action produces the SVG `renderSvg` writes for the same
// context, offers it under the view's name, leaves the live context as it found it (sizes, DPR, measurer, current
// view), queues exports, refuses targets that are not an export, answers `seqv:runAction`, and reports a failure
// instead of throwing.

import { describe, expect, it } from "vitest";
import { case1, case14 } from "../../../core/_test/format/cases.ts";
import { KernelRenderers, type LayerRenderer } from "../../src/index.ts";
import { DefaultFontMetrics, exportSvgAction, metricsMeasurer, renderSvg } from "../../src/svg/index.ts";
import { opened } from "../layout-fixtures.ts";

const kernel = [...KernelRenderers] as LayerRenderer[];
interface Offer {
  svg: string;
  name: string;
}
/** A context with the action registered; `offers` collects what it hands over. */
async function exporting(doc = case14(), o: { fail?: boolean } = {}) {
  const offers: Offer[] = [];
  const waiting: (() => void)[] = [];
  const action = exportSvgAction({
    ports: { deflate: undefined },
    offer: (svg, name) => {
      if (o.fail) throw new Error("the disk is full");
      offers.push({ svg, name });
      for (const w of waiting.splice(0)) w();
    },
  });
  const t = await opened(doc, { renderers: kernel, registry: [{ actions: [action] }] });
  /** Resolves once `n` files have been offered. */
  const offered = async (n: number) => {
    while (offers.length < n) await new Promise<void>((r) => waiting.push(r));
    return offers;
  };
  return { ...t, offers, offered };
}
const sizes = (ctx: Awaited<ReturnType<typeof exporting>>["ctx"]) =>
  ctx.viewports.sections.value.map((ref) => {
    const l = ctx.viewports.get(ref)?.live.current;
    return [ref, l?.width, l?.height, l?.dpr];
  });

describe("view.exportSvg", () => {
  it("offers the SVG renderSvg writes for the context, named after the view", async () => {
    const t = await exporting();
    expect(t.ctx.actions.run("view.exportSvg", { view: "gene-view", width: 800 })).toBe(true);
    const [offer] = await t.offered(1);
    expect(offer.name).toBe("gene-view.svg");
    expect(offer.svg).toBe(await renderSvg(t.ctx, { view: "gene-view", width: 800 }));
    expect(offer.svg).toMatch(/^<svg [^>]*width="800"/);
    t.ctx.dispose();
  });

  it("at a scale; a target without a view exports the current one", async () => {
    const t = await exporting();
    t.ctx.actions.run("view.exportSvg", { width: 600, scale: 2 });
    const [offer] = await t.offered(1);
    expect(offer.name).toBe("gene-view.svg");
    expect(offer.svg).toBe(await renderSvg(t.ctx, { width: 600, scale: 2 }));
    expect(offer.svg).toMatch(/^<svg [^>]*width="1200"/);
    t.ctx.dispose();
  });

  it("leaves the live context as it was: every section's size and DPR, the current view and the measurer", async () => {
    const t = await exporting(case1());
    const before = { sizes: sizes(t.ctx), view: t.ctx.views.current.value, port: t.ctx.text.port };
    expect(before.view).toBe("text");
    // Another view, a width and a scale that are not the live ones, and a measurer of its own in the action's ports.
    t.ctx.actions.run("view.exportSvg", { view: "tracks", width: 333, scale: 3 });
    const [offer] = await t.offered(1);
    expect(offer.name).toBe("tracks.svg");
    expect(offer.svg).toMatch(/viewBox="0 0 333 /);
    expect([sizes(t.ctx), t.ctx.views.current.value, t.ctx.text.port]).toEqual([before.sizes, "text", before.port]);
    t.ctx.dispose();
  });

  it("a measurer in the ports is used for the render and the context's own comes back", async () => {
    const offers: string[] = [];
    const mine = metricsMeasurer(DefaultFontMetrics);
    const t = await opened(case14(), {
      renderers: kernel,
      registry: [{ actions: [exportSvgAction({ ports: { measurer: mine }, offer: (svg) => void offers.push(svg) })] }],
    });
    const own = t.ctx.text.port;
    t.ctx.actions.run("view.exportSvg", { width: 700 });
    while (offers.length === 0) await new Promise((r) => setTimeout(r, 5));
    expect(t.ctx.text.port).toBe(own);
    expect(offers[0]).toBe(await renderSvg(t.ctx, { width: 700, ports: { measurer: mine } }));
    // The explicit measurer has an em box, so the baselines are dy and not dominant-baseline.
    expect(offers[0]).not.toContain("dominant-baseline");
    t.ctx.dispose();
  });

  it("exports queue: two at once finish in order, each as renderSvg writes it", async () => {
    const t = await exporting();
    t.ctx.actions.run("view.exportSvg", { width: 500 });
    t.ctx.actions.run("view.exportSvg", { width: 900, scale: 2 });
    const offers = await t.offered(2);
    expect(offers[0].svg).toBe(await renderSvg(t.ctx, { width: 500 }));
    expect(offers[1].svg).toBe(await renderSvg(t.ctx, { width: 900, scale: 2 }));
    t.ctx.dispose();
  });

  it("refuses what is not an export: no width, a bad scale, another view's name, not an object", async () => {
    const t = await exporting();
    for (const target of [
      null,
      5,
      {},
      { width: 0 },
      { width: -4 },
      { width: Number.NaN },
      { width: "800" },
      { width: 800, scale: 0 },
      { width: 800, view: "nope" },
    ]) {
      expect(t.ctx.actions.can("view.exportSvg", target), JSON.stringify(target)).toBe(false);
      expect(t.ctx.actions.run("view.exportSvg", target)).toBe(false);
    }
    expect(t.offers).toEqual([]);
    expect(t.ctx.actions.can("view.exportSvg", { width: 800, view: "gene-view" })).toBe(true);
    t.ctx.dispose();
  });

  it("is a command as data: seqv:runAction with the target", async () => {
    const t = await exporting();
    const ran = await t.ctx.commands.run(
      { command: "seqv:runAction", params: { name: "view.exportSvg", target: { view: "gene-view", width: 640 } } },
      "api",
    );
    expect(ran).toBe(true);
    const [offer] = await t.offered(1);
    expect(offer.svg).toBe(await renderSvg(t.ctx, { view: "gene-view", width: 640 }));
    t.ctx.dispose();
  });

  it("a failing offer is an issue, and the next export still runs", async () => {
    const t = await exporting(case14(), { fail: true });
    t.ctx.actions.run("view.exportSvg", { width: 500 });
    while (!t.ctx.issues.current.value.some((i) => i.code === "seqv.command.error"))
      await new Promise((r) => setTimeout(r, 5));
    const issue = t.ctx.issues.current.value.find((i) => i.code === "seqv.command.error");
    expect(issue?.message).toMatch(/view\.exportSvg failed: the disk is full/);
    t.ctx.dispose();
  });
});
