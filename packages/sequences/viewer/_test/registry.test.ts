/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer owns `seqv` (seqspec:viewer#the-viewer-spec; seqspec:format#extensions): `DefaultRegistry.svs` holds its
// issue codes and `seqv:map`, and an extension still cannot register under `seqv`.

import { SimpleParamsSchema } from "@molstar/seq-core/format";
import { describe, expect, it, vi } from "vitest";
import { createSeqViewer, SeqRegistrationError, ViewerIssueCodes, ViewerLayerTypes } from "../src/index.js";

describe("seqv", () => {
  it("DefaultRegistry registers the viewer's issue codes and seqv:map", () => {
    const ctx = createSeqViewer();
    expect(ctx.registry.svs.get("layerTypes", "seqv:map")).toBe(ViewerLayerTypes.layerTypes?.[0]);
    expect(ctx.registry.svs.get("issueCodes", "seqv.patch.skipped")?.severity).toBe("warning");
    expect(ctx.registry.svs.list("issueCodes").filter((c) => c.code.startsWith("seqv.")).length).toBe(
      ViewerIssueCodes.issueCodes?.length,
    );
  });

  it("an extension under seqv is refused, nothing registered", () => {
    const ctx = createSeqViewer();
    const version = ctx.registry.svs.version;
    const ext = {
      svs: {
        issueCodes: [{ code: "seqv.mine", severity: "error" as const }],
        custom: [{ namespace: "seqv", kind: "x", schema: SimpleParamsSchema({}) }],
      },
    };
    expect(() => ctx.register(ext)).toThrow(SeqRegistrationError);
    expect(() => ctx.register(ext)).toThrow('issueCodes: "seqv.mine" is reserved');
    expect([ctx.registry.svs.version, ctx.registry.svs.get("issueCodes", "seqv.mine")]).toEqual([version, undefined]);
  });

  it("a spec from another copy of the viewer in the page (an HMR update, a host bundling it twice) is the viewer's own", async () => {
    vi.resetModules();
    const a = await import("../src/index.js");
    vi.resetModules();
    const b = await import("../src/index.js");
    expect(a.ViewerIssueCodes).not.toBe(b.ViewerIssueCodes);
    const ctx = b.createSeqViewer(a.DefaultSeqViewerSpec());
    expect(ctx.registry.svs.get("layerTypes", "seqv:map")?.name).toBe("seqv:map");
    ctx.dispose();
    // A copy of its entries is not: the brand does not spread.
    const other = createSeqViewer();
    const copy = { ...ViewerIssueCodes, issueCodes: [{ code: "seqv.copied", severity: "error" as const }] };
    expect(() => other.register({ svs: copy })).toThrow(SeqRegistrationError);
    other.dispose();
  });
});
