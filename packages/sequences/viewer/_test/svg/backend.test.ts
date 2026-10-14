/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The SVG backend (SV1a; seqspec:viewer#the-render-backend) in Node: every batch kind written and read back by a small
// XML parser, grouping by (color, alpha), the numbers' format, clips, passes, images and layers, the PNG, escaping, the
// DPR's division, and the entry's freedom from the DOM and Node.

import { readdirSync, readFileSync } from "node:fs";
import { deflateSync, inflateSync } from "node:zlib";
import { describe, expect, it } from "vitest";
import { itWithin } from "../../../core/_test/model/timing.ts";
import {
  GlyphBatch,
  GlyphShape,
  NO_COLOR,
  Pass,
  PathBatch,
  PathMode,
  PathVerb,
  RectBatch,
  SectorBatch,
  SHAPE_COUNT,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "../../src/draw/api.ts";
import { DefaultSkin } from "../../src/index.ts";
import { encodePng, SvgBackend, type SvgOptions } from "../../src/svg/index.ts";
import { all, type El, parseXml, rectsOf, tokens } from "./xml.ts";

const GLYPHS = [
  { char: "W", path: null, box: [1, 1] as const, fallback: 0 },
  { char: null, path: 'M0 0 L10 0 L0 10 Z "x"', box: [10, 20] as const, fallback: 0 },
];
const TOKENS: Record<number, string> = {
  [SkinToken.ink]: DefaultSkin.ink,
  [SkinToken.selection]: DefaultSkin.selection,
  [SkinToken.surface]: DefaultSkin.surface,
};
const resources = {
  font: (id: number) => (id === 1 ? '600 16px "Inter", serif' : "400 12px ui-sans-serif, sans-serif"),
  fontPx: (id: number) => (id === 1 ? 16 : 12),
  token: (t: number) => TOKENS[t] ?? "red",
  glyph: (id: number) => GLYPHS[id - SHAPE_COUNT] ?? null,
};
/** A backend with one main pass of w × h device px at `dpr`, drawn by `f`; the document parsed. */
async function doc(f: (be: SvgBackend) => void, o: SvgOptions = {}, w = 200, h = 100, dpr = 1) {
  const be = new SvgBackend(resources, o);
  be.beginPass(Pass.main, w, h, dpr);
  f(be);
  be.endPass();
  const xml = await be.toString();
  return { be, xml, svg: parseXml(xml) };
}
const paths = (svg: El) => all(svg, "path").filter((p) => p.attrs.d !== undefined);

describe("rectangles", () => {
  it("one path per (color, alpha) in order of first appearance, each rectangle a subpath", async () => {
    const { svg, be } = await doc((be) => {
      const b = new RectBatch();
      b.begin(0xff0000);
      b.push(1, 2, 3, 4);
      b.begin(0x0000ff, 128);
      b.push(10, 20, 5, 5);
      b.begin(0xff0000);
      b.push(30, 2, 7, 4);
      b.begin(NO_COLOR);
      b.push(0, 0, 9, 9);
      b.begin(0xff0000, 128);
      b.push(50, 50, 1, 1);
      be.rects(b);
    });
    const p = paths(svg);
    expect(p.map((e) => [e.attrs.fill, e.attrs["fill-opacity"]])).toEqual([
      ["#ff0000", undefined],
      ["#0000ff", "0.502"],
      ["#ff0000", "0.502"],
    ]);
    expect(p.every((e) => e.attrs["shape-rendering"] === "crispEdges")).toBe(true);
    expect(rectsOf(p[0].attrs.d)).toEqual([
      [1, 2, 3, 4],
      [30, 2, 7, 4],
    ]);
    expect(rectsOf(p[1].attrs.d)).toEqual([[10, 20, 5, 5]]);
    expect(be.stats.rects).toBe(4);
  });

  it("an outline is four edges inside the rectangle, grouped with fills of its color; a thick one is a fill", async () => {
    const { svg } = await doc((be) => {
      const b = new RectBatch();
      b.begin(0x00ff00, 255, 2);
      b.push(10, 10, 20, 10);
      b.push(50, 10, 4, 10);
      be.rects(b);
    });
    expect(rectsOf(paths(svg)[0].attrs.d)).toEqual([
      [10, 10, 20, 2],
      [10, 18, 20, 2],
      [10, 12, 2, 6],
      [28, 12, 2, 6],
      [50, 10, 4, 10],
    ]);
  });

  it("skin tokens: the CSS color read to #rrggbb, its own opacity times the alpha byte", async () => {
    const { svg } = await doc((be) => {
      const b = new RectBatch();
      b.begin(SkinToken.selection, 255);
      b.push(0, 0, 5, 5);
      b.begin(SkinToken.ink, 255);
      b.push(0, 0, 5, 5);
      be.rects(b);
    });
    const [sel, ink] = paths(svg);
    expect([sel.attrs.fill, sel.attrs["fill-opacity"], ink.attrs.fill]).toEqual(["#2b7fff", "0.1", "#1b1f24"]);
    expect(ink.attrs["fill-opacity"]).toBeUndefined();
  });

  it("device px over the DPR: CSS px in the viewBox, the file at the DPR; neighbors share edges exactly", async () => {
    const { svg } = await doc(
      (be) => {
        const b = new RectBatch();
        b.begin(0x123456);
        // Thirds at DPR 3: each edge rounded on its own, the width the difference of two edges.
        for (let i = 0; i < 6; i++) b.push(i, 0, 1, 3);
        be.rects(b);
      },
      {},
      30,
      9,
      3,
    );
    expect([svg.attrs.width, svg.attrs.height, svg.attrs.viewBox]).toEqual(["30", "9", "0 0 10 3"]);
    const r = rectsOf(paths(svg)[0].attrs.d);
    expect(r.map((x) => [x[0], Math.round((x[0] + x[2]) * 100) / 100])).toEqual([
      [0, 0.33],
      [0.33, 0.67],
      [0.67, 1],
      [1, 1.33],
      [1.33, 1.67],
      [1.67, 2],
    ]);
    expect(r.every((x) => x[3] === 1)).toBe(true);
  });
});

describe("size and time", () => {
  itWithin(
    "the writer is linear: 200,000 rectangles in two colors, a few seconds at most, a constant size each",
    async (within) => {
      const sizes: number[] = [];
      for (const n of [25_000, 200_000]) {
        const { xml } = await within(n === 25_000 ? 1000 : 4000, () =>
          doc((be) => {
            const b = new RectBatch(n);
            for (let i = 0; i < n; i++) {
              b.begin(i & 1 ? 0xff0000 : 0x0000ff, 200);
              b.push(i % 1000, i >> 10, 2, 3);
            }
            be.rects(b);
          }),
        );
        sizes.push(xml.length / n);
      }
      expect(Math.abs(sizes[1] - sizes[0]) / sizes[0]).toBeLessThan(0.15);
      expect(sizes[1]).toBeLessThan(30);
    },
  );
});

describe("sectors", () => {
  it("an arc of a ring: outer clockwise, a line, the inner back; a full ring as two halves per circle", async () => {
    const { svg, be } = await doc((be) => {
      const s = new SectorBatch();
      s.begin(0x884400);
      s.push(50, 50, 10, 30, 0, Math.PI / 2);
      s.push(50, 50, 10, 30, 0, (3 * Math.PI) / 2);
      s.push(50, 50, 10, 30, 0, 2 * Math.PI);
      s.push(50, 50, 10, 30, 1, 1);
      be.sectors(s);
    });
    const d = paths(svg)[0].attrs.d.split("M").filter(Boolean);
    expect(d).toHaveLength(3);
    // A quarter: from (80, 50) to (50, 80), then to the inner radius at (50, 60) and back to (60, 50).
    expect(tokens(`M${d[0]}`)).toEqual([
      "M",
      80,
      50,
      "A",
      30,
      30,
      0,
      0,
      1,
      50,
      80,
      "L",
      50,
      60,
      "A",
      10,
      10,
      0,
      0,
      0,
      60,
      50,
      "Z",
    ]);
    // Three quarters: the large-arc flag.
    expect(tokens(d[1])[6]).toBe(1);
    // A full ring: both circles in two arcs each.
    expect(d[2].match(/A/g)).toHaveLength(4);
    expect(be.stats.sectors).toBe(4);
  });
});

describe("paths", () => {
  it("consecutive subpaths of one style are one path; fills and strokes; every verb", async () => {
    const { svg } = await doc((be) => {
      const p = new PathBatch();
      p.begin(0x00aa00, 255, PathMode.fill);
      p.moveTo(1, 2);
      p.lineTo(3, 4);
      p.close();
      p.begin(0x00aa00, 255, PathMode.fill);
      p.moveTo(5, 6);
      p.quadTo(7, 8, 9, 10);
      p.begin(0x0000aa, 128, PathMode.stroke, 3);
      p.moveTo(0, 0);
      let o = p.verb(PathVerb.cubic);
      p.coords.set([1, 1, 2, 2, 3, 0], o);
      o = p.verb(PathVerb.arc);
      p.coords.set([10, 10, 5, 0, Math.PI, 0], o);
      p.begin(NO_COLOR);
      p.moveTo(0, 0);
      p.begin(0xaa0000);
      o = p.verb(PathVerb.roundRect);
      p.coords.set([10, 10, 20, 10, 2], o);
      be.paths(p);
    });
    const [a, b, c] = paths(svg);
    expect([a.attrs.fill, tokens(a.attrs.d)]).toEqual([
      "#00aa00",
      ["M", 1, 2, "L", 3, 4, "Z", "M", 5, 6, "Q", 7, 8, 9, 10],
    ]);
    expect([b.attrs.stroke, b.attrs["stroke-opacity"], b.attrs.fill, b.attrs["stroke-width"]]).toEqual([
      "#0000aa",
      "0.502",
      "none",
      "1.5",
    ]);
    // The cubic, then the arc: a line to its start (15, 10), then a half circle clockwise to (5, 10).
    expect(tokens(b.attrs.d)).toEqual(["M", 0, 0, "C", 1, 1, 2, 2, 3, 0, "L", 15, 10, "A", 5, 5, 0, 0, 1, 5, 10]);
    expect(c.attrs.d.startsWith("M12 10H28A2 2 0 0 1 30 12V18")).toBe(true);
    expect(paths(svg)).toHaveLength(3);
  });
});

describe("glyphs", () => {
  it("shapes join their color's path, characters are text 2r high, a path glyph is a def and a use per mark", async () => {
    const { svg, be } = await doc((be) => {
      const g = new GlyphBatch();
      g.font = 0;
      g.begin(0x222222);
      for (let k = 0; k <= GlyphShape.tick; k++) g.push(10 + 20 * k, 20, 6, k);
      g.push(30, 60, 8, SHAPE_COUNT);
      g.push(60, 60, 5, SHAPE_COUNT + 1);
      g.push(90, 60, 5, SHAPE_COUNT + 1);
      be.glyphs(g);
    });
    const shapes = paths(svg).filter((p) => p.attrs.id === undefined);
    expect(shapes).toHaveLength(1);
    expect(shapes[0].attrs.d.match(/M/g)?.length).toBe(9);
    const text = all(svg, "text");
    expect([text[0].text, text[0].attrs.x, text[0].attrs.y, text[0].attrs["font-size"]]).toEqual([
      "W",
      "30",
      "60",
      "16",
    ]);
    expect(text[0].attrs.x).toBe("30");
    // The path glyph: one def (its data escaped), two uses scaled so the larger side (20) is 2r.
    const def = paths(svg).filter((p) => p.attrs.id === "p9");
    expect(def).toHaveLength(1);
    expect(def[0].attrs.d).toBe('M0 0 L10 0 L0 10 Z "x"');
    const uses = all(svg, "use");
    expect(uses.map((u) => u.attrs["xlink:href"])).toEqual(["#p9", "#p9"]);
    expect(uses[0].attrs.transform).toBe("translate(57.5 55) scale(0.5)");
    expect(be.stats.glyphs).toBe(11);
  });
});

describe("text", () => {
  it("font split into family, size and weight; anchor and baseline; an unrotated run is x and y", async () => {
    const { svg } = await doc((be) => {
      const t = new TextBatch();
      t.style(1, TextAlign.right, TextBaseline.middle);
      t.begin(0x112233);
      t.push(100.4, 50, "a < b & c");
      t.begin(0x112233, 128);
      t.push(5, 6, "  two  spaces ");
      be.text(t);
      const u = new TextBatch();
      u.style(0, TextAlign.left, TextBaseline.alphabetic);
      u.begin(SkinToken.ink);
      u.push(1, 2, "Z");
      be.text(u);
    });
    const [g, g2, g3] = all(svg, "g").filter((e) => e.attrs["font-family"] !== undefined);
    expect([
      g.attrs["font-family"],
      g.attrs["font-size"],
      g.attrs["font-weight"],
      g.attrs["text-anchor"],
      g.attrs["dominant-baseline"],
    ]).toEqual(['"Inter", serif', "16", "600", "end", "central"]);
    expect(g.attrs["xml:space"]).toBe("preserve");
    expect([g.children[0].text, g.children[0].attrs.x, g.children[0].attrs.y]).toEqual(["a < b & c", "100.4", "50"]);
    expect([g2.attrs["fill-opacity"], g2.children[0].text]).toEqual(["0.502", "  two  spaces "]);
    // The defaults write nothing: left, alphabetic, regular.
    expect([g3.attrs["text-anchor"], g3.attrs["dominant-baseline"], g3.attrs["font-weight"], g3.attrs.fill]).toEqual([
      undefined,
      undefined,
      undefined,
      "#1b1f24",
    ]);
    expect(g3.attrs["font-family"]).toBe("ui-sans-serif, sans-serif");
  });

  it("a rotated or scaled item is translate(x y) rotate(a) scale(sx sy), the scale about the anchor before the angle", async () => {
    const { svg } = await doc((be) => {
      const t = new TextBatch();
      t.style(0, TextAlign.center, TextBaseline.alphabetic, -Math.PI / 2);
      t.begin(0);
      t.push(10, 20, "up");
      const s = new TextBatch();
      s.style(0, TextAlign.center, TextBaseline.alphabetic);
      s.begin(0);
      let k = s.add("A");
      s.x[k] = 30;
      s.y[k] = 40;
      s.sx[k] = 1.5;
      s.sy[k] = 0.75;
      k = s.add("B");
      s.x[k] = 50;
      s.y[k] = 40;
      be.text(t);
      be.text(s);
    });
    const t = all(svg, "text");
    expect(t.map((e) => e.attrs.transform ?? `${e.attrs.x} ${e.attrs.y}`)).toEqual([
      "translate(10 20) rotate(-90)",
      "translate(30 40) scale(1.5 0.75)",
      "50 40",
    ]);
  });

  it("characters XML cannot carry are dropped, markup escaped", async () => {
    const { xml, svg } = await doc((be) => {
      const t = new TextBatch();
      t.style(0, TextAlign.left, TextBaseline.alphabetic);
      t.begin(0);
      t.push(0, 0, 'a\u0000b\u0007<i>&"cé\u{1f9ec}');
      be.text(t);
    });
    expect(all(svg, "text")[0].text).toBe('ab<i>&"cé\u{1f9ec}');
    expect(xml).toContain("ab&lt;i&gt;&amp;&quot;c");
  });
});

describe("clips and passes", () => {
  it("clips are clipPaths defined once per rectangle, nested groups, closed at the end of the pass", async () => {
    const { svg, be } = await doc((be) => {
      const r = new RectBatch();
      r.begin(0);
      r.push(0, 0, 1, 1);
      for (let k = 0; k < 3; k++) {
        be.pushClip(10, 10, 50, 50);
        be.pushClip(20, 20, 10, 10);
        be.rects(r);
        be.popClip();
        be.popClip();
      }
      be.pushClip(1, 2, 3, 4);
      be.rects(r);
    });
    expect(all(svg, "clipPath").map((c) => Object.values(c.children[0].attrs).join(" "))).toEqual([
      "10 10 50 50",
      "20 20 10 10",
      "1 2 3 4",
    ]);
    const groups = all(svg, "g").filter((g) => g.attrs["clip-path"]);
    expect(groups).toHaveLength(7);
    expect(groups[0].attrs["clip-path"]).toBe("url(#c0)");
    expect(be.stats.clips).toBe(7);
  });

  it("the overlay and the navigator are groups only on request; their calls are still seen", async () => {
    const f = (be: SvgBackend) => {
      const r = new RectBatch();
      r.begin(0);
      r.push(0, 0, 1, 1);
      be.rects(r);
      for (const pass of [Pass.overlay, Pass.nav, Pass.navOverlay]) {
        be.endPass();
        be.beginPass(pass, 200, 100, 1);
        be.pushClip(0, 0, 5, 5);
        be.rects(r);
      }
    };
    const none = await doc(f, { trace: true });
    expect(all(none.svg, "g").map((g) => g.attrs.class)).toEqual(["seqv-main"]);
    expect(none.be.trace?.filter((t) => t === "rects 1")).toHaveLength(4);
    expect(none.be.stats.rects).toBe(1);
    const some = await doc(f, { overlay: true });
    expect(
      all(some.svg, "g")
        .map((g) => g.attrs.class)
        .filter(Boolean),
    ).toEqual(["seqv-main", "seqv-overlay"]);
    const every = await doc(f, { overlay: true, nav: true });
    expect(
      all(every.svg, "g")
        .map((g) => g.attrs.class)
        .filter(Boolean),
    ).toEqual(["seqv-main", "seqv-overlay", "seqv-nav", "seqv-nav-overlay"]);
  });

  it("a clear of the whole pass drops what was drawn; copy is an error (no Copy cap)", async () => {
    const { svg, be } = await doc((be) => {
      const r = new RectBatch();
      r.begin(0);
      r.push(0, 0, 1, 1);
      be.rects(r);
      be.clear(10, 10, 5, 5);
      be.rects(r);
      be.clear(0, 0, 200, 100);
      r.push(2, 2, 1, 1);
      be.rects(r);
    });
    expect(paths(svg)).toHaveLength(1);
    expect(rectsOf(paths(svg)[0].attrs.d)).toHaveLength(2);
    expect(() => be.copy(1, 0, 0, 0, 5, 5)).toThrow("no Copy");
  });

  it("section offsets, the document's size, a background, a title and embedded fonts", async () => {
    const be = new SvgBackend(resources, { background: true, title: "A & B", scale: 2 });
    for (const y of [0, 110]) {
      be.beginSection(0, y);
      be.beginPass(Pass.main, 400, 200, 2);
      be.endPass();
      be.endSection();
    }
    const svg = parseXml(
      await be.toString({
        fonts: [{ family: "Inter", weight: 600, style: "normal", data: new Uint8Array([1, 2, 3, 4]), format: "woff2" }],
      }),
    );
    expect([svg.attrs.width, svg.attrs.height, svg.attrs.viewBox]).toEqual(["400", "420", "0 0 200 210"]);
    expect(all(svg, "title")[0].text).toBe("A & B");
    expect(all(svg, "style")[0].text).toBe(
      '@font-face{font-family:"Inter";font-weight:600;font-style:normal;src:url(data:font/woff2;base64,AQIDBA==) format("woff2")}',
    );
    expect(
      all(svg, "g")
        .filter((g) => g.attrs.transform)
        .map((g) => g.attrs.transform),
    ).toEqual(["translate(0 0)", "translate(0 110)"]);
    expect(all(svg, "rect").find((r) => r.attrs.fill)?.attrs).toMatchObject({
      width: "200",
      height: "210",
      fill: "#ffffff",
    });
  });
});

describe("numbers", () => {
  it("shortest decimals at two places; a number that is not finite is 0 and counted", async () => {
    const { xml, be } = await doc((be) => {
      const t = new TextBatch();
      t.style(0, TextAlign.left, TextBaseline.alphabetic);
      t.begin(0);
      for (const [x, y] of [
        [1 / 3, 2 / 3],
        [0.005, -0.004],
        [123456.789, 1e-9],
        [Number.NaN, Number.POSITIVE_INFINITY],
      ])
        t.push(x, y, "n");
      be.text(t);
    });
    for (const m of ['x="0.33" y="0.67"', 'x="0.01" y="0"', 'x="123456.79" y="0"', 'x="0" y="0"'])
      expect(xml).toContain(m);
    expect(xml).not.toMatch(/NaN|Infinity|e-|undefined/);
    expect(be.stats.nonFinite).toBe(2);
  });
});

describe("images and layers", () => {
  /** The pixels of a PNG data URI: [w, h, rgba bytes], its chunks' CRC-32s checked. */
  function decode(uri: string) {
    const png = Uint8Array.from(Buffer.from(uri.replace("data:image/png;base64,", ""), "base64"));
    const v = new DataView(png.buffer);
    expect([...png.subarray(0, 8)]).toEqual([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]);
    const idat: Uint8Array[] = [];
    let w = 0;
    let h = 0;
    for (let at = 8; at < png.length; ) {
      const n = v.getUint32(at);
      const type = String.fromCharCode(...png.subarray(at + 4, at + 8));
      const body = png.subarray(at + 8, at + 8 + n);
      if (type === "IHDR") [w, h] = [v.getUint32(at + 8), v.getUint32(at + 12)];
      if (type === "IDAT") idat.push(body);
      expect(v.getUint32(at + 8 + n), `${type}'s CRC`).toBe(crc(png.subarray(at + 4, at + 8 + n)));
      at += 12 + n;
    }
    const raw = inflateSync(Buffer.concat(idat));
    const out = new Uint8Array(w * h * 4);
    for (let y = 0; y < h; y++) {
      expect(raw[y * (1 + 4 * w)]).toBe(0);
      out.set(raw.subarray(y * (1 + 4 * w) + 1, (y + 1) * (1 + 4 * w)), y * 4 * w);
    }
    return { w, h, rgba: out };
  }
  const crc = (b: Uint8Array) => {
    let c = 0xffffffff;
    for (const x of b) {
      c ^= x;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    }
    return (c ^ 0xffffffff) >>> 0;
  };

  it("pixels become one PNG data URI per image drawn, a <use> each; a source rectangle is cut by a clip", async () => {
    const px = new Uint32Array(3 * 2);
    new Uint8Array(px.buffer).set([
      1, 2, 3, 255, 4, 5, 6, 128, 7, 8, 9, 0, 10, 11, 12, 255, 13, 14, 15, 255, 16, 17, 18, 255,
    ]);
    const { svg, be } = await doc((be) => {
      const id = be.images.fromPixels(3, 2, px);
      const unused = be.images.fromPixels(3, 2, px);
      be.images.dispose(unused);
      be.image(id, Float64Array.of(0, 0, 3, 2, 10, 10, 30, 20), false);
      be.image(id, Float64Array.of(1, 0, 1, 2, 50, 0, 10, 20), true);
      be.images.dispose(id);
    });
    const imgs = all(svg, "image");
    expect(imgs).toHaveLength(1);
    expect([imgs[0].attrs.id, imgs[0].attrs.width, imgs[0].attrs.height]).toEqual(["i0", "3", "2"]);
    const { w, h, rgba } = decode(imgs[0].attrs["xlink:href"]);
    expect([w, h, [...rgba.subarray(0, 8)]]).toEqual([3, 2, [1, 2, 3, 255, 4, 5, 6, 128]]);
    expect([...rgba.subarray(12)]).toEqual([10, 11, 12, 255, 13, 14, 15, 255, 16, 17, 18, 255]);
    const uses = all(svg, "use");
    expect(uses.map((u) => [u.attrs.transform, u.attrs["image-rendering"]])).toEqual([
      ["translate(10 10) scale(10 10)", "pixelated"],
      ["translate(40 0) scale(10 10)", undefined],
    ]);
    // Only the second cuts: its destination rectangle.
    expect(all(svg, "clipPath")).toHaveLength(1);
    expect(be.stats.images).toBe(2);
  });

  it("a layer stays vector: a group in defs, clipped to its bounds, used scaled; DPR divides its coordinates", async () => {
    const { svg } = await doc(
      (be) => {
        const layer = be.images.beginLayer(40, 20);
        const r = new RectBatch();
        r.begin(0xff00ff);
        r.push(4, 2, 8, 4);
        layer.rects(r);
        const id = be.images.endLayer();
        be.image(id, Float64Array.of(0, 0, 40, 20, 100, 20, 80, 40), true);
      },
      {},
      200,
      100,
      2,
    );
    const g = all(svg, "g").find((e) => e.attrs.id === "i0");
    expect(g?.attrs["clip-path"]).toBeUndefined();
    expect(rectsOf(all(g as El, "path")[0].attrs.d)).toEqual([[2, 1, 4, 2]]);
    // The layer's CSS px (20 × 10) drawn at 2x over (50, 10): scale 2.
    expect(all(svg, "use")[0].attrs.transform).toBe("translate(50 10) scale(2 2)");
    expect(all(svg, "image")).toHaveLength(0);
  });

  it("a stored PNG decodes, and so does one through a deflate port (sync or async)", async () => {
    const px = Uint32Array.from({ length: 300 * 300 }, (_, i) => (i * 2654435761) >>> 0);
    for (const deflate of [
      undefined,
      (b: Uint8Array) => new Uint8Array(deflateSync(b)),
      async (b: Uint8Array) => new Uint8Array(deflateSync(b)),
    ]) {
      const png = await encodePng(300, 300, px, deflate);
      const { w, h, rgba } = decode(`data:image/png;base64,${Buffer.from(png).toString("base64")}`);
      expect([w, h]).toEqual([300, 300]);
      // Compared as bytes: a deep match of 90,000 pixels was most of the test.
      expect(Buffer.from(rgba.buffer, rgba.byteOffset, rgba.byteLength).equals(Buffer.from(px.buffer))).toBe(true);
    }
    // Stored blocks are at most 65,535 bytes: this one needs several, and stays valid.
    expect((await encodePng(300, 300, px)).length).toBeGreaterThan(300 * 300 * 4);
  });
});

describe("platform-free", () => {
  const dir = new URL("../../src/svg/", import.meta.url);
  it("the entry imports nothing of the DOM or Node, and names none of their globals", () => {
    const files = readdirSync(dir).filter((f) => f.endsWith(".ts"));
    expect(files).toContain("backend.ts");
    for (const f of files) {
      const code = readFileSync(new URL(f, dir), "utf8")
        .replace(/\/\*[\s\S]*?\*\//g, "")
        .replace(/\/\/.*$/gm, "");
      expect(code, f).not.toMatch(
        /from "node:|require\(|\b(?:document|window|navigator|OffscreenCanvas|HTMLCanvasElement|ImageData|Buffer|process|btoa|atob|Blob|Image)\b/,
      );
      expect(code, f).not.toMatch(/from "\.\.\/(?:canvas|dom|react)\//);
    }
  });
  it("runs here, where there is no DOM", () => {
    expect(typeof (globalThis as { document?: unknown }).document).toBe("undefined");
  });
});
