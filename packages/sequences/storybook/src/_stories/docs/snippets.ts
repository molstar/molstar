/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The code snippets of the Seq* docs pages, kept here so their indentation survives (MDX trims it in attributes).

export const BUILDER = `import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";

const b = createSeqViewSpec({ title: "A protein following its coding DNA" });
const gene = b.data.sequence({ ref: "gene", alphabet: "dna", value: GENE });
const protein = translate(b.data, { ref: "protein", from: gene, cds: [[12, 46], [76, 105]] });
const v = b.view({ ref: "pair" });
const dna = v.section({ ref: "pair.dna" }).block({ sequence: gene });
dna.track({ label: "DNA" }).letters();
const p = v.section({ ref: "pair.protein", follow: dna.handle }).block({ sequence: protein });
p.track({ label: "Protein" }).letters();
const doc = b.getDocument();`;

export const STARS = `// org.example:stars, abridged (@molstar/seq-viewer/testing)
export const StarsFormat: SvsRegistryEntry = {
  fields: [AffinityField],
  layerTypes: [{
    name: "org.example:stars",
    options: SimpleParamsSchema({ star_points: OptionalField(int, 5, "the points of each star") }),
    accepts: ["ranges"], layouts: ["tracks", "text"], axes: ["blocks"], sizing: "free",
    fallback: "points",
    channels: [{ name: "color", kind: "color", unit: null, default: null }],
  }],
};
export const StarsRenderer: LayerRenderer<Stars> = {
  name: "org.example:stars", geometries: ["linear", "wrapped", "circular"], stable: false, minHeight: 8, content: false,
  height: () => 24,
  prepare(input, rc) { /* points' features, plus a radius per item from the affinity field */ },
  draw(v, s, p) { /* a star path per item at mark(v, block, index, y, out), colored by v.color */ },
  hit(v, s, p, out) { /* the star under out.point: hitItem(out, set, annotation, row) */ },
  overlay(v, s) {
    // The hovered star (v.marks.item(v.layer), v.marks.itemAnnotation(v.layer)), ringed in one call: the ink derived
    // from its color for the skin, inside its cell; nothing when the section turns hover emphasis off.
    ringItem(v, color, GlyphShape.circle, x, y, radius, half);
  },
  label(h, lc, out, s) { out.line("affinity", …); },
  legend(p, themes, out, s) { legendOf(p, themes, out); out.size("affinity", 10, SizeUnit.px); },
};
// registered as { svs: StarsFormat, renderers: [StarsRenderer] }`;

export const CONTRACT = `// @molstar/seq-core/format: the port a host implements
interface Loaders {
  load(request: LoadRequest, signal: AbortSignal): Promise<Uint8Array>;
  stream?(source: StreamSource): StreamProvider | null;   // by \`stream\` id or \`protocol\`; null: svs.provider.unavailable
}
interface StreamSource {
  ref: string; stream: string | null; protocol: string | null; uri: string | null;
  data_kind: "letters" | "items" | "values" | "alignment";
  parent: string | null;                      // a nested source's parent; null: the document's
  declared: Record<string, unknown>;          // its readers' row_count, column_count, length, target; a statistic's measure
}
interface StreamProvider {
  // Every request, the manifest first. Rejects with StreamNotServed, any other error (a failure), or signal.reason.
  fetch(request: StreamRequest, signal: AbortSignal): Promise<DataPatch>;
}
class StreamNotServed extends Error {}        // reported once, never retried`;
