/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `org.example:stars` example (seqspec:viewer#extension-points): a namespaced layer type end to end on the renderer
// API alone. Its format half: the type (params table, kernel fallback `points`) and a field provider,
// `org.example:affinity` (pKd from a `kd_nm` column, on annotations carrying `custom["org.example"]`). Its renderer
// reads `custom` and the provider's values in prepare only, as data; draws a star per item sized by its affinity
// through `mark`, colored by its color channel; hits its stars; lists its colors and a size sample; rings the hovered
// star in one call (`ringItem`, the hover emphasis: seqspec:interaction#hover).

import type { FieldProviderEntry, SvsRegistryEntry } from '@molstar/seq-core/format';
import { int, numberColumn, OptionalField, SimpleParamsSchema } from '@molstar/seq-core/format';
import {
  type AxisView,
  type DrawView,
  type Features,
  type Hit,
  hitItem,
  type LayerRenderer,
  legendOf,
  mark,
  NO_COLOR,
  PathBatch,
  PathVerb,
  PreparedParts,
  ringItem,
  SizeUnit,
  ThemeLocation,
} from '../../draw/api.js';

const NS = 'org.example';
export const AffinityField: FieldProviderEntry = {
  name: `${NS}:affinity`,
  type: 'number',
  applies: (a) => NS in a.custom.doc,
  read(a, out, strings) {
    const kd = a.fields.find((f) => f.name === strings.find('kd_nm') && f.kind === 'number');
    for (let r = 0; kd && r < a.count; r++) out[r] = 9 - Math.log10(kd.data[r] as number);
  },
};
export const StarsFormat: SvsRegistryEntry = {
  fields: [AffinityField],
  layerTypes: [
    {
      name: `${NS}:stars`,
      options: SimpleParamsSchema({ star_points: OptionalField(int, 5, 'the points of each star') }),
      accepts: ['ranges'],
      layouts: ['tracks', 'text'],
      axes: ['blocks'],
      sizing: 'free',
      fallback: 'points',
      channels: [{ name: 'color', kind: 'color', unit: null, default: null }],
    },
  ],
};

export interface Stars extends Features {
  /** Per item, its radius (CSS px) from its affinity; the unit `custom["org.example"]` names. */
  readonly radius: Float64Array;
  /** Per star, its affinity (pKd) and its item's id in its annotation (what a hover names). */
  readonly pkd: Float64Array;
  readonly item: Int32Array;
  readonly unit: string;
}
const paths = new PathBatch();
const at = new Float64Array(2);
export const StarsRenderer: LayerRenderer<Stars> = {
  name: `${NS}:stars`,
  // Drawn on the batch API, so the engine maps it onto rings too (seqspec:viewer#the-render-backend).
  geometries: ['linear', 'wrapped', 'circular'],
  stable: false,
  minHeight: 8,
  content: false,
  height: () => 24,
  prepare(input, rc) {
    const f = PreparedParts.points.prepare(input, rc) as Features;
    const radius = new Float64Array(f.count);
    const pkd = new Float64Array(f.count);
    const item = new Int32Array(f.count);
    let unit = '';
    for (let k = 0; k < f.count; k++) {
      const a = input.set.annotations[f.annotation[k]];
      const v = a ? (numberColumn(a, input.set, AffinityField.name, rc.registry.svs)?.[f.row[k]] ?? 6) : 6;
      radius[k] = Math.max(3, Math.min(10, v - 3));
      pkd[k] = v;
      item[k] = a ? a.item[f.row[k]] : -1;
      unit ||= String((a?.custom.doc[NS] as { unit?: string } | undefined)?.unit ?? '');
    }
    return { ...f, radius, pkd, item, unit };
  },
  draw(v: DrawView, s: Stars, p: unknown) {
    const n = Math.max(3, Number((p as { options?: { star_points?: number } }).options?.star_points ?? 5));
    const axis = v.axis;
    paths.reset();
    for (let k = 0; axis && k < s.count; k++) {
      if (!mark(v, s.block[k], s.start[k], v.top + v.height / 2, at)) continue;
      const loc = ThemeLocation.at(v.loc, axis, v.layer, s.block[k], s.start[k]);
      loc.annotation = s.annotation[k];
      loc.row = s.row[k];
      paths.begin(v.color ? v.color.value(loc) : NO_COLOR);
      for (let i = 0; i <= 2 * n; i++) {
        const r = (i % 2 === 0 ? s.radius[k] : s.radius[k] / 2.5) * v.dpr;
        const a = -Math.PI / 2 + (i * Math.PI) / n;
        const o = paths.verb(i === 0 ? PathVerb.move : i === 2 * n ? PathVerb.close : PathVerb.line);
        if (i === 2 * n) continue;
        paths.coords[o] = at[0] + r * Math.cos(a);
        paths.coords[o + 1] = at[1] + r * Math.sin(a);
      }
    }
    v.backend.paths(paths);
  },
  hit(v: AxisView, s: Stars, _p: unknown, out: Hit): boolean {
    const y = v.top + v.height / 2;
    for (let k = 0; v.axis && k < s.count; k++) {
      const sk = v.shownOf(s.block[k], s.start[k]);
      if (s.block[k] !== out.block || sk < 0) continue;
      const x = (v.x(s.block[k], sk) + v.x(s.block[k], (sk + 1) as never)) / 2;
      if (Math.hypot(out.point[0] - x, out.point[1] - y) <= s.radius[k])
        return hitItem(out, v.axis.set, s.annotation[k], s.row[k]);
    }
    return false;
  },
  // The hovered star (`marks.item` of this layer): ringed in its color's hover ink, inside the entry's band.
  overlay(v: DrawView, s: Stars) {
    const item = v.marks.item(v.layer);
    const ann = v.marks.itemAnnotation(v.layer);
    for (let k = 0; item >= 0 && v.axis && k < s.count; k++) {
      if (s.item[k] !== item || s.annotation[k] !== ann || !mark(v, s.block[k], s.start[k], v.top + v.height / 2, at))
        continue;
      const loc = ThemeLocation.at(v.loc, v.axis, v.layer, s.block[k], s.start[k]);
      loc.annotation = s.annotation[k];
      loc.row = s.row[k];
      const half = (v.height * v.dpr) / 2;
      ringItem(v, v.color ? v.color.value(loc) : NO_COLOR, 0, at[0], at[1], s.radius[k] * v.dpr, half);
    }
  },
  // `custom["org.example"]`, read in prepare, labels the star under the pointer (seqspec:viewer#labels).
  label(h, _lc, out, s) {
    for (let k = 0; k < s.count; k++)
      if (s.annotation[k] === h.item.annotation && s.item[k] === h.item.item) {
        out.line('affinity', `${s.pkd[k].toFixed(1)} ${s.unit}`);
        return;
      }
  },
  legend(p, themes, out, s) {
    legendOf(p, themes, out);
    out.size(`affinity (${s.unit})`, 10, SizeUnit.px);
  },
};
