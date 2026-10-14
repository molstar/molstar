/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's `style` (seqspec:format#sections): a preset's name, or an inline preset whose `preset` names the base and
// whose other keys, a curated list in the format's own names (snake_case, JSON base types), override that preset's
// render props. Resolve maps them onto the viewer's names (seqspec:viewer#render-props) here, so the viewer may rename
// its props without breaking a document. A bad key or value is a warning (`svs.section.style`) naming what is allowed;
// the section is drawn with the rest.

import { Issue } from "@molstar/seq-core/model";
import { type CheckContext, type FieldType, isPlainObject, ptr, show, typeError } from "../schema/field-schema.ts";

/** An inline preset: `preset` the base (absent or null: the section's default for its layout and axis). */
export interface SectionStyle {
  preset?: string | null;
  column_heads?: "shared" | "own_row";
  data_floor?: number;
  row_height?: number;
  header_width?: number;
  vertical_nav?: "scrollbar" | "overview" | "both" | "none";
  wrapped_group_labels?: "above_end" | "inline_start" | "none";
  row_borders?: boolean;
  hover_emphasis?: boolean;
}

/** What an inline style sets, in the viewer's render props' names and values. */
export interface StyleRenderProps {
  dataColumnHeads?: "shared" | "ownRow";
  dataFloor?: number;
  rowHeight?: number;
  headerWidth?: number;
  verticalNav?: "scrollbar" | "overview" | "both" | "none";
  wrappedGroupLabels?: "aboveEnd" | "inlineStart" | "none";
  rowBorders?: boolean;
  hoverEmphasis?: boolean;
}

/** One key: the render prop it sets, and its values (an enumeration, the format's value to the prop's), px or a bool. */
export interface SectionStyleKey {
  readonly prop: keyof StyleRenderProps;
  readonly values?: Readonly<Record<string, string>>;
  /** A number of CSS px, at least this. */
  readonly min?: number;
  readonly bool?: true;
}
const same = (...v: string[]) => Object.fromEntries(v.map((x) => [x, x]));

/** The keys besides `preset`, in the order the format's table lists them. */
export const SectionStyleKeys: Readonly<Record<Exclude<keyof SectionStyle, "preset">, SectionStyleKey>> = {
  column_heads: { prop: "dataColumnHeads", values: { shared: "shared", own_row: "ownRow" } },
  data_floor: { prop: "dataFloor", min: 0 },
  row_height: { prop: "rowHeight", min: 1 },
  header_width: { prop: "headerWidth", min: 0 },
  vertical_nav: { prop: "verticalNav", values: same("scrollbar", "overview", "both", "none") },
  wrapped_group_labels: {
    prop: "wrappedGroupLabels",
    values: { above_end: "aboveEnd", inline_start: "inlineStart", none: "none" },
  },
  row_borders: { prop: "rowBorders", bool: true },
  hover_emphasis: { prop: "hoverEmphasis", bool: true },
};
const KEYS = Object.keys(SectionStyleKeys) as (keyof typeof SectionStyleKeys)[];
const keyOf = (k: string) => (Object.hasOwn(SectionStyleKeys, k) ? SectionStyleKeys[k as (typeof KEYS)[number]] : null);

/** What a key accepts, for messages. */
function allowed(key: SectionStyleKey): string {
  if (key.values)
    return `one of ${Object.keys(key.values)
      .map((v) => `"${v}"`)
      .join(", ")}`;
  return key.bool ? "true or false" : `a number of px, at least ${key.min}`;
}
const valid = (key: SectionStyleKey, v: unknown) =>
  key.values
    ? typeof v === "string" && Object.hasOwn(key.values, v)
    : key.bool
      ? typeof v === "boolean"
      : typeof v === "number" && Number.isFinite(v) && v >= (key.min ?? 0);

function warn(ctx: CheckContext, path: string, message: string) {
  ctx.issues.push(Issue.warning("svs.section.style", path, `${message}; it is ignored, the rest of the style applies`));
}

/** `style`'s type: a preset name, or an inline preset; its bad keys and values are warnings. */
export const sectionStyle: FieldType<string | SectionStyle> = {
  name: "str | SectionStyle",
  check(v, path, ctx) {
    if (typeof v === "string") return;
    if (!isPlainObject(v)) return typeError(ctx, path, "a preset name or an inline preset {preset, ...keys}", v);
    for (const [k, x] of Object.entries(v)) {
      const at = ptr(path, k);
      if (k === "preset") {
        if (x !== null && typeof x !== "string") warn(ctx, at, `preset ${show(x)} is not a preset name or null`);
        continue;
      }
      const key = keyOf(k);
      if (!key) {
        const message = `unknown style key "${k}"; the keys: preset, ${KEYS.join(", ")}`;
        if (ctx.lenient)
          ctx.issues.push(Issue.warning("svs.version.newer_minor", at, `${message} (a newer minor version)`));
        else warn(ctx, at, message);
      } else if (!valid(key, x)) warn(ctx, at, `${k} ${show(x)} is not ${allowed(key)}`);
    }
  },
};

/** A section's resolved style: the preset named (null: the default) and the valid keys as render props (null: none). */
export function resolveStyle(v: unknown): { preset: string | null; render: StyleRenderProps | null } {
  if (typeof v === "string") return { preset: v, render: null };
  if (!isPlainObject(v)) return { preset: null, render: null };
  let render: Record<string, unknown> | null = null;
  for (const k of KEYS) {
    const key = SectionStyleKeys[k];
    const x = v[k];
    if (x === undefined || !valid(key, x)) continue;
    render ??= {};
    render[key.prop] = key.values ? key.values[x as string] : x;
  }
  return { preset: typeof v.preset === "string" ? v.preset : null, render: render as StyleRenderProps | null };
}
