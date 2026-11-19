/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type DrawLists,
  NoColor,
  PathCommandPoints,
  PathList,
  Passes,
  type PassVisitor,
  walkPass,
} from './draw-lists.js';
import { fmt } from './record.js';
import { pngDataUrl } from './png.js';
import { squirclePath } from './squircle.js';
import { DefaultFonts, type FontStacks } from './tokens.js';

// The SVG backend (_spec/README.md "Render", Backends): draw lists as an SVG string, no DOM. Passes are groups in
// order; a clip is a group with a clip path; squircles are the core's curve as a path; an image is an `<image>` of its
// pixels as a PNG data URL; a text past its max width fades through a gradient mask.
// Every string from data or a spec is escaped (security: labels from data never become markup).

const XmlEscapes: Readonly<Record<string, string>> = {
  '&': '&amp;',
  '<': '&lt;',
  '>': '&gt;',
  '"': '&quot;',
  "'": '&#39;',
};

/**
 * The characters XML 1.0 forbids: controls but tab, newline and return; U+FFFE and U+FFFF; a surrogate with no
 * partner (a text is UTF-16 code units: a lone half is no character, and a serializer would break on it).
 */
const Forbidden = new RegExp(
  '[\\u0000-\\u0008\\u000b\\u000c\\u000e-\\u001f\\ufffe\\uffff]|[\\ud800-\\udbff](?![\\udc00-\\udfff])|(?<![\\ud800-\\udbff])[\\udc00-\\udfff]',
  'g',
);

/**
 * Text and attribute values escaped; a character XML forbids (a control but tab, newline, return, a lone surrogate)
 * becomes U+FFFD, one for one, so a glyph run stays aligned with its x list.
 */
export function escapeXml(s: string): string {
  return s.replace(Forbidden, '\ufffd').replace(/[&<>"']/g, (c) => XmlEscapes[c]);
}

/** `#rrggbb` and an opacity attribute when not opaque. */
function paint(attr: 'fill' | 'stroke', c: number): string {
  const rgb = `#${((c >>> 8) & 0xffffff).toString(16).padStart(6, '0')}`;
  const a = (c & 0xff) / 255;
  return a < 1 ? `${attr}="${rgb}" ${attr}-opacity="${fmt(a)}"` : `${attr}="${rgb}"`;
}

/** A font's family role (`sans`, `mono`) written as its stack; the rest of the shorthand as attributes. */
function fontAttrs(font: string, stacks: FontStacks): string {
  const size = /(\d+(?:\.\d+)?)px/.exec(font);
  const weight = /\b(bold|[1-9]00)\b/.exec(font);
  const family = /\bmono(space)?\b/i.test(font) ? stacks.mono : /\bsans\b/i.test(font) ? stacks.sans : font;
  let s = `font-family="${escapeXml(family)}" font-size="${size ? size[1] : '12'}"`;
  if (weight) s += ` font-weight="${weight[1]}"`;
  return s;
}

const Anchors = ['start', 'middle', 'end'];
const Commands = ['M', 'L', 'C', 'Z'];

function pathData(l: PathList, i: number): string {
  const d: string[] = [];
  let pt = l.pointStart[i];
  for (let c = l.cmdStart[i], end = c + l.cmdCount[i]; c < end; c++) {
    const cmd = l.commands[c];
    d.push(Commands[cmd]);
    for (let k = 0; k < PathCommandPoints[cmd]; k++, pt++)
      d.push(`${fmt(l.points[2 * pt])} ${fmt(l.points[2 * pt + 1])}`);
  }
  return d.join(' ');
}

export interface SvgOptions {
  /** The font stacks the family roles name (the token table's). */
  readonly fonts?: FontStacks;
  /** A title for the document (escaped). */
  readonly title?: string;
  /** Metadata as `data-*` attributes of the root (the open domains an export computed, whether it settled). */
  readonly metadata?: Readonly<Record<string, string>>;
  /** How far a text past its max width fades out, px (the layout's `fade`). */
  readonly fade?: number;
  /** Prefixes the document's ids (clips, masks), so two documents inlined in one page do not collide. */
  readonly idPrefix?: string;
  /**
   * The sections drawn: each an empty group with its address, direction, pieces and band (`data-section`,
   * `data-direction`, `data-pieces`, `data-y`, `data-height`), so the export says what the axis is.
   */
  readonly sections?: readonly {
    readonly address: string;
    readonly direction: string;
    readonly pieces: string;
    readonly y: number;
    readonly height: number;
  }[];
}

/** The draw lists as an SVG document of their size. */
export function drawListsToSvg(lists: DrawLists, options: SvgOptions = {}): string {
  const fonts = options.fonts ?? DefaultFonts;
  const prefix = (options.idPrefix ?? '').replace(/[^A-Za-z0-9_-]/g, '');
  const defs: string[] = [];
  const body: string[] = [];
  let ids = 0;
  const scratch = new PathList();
  const v: PassVisitor = {
    rect(l, i) {
      const fill = l.fill[i] === NoColor ? 'fill="none"' : paint('fill', l.fill[i]);
      const stroke =
        l.stroke[i] === NoColor ? '' : ` ${paint('stroke', l.stroke[i])} stroke-width="${fmt(l.strokeWidth[i])}"`;
      if (l.radius[i] > 0) {
        scratch.n = 0;
        scratch.commandCount = 0;
        scratch.pointCount = 0;
        scratch.begin();
        squirclePath(scratch, l.x[i], l.y[i], l.w[i], l.h[i], l.radius[i], l.round[i]);
        scratch.end(l.fill[i]);
        body.push(`<path d="${pathData(scratch, 0)}" ${fill}${stroke}/>`);
      } else
        body.push(
          `<rect x="${fmt(l.x[i])}" y="${fmt(l.y[i])}" width="${fmt(l.w[i])}" height="${fmt(l.h[i])}" ${fill}${stroke}/>`,
        );
    },
    path(l, i) {
      const fill = l.fill[i] === NoColor ? 'fill="none"' : paint('fill', l.fill[i]);
      let s = `<path d="${pathData(l, i)}" ${fill}`;
      if (l.stroke[i] !== NoColor) s += ` ${paint('stroke', l.stroke[i])} stroke-width="${fmt(l.strokeWidth[i])}"`;
      if (l.dash[i] > 0) s += ` stroke-dasharray="${fmt(l.dash[i])}"`;
      body.push(`${s}/>`);
    },
    glyphs(l, i) {
      const font = fontAttrs(lists.string(l.font[i]), fonts);
      const y = fmt(l.y[i]);
      const start = l.start[i];
      const end = start + l.count[i];
      if (l.placed[i] === 1) {
        // a placed run: a text a letter, moved to its place, turned and scaled about it
        for (let j = start; j < end; j++) {
          const color = l.color[i] !== NoColor ? l.color[i] : l.colors[j];
          const turn = l.angle[j] !== 0 ? ` rotate(${fmt((l.angle[j] * 180) / Math.PI)})` : '';
          const scale = l.vs[j] !== 1 ? ` scale(1 ${fmt(l.vs[j])})` : '';
          body.push(
            `<text transform="translate(${fmt(l.x[j])} ${fmt(l.ly[j])})${turn}${scale}" text-anchor="middle" ${font} ${paint('fill', color)}>${escapeXml(String.fromCharCode(l.codes[j]))}</text>`,
          );
        }
        return;
      }
      if (l.color[i] !== NoColor) {
        const xs: string[] = [];
        let text = '';
        for (let j = start; j < end; j++) {
          xs.push(fmt(l.x[j]));
          text += String.fromCharCode(l.codes[j]);
        }
        body.push(
          `<text x="${xs.join(' ')}" y="${y}" text-anchor="middle" ${font} ${paint('fill', l.color[i])}>${escapeXml(text)}</text>`,
        );
        return;
      }
      const spans: string[] = [];
      for (let j = start; j < end; j++)
        spans.push(
          `<tspan x="${fmt(l.x[j])}" ${paint('fill', l.colors[j])}>${escapeXml(String.fromCharCode(l.codes[j]))}</tspan>`,
        );
      body.push(`<text y="${y}" text-anchor="middle" ${font}>${spans.join('')}</text>`);
    },
    text(l, i) {
      const anchor = Anchors[l.anchor[i]];
      let mask = '';
      if (l.width[i] > l.maxWidth[i]) {
        // faded past its max width, never ellipsized: a gradient from opaque to clear over the last `fade` px
        const id = `${prefix}f${ids++}`;
        const x0 = anchor === 'start' ? l.x[i] : anchor === 'end' ? l.x[i] - l.maxWidth[i] : l.x[i] - l.maxWidth[i] / 2;
        const x1 = x0 + l.maxWidth[i];
        const fade = Math.min(options.fade ?? 12, l.maxWidth[i] / 2);
        defs.push(
          `<linearGradient id="${id}g" gradientUnits="userSpaceOnUse" x1="${fmt(x1 - fade)}" y1="0" x2="${fmt(x1)}" y2="0">` +
            '<stop offset="0" stop-color="#fff"/><stop offset="1" stop-color="#fff" stop-opacity="0"/></linearGradient>' +
            `<mask id="${id}" maskUnits="userSpaceOnUse"><rect x="${fmt(x0)}" y="${fmt(l.y[i] - 1000)}" width="${fmt(l.maxWidth[i])}" height="2000" fill="url(#${id}g)"/></mask>`,
        );
        mask = ` mask="url(#${id})"`;
      }
      body.push(
        `<text x="${fmt(l.x[i])}" y="${fmt(l.y[i])}" text-anchor="${anchor}" ${fontAttrs(lists.string(l.font[i]), fonts)} ${paint('fill', l.color[i])}${mask}>${escapeXml(lists.string(l.str[i]))}</text>`,
      );
    },
    image(l, i) {
      // one `<image>` a tile, its pixels a PNG data URL, never smoothed (a sample a device pixel or more)
      const img = l.images[i];
      body.push(
        `<image x="${fmt(l.x[i])}" y="${fmt(l.y[i])}" width="${fmt(l.w[i])}" height="${fmt(l.h[i])}" preserveAspectRatio="none" style="image-rendering:pixelated" href="${pngDataUrl(img.width, img.height, img.data)}"/>`,
      );
    },
    clip(l, i) {
      const id = `${prefix}c${ids++}`;
      defs.push(
        `<clipPath id="${id}"><rect x="${fmt(l.x[i])}" y="${fmt(l.y[i])}" width="${fmt(l.w[i])}" height="${fmt(l.h[i])}"/></clipPath>`,
      );
      body.push(`<g clip-path="url(#${id})">`);
    },
    unclip() {
      body.push('</g>');
    },
  };
  for (const sec of options.sections ?? [])
    body.push(
      `<g data-section="${escapeXml(sec.address)}" data-direction="${escapeXml(sec.direction)}" data-pieces="${escapeXml(sec.pieces)}" data-y="${fmt(sec.y)}" data-height="${fmt(sec.height)}"/>`,
    );
  for (const name of Passes) {
    const pass = lists.passes[name];
    if (!pass.size) continue;
    body.push(`<g data-pass="${name}">`);
    walkPass(pass, v);
    for (let d = pass.depth; d > 0; d--) body.push('</g>');
    body.push('</g>');
  }
  const w = fmt(lists.width);
  const h = fmt(lists.height);
  let meta = '';
  for (const [k, value] of Object.entries(options.metadata ?? {}))
    if (/^[a-z][a-z0-9-]*$/.test(k)) meta += ` data-${k}="${escapeXml(value)}"`;
  const title = options.title === undefined ? '' : `<title>${escapeXml(options.title)}</title>`;
  const defsBlock = defs.length ? `<defs>${defs.join('')}</defs>` : '';
  return `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}"${meta}>${title}${defsBlock}${body.join('')}</svg>`;
}
