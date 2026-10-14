/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small XML reader for the SVG tests: elements with their attributes, children and text, entities decoded. It throws
// on what a well-formed document cannot have (a mismatched or unclosed tag, an unknown entity, an unquoted attribute).

export interface El {
  name: string;
  attrs: Record<string, string>;
  children: El[];
  text: string;
}
const ENTITIES: Record<string, string> = { amp: '&', lt: '<', gt: '>', quot: '"', apos: "'" };
const decode = (s: string) =>
  s.replace(/&(\w+);/g, (_, n) => {
    if (!(n in ENTITIES)) throw new Error(`unknown entity &${n};`);
    return ENTITIES[n];
  });

export function parseXml(xml: string): El {
  const root: El = { name: '#root', attrs: {}, children: [], text: '' };
  const stack = [root];
  const tag = /<(\/?)([\w:-]+)((?:\s+[\w:-]+="[^"]*")*)\s*(\/?)>|([^<]+)/gy;
  let at = 0;
  for (let m = tag.exec(xml); m; m = tag.exec(xml)) {
    at = tag.lastIndex;
    const top = stack[stack.length - 1];
    if (m[5] !== undefined) top.text += decode(m[5]);
    else if (m[1]) {
      if (top.name !== m[2]) throw new Error(`</${m[2]}> closes <${top.name}>`);
      stack.pop();
    } else {
      const el: El = { name: m[2], attrs: {}, children: [], text: '' };
      for (const a of m[3].matchAll(/([\w:-]+)="([^"]*)"/g)) el.attrs[a[1]] = decode(a[2]);
      top.children.push(el);
      if (!m[4]) stack.push(el);
    }
  }
  if (at !== xml.length) throw new Error(`not well formed at ${at}: ${xml.slice(at, at + 40)}`);
  if (stack.length !== 1) throw new Error(`<${stack[stack.length - 1].name}> is not closed`);
  if (root.children.length !== 1) throw new Error('one root element expected');
  return root.children[0];
}
/** Every descendant (and the element itself) named `name`, in document order. */
export function all(el: El, name: string, out: El[] = []): El[] {
  if (el.name === name) out.push(el);
  for (const c of el.children) all(c, name, out);
  return out;
}
/** The numbers of a path or list attribute, with the commands between them: ["M", 1, 2, "h", 3]. */
export function tokens(d: string): (string | number)[] {
  return [...d.matchAll(/[A-Za-z]|-?\d*\.?\d+/g)].map((m) => (/[A-Za-z]/.test(m[0]) ? m[0] : Number(m[0])));
}
/** The rectangles of a path made of `M x y h w v h h -w z` subpaths: [x, y, w, h]. */
export function rectsOf(d: string): number[][] {
  return [...d.matchAll(/M(-?[\d.]+) (-?[\d.]+)h(-?[\d.]+)v(-?[\d.]+)h-?[\d.]+z/g)].map((m) => m.slice(1).map(Number));
}
