/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import { RawAttach } from './scoped-attach.offenders.js';

// What a component attaches imperatively is released by construction (rendering.md §2.4): the scoped helpers
// (dom/scoped.ts: `on`, `observeResize`, `observeIntersection`, `frame`, `timeout`) are the way. A raw
// `addEventListener`, `new ResizeObserver` or `IntersectionObserver`, `setInterval` or `requestAnimationFrame` in a
// component file (the `.tsx` of ui, viewer and assistant, and the extensions' Prism UI) counts; the files that still
// have them are listed with their counts (scoped-attach.offenders.ts), a new one or more in a listed one fails, and a
// listed one with fewer must lower its entry.

const prism = path.resolve(import.meta.dirname, '../../..');
const repo = path.resolve(prism, '../..');
const Raw = /\.addEventListener\(|new (Resize|Intersection)Observer\(|\bsetInterval\(|\brequestAnimationFrame\(/g;

function* tsx(dir: string): Generator<string> {
  if (!fs.existsSync(dir)) return;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (['_test', '_stories', '_generated', 'testing', 'node_modules'].includes(entry.name)) continue;
      yield* tsx(full);
    } else if (entry.name.endsWith('.tsx')) yield full;
  }
}

function componentFiles(): string[] {
  const out: string[] = [];
  for (const pkg of ['ui', 'viewer', 'assistant']) out.push(...tsx(path.join(prism, pkg, 'src')));
  const extensions = path.join(repo, 'extensions');
  for (const ext of fs.readdirSync(extensions))
    for (const f of tsx(path.join(extensions, ext, 'src')))
      if (fs.readFileSync(f, 'utf8').includes('@molstar/prism-ui')) out.push(f);
  return out;
}

export function rawAttach(source: string): number {
  return [...source.replace(/(^|[^:])\/\/.*$/gm, '$1').matchAll(Raw)].length;
}

describe('scoped attachments', () => {
  it('counts the raw ones', () => {
    expect(
      rawAttach("el.addEventListener('x', f);\nrequestAnimationFrame(g);\n// setInterval(h)\non(el, 'x', f);"),
    ).toBe(2);
  });

  it('no new raw attachment in a component; the listed ones only shrink', () => {
    const found: Record<string, number> = {};
    for (const f of componentFiles()) {
      const n = rawAttach(fs.readFileSync(f, 'utf8'));
      if (n > 0) found[path.relative(repo, f).split(path.sep).join('/')] = n;
    }
    const grown: string[] = [];
    const shrunk: string[] = [];
    for (const [file, n] of Object.entries(found)) {
      const allowed = RawAttach[file] ?? 0;
      if (n > allowed) grown.push(`${file}: ${n} (allowed ${allowed})`);
      else if (n < allowed) shrunk.push(`${file}: ${n} (listed ${allowed})`);
    }
    for (const [file, n] of Object.entries(RawAttach)) if (!(file in found)) shrunk.push(`${file}: 0 (listed ${n})`);
    expect(grown, 'a raw attachment in a component: use dom/scoped.ts (rendering.md §2.4)').toEqual([]);
    expect(shrunk, 'fewer than listed: lower or remove the entry in scoped-attach.offenders.ts').toEqual([]);
  });
});
