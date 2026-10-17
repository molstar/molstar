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

// What a component attaches imperatively is released by construction (rendering.md §2.4): the scoped helpers
// (dom/scoped.ts: `on`, `observeResize`, `observeIntersection`, `frame`, `frames`, `timeout`, `drag`) are the way. A raw
// `addEventListener`, `new ResizeObserver` or `IntersectionObserver`, `setInterval`, `requestAnimationFrame` or a frame
// callback asked for by hand (`requestFrameCallback`: pooled, but not released with its scope) in a
// component file (the `.tsx` of ui, viewer and assistant, and the extensions' Prism UI) fails. Strict (2026-10-17: 93
// in 30 files at the start, all moved to the helpers).

const prism = path.resolve(import.meta.dirname, '../../..');
const repo = path.resolve(prism, '../..');
const Raw =
  /\.addEventListener\(|new (Resize|Intersection)Observer\(|\bsetInterval\(|\brequestAnimationFrame\(|\brequestFrameCallback\(/g;

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

  it('no raw attachment in a component', () => {
    const found: string[] = [];
    for (const f of componentFiles()) {
      const n = rawAttach(fs.readFileSync(f, 'utf8'));
      if (n > 0) found.push(`${path.relative(repo, f).split(path.sep).join('/')}: ${n}`);
    }
    expect(found, 'a raw attachment in a component: use dom/scoped.ts (rendering.md §2.4)').toEqual([]);
  });
});
