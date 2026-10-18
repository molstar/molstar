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
// callback asked for by hand (`requestFrameCallback`: pooled, but not released with its scope) in a component file
// fails: the `.tsx` of ui, viewer and assistant and the extensions' Prism UI, and the `.ts` view code beside them (ui's
// `components/` and `views/`, the viewer's and the assistant's `views/`, and the helpers views call in a scope:
// `dom/tree-keys.ts`, `dom/focus.ts`). Strict (2026-10-17: 93 in 30 files at the start, all moved to the helpers); a
// `.ts` that attaches for a reason is in `Exempt`, with it.

const prism = path.resolve(import.meta.dirname, '../../..');
const repo = path.resolve(prism, '../..');
const Raw =
  /\.addEventListener\(|new (Resize|Intersection)Observer\(|\bsetInterval\(|\brequestAnimationFrame\(|\brequestFrameCallback\(/g;

/** `.ts` view code that attaches by hand, and why (paths from the repository root). */
const Exempt: Record<string, string> = {
  'packages/prism/ui/src/components/markdown/markdown.ts':
    'listeners on the links it creates, collected with them; a link removed under the pointer is left by `onCleanup`',
  'packages/prism/ui/src/components/reorder/drag-to-reorder.ts':
    '`dragToReorder` is the press-attached form, for code with no scope (stories); components use `reorderRows`',
  'packages/prism/ui/src/views/shell/prism-field.ts':
    "a simulation's frame loop, started, paused and disposed by its component's scope (prism-splash.tsx)",
  'packages/prism/ui/src/views/shell/quantum-field.ts':
    "a simulation's frame loop, started, paused and disposed by its component's scope (prism-splash.tsx)",
};

function* files(dir: string, ext: '.tsx' | '.ts'): Generator<string> {
  if (!fs.existsSync(dir)) return;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (['_test', '_stories', '_generated', 'testing', 'node_modules'].includes(entry.name)) continue;
      yield* files(full, ext);
    } else if (entry.name.endsWith(ext) && !entry.name.endsWith('.d.ts')) yield full;
  }
}

function componentFiles(): string[] {
  const out: string[] = [];
  for (const pkg of ['ui', 'viewer', 'assistant']) out.push(...files(path.join(prism, pkg, 'src'), '.tsx'));
  for (const dir of ['ui/src/components', 'ui/src/views', 'viewer/src/views', 'assistant/src/views'])
    out.push(...files(path.join(prism, dir), '.ts'));
  out.push(path.join(prism, 'ui/src/dom/tree-keys.ts'), path.join(prism, 'ui/src/dom/focus.ts'));
  const extensions = path.join(repo, 'extensions');
  for (const ext of fs.readdirSync(extensions))
    for (const f of files(path.join(extensions, ext, 'src'), '.tsx'))
      if (fs.readFileSync(f, 'utf8').includes('@molstar/prism-ui')) out.push(f);
  return out;
}

const relative = (f: string) => path.relative(repo, f).split(path.sep).join('/');

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
      if (Exempt[relative(f)]) continue;
      const n = rawAttach(fs.readFileSync(f, 'utf8'));
      if (n > 0) found.push(`${relative(f)}: ${n}`);
    }
    expect(found, 'a raw attachment in a component: use dom/scoped.ts (rendering.md §2.4)').toEqual([]);
  });

  it('scans the .ts view code too, and every exemption is still needed', () => {
    const scanned = new Set(componentFiles().map(relative));
    expect(scanned.has('packages/prism/ui/src/components/tooltip/tooltip.ts')).toBe(true);
    expect(scanned.has('packages/prism/ui/src/dom/tree-keys.ts')).toBe(true);
    for (const f of Object.keys(Exempt)) {
      expect(scanned.has(f), `${f}: exempt but not scanned`).toBe(true);
      expect(rawAttach(fs.readFileSync(path.join(repo, f), 'utf8')), `${f}: exempt, attaches nothing`).toBeGreaterThan(
        0,
      );
    }
  });
});
