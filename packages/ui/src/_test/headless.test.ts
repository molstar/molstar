/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs';
import { join, relative, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';

// Headless (the models drive the app without a document: an agent, a script, a headless plugin): every model file of
// the package is DOM-free. It imports no view (a `.tsx` module) and no runtime from `dom/` (types only), and touches no
// DOM global; and it loads in plain Node (this file runs without jsdom).

const src = resolve(process.cwd(), 'packages/ui/src');

/** The model files: `*-model.ts`, `*-commands.ts`, the shell's layout and pane map, the stage's ink. */
const isModel = (path: string) =>
  /(^|[/\\])[\w-]*model\.ts$/.test(path) ||
  /-commands\.ts$/.test(path) ||
  /shell[/\\](layout|pane-map)\.ts$/.test(path) ||
  /stage[/\\]canvas-ink\.ts$/.test(path);

function walk(dir: string, out: string[] = []): string[] {
  for (const name of readdirSync(dir)) {
    if (name === 'node_modules' || name.startsWith('_')) continue;
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walk(p, out);
    else if (p.endsWith('.ts') && isModel(p)) out.push(p);
  }
  return out;
}

/** The source without comments and strings (so a word in a doc comment is not a use). */
function code(text: string): string {
  return text
    .replace(/\/\*[\s\S]*?\*\//g, '')
    .replace(/\/\/.*$/gm, '')
    .replace(/'(?:[^'\\\n]|\\.)*'|"(?:[^"\\\n]|\\.)*"|`(?:[^`\\]|\\.)*`/g, "''");
}

/** A DOM global used (not a property of the same name: `place.window` is the layout's). */
const DomGlobals =
  /(?<![.\w$])(document|window|navigator|getComputedStyle|requestAnimationFrame|cancelAnimationFrame|matchMedia)\s*[.(]|\bnew\s+(ResizeObserver|IntersectionObserver|MutationObserver|Path2D)\b/;

const models = walk(src);

describe('headless: the model files are DOM-free', () => {
  it('finds the models', () => {
    expect(models.length).toBeGreaterThan(3);
  });

  it('catches a DOM global (positive control), not a property or a comment of the same name', () => {
    expect(code('const w = window.innerWidth;').match(DomGlobals)).not.toBeNull();
    expect(code('new ResizeObserver(f)').match(DomGlobals)).not.toBeNull();
    expect(code('getComputedStyle(el).color').match(DomGlobals)).not.toBeNull();
    expect(code('if (place.window) x = place.window.rect;').match(DomGlobals)).toBeNull();
    expect(code('// the document. and window.\nconst a = 1;').match(DomGlobals)).toBeNull();
  });

  for (const file of models) {
    const name = relative(src, file);
    it(`${name}: no view, no dom runtime, no DOM global`, async () => {
      const text = readFileSync(file, 'utf8');
      for (const m of text.matchAll(/^import\s+(type\s+)?[^;]*?from\s+'([^']+)';/gms)) {
        const typeOnly = !!m[1];
        const target = m[2];
        if (!target.startsWith('.')) continue;
        const resolved = resolve(file, '..', target.replace(/\.js$/, ''));
        // A view (a `.tsx` module) is never imported by a model, not even its types' module at runtime.
        if (!typeOnly) expect(existsSync(`${resolved}.tsx`), `${name} imports the view ${target}`).toBe(false);
        // The dom runtime: types only.
        if (/(^|\/)dom(\/|$)/.test(target)) expect(typeOnly, `${name} imports the dom runtime (${target})`).toBe(true);
      }
      const used = code(text).match(DomGlobals);
      expect(used?.[0], `${name} uses a DOM global`).toBeUndefined();
      // And it loads in plain Node.
      await import(file);
    });
  }
});
