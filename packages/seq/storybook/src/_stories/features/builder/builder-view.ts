/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { exampleStory } from '../../examples/example-root.js';
import type { BuilderCase, Panel } from './builder-compose.js';

// The Builder stories' page: the code, the view it builds and any panels (not a stories file: every export of one is a
// story).

const Pre =
  'margin: 0; padding: var(--msui-space-3); overflow: auto; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line); color: var(--msui-text);';

function titled(title: string, text: string): HTMLElement {
  const box = document.createElement('div');
  box.style.cssText = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
  const h = document.createElement('h4');
  h.style.cssText =
    'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';
  h.textContent = title;
  const pre = document.createElement('pre');
  pre.style.cssText = Pre;
  pre.textContent = text;
  box.append(h, pre);
  return box;
}

/** The code, the view the code builds, and any panels from the same spec. */
export function builderStory(c: BuilderCase, panels: readonly Panel[] = c.panels()): HTMLElement {
  const root = document.createElement('div');
  root.style.cssText =
    'display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5); align-items: start; background: var(--msui-bg);';
  const left = document.createElement('div');
  left.style.cssText = 'display: grid; gap: var(--msui-space-5); padding: var(--msui-space-5); min-width: 0;';
  left.append(titled('The code', c.source), ...panels.map((p) => titled(p.title, p.text)));
  root.append(left, exampleStory(c));
  return root;
}
