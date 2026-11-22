/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { BuilderCases, type BuilderCase, type Panel } from './builder-compose.js';

// The builder (`createSeqSpec`), each story the code beside the view it builds: the code is the source of the module
// the view is built from (Vite's `?raw`), so it cannot drift. Light and dark by the toolbar's theme.

const meta: Meta = { title: 'Features/Builder' };
export default meta;

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

const story = (id: string): StoryObj => {
  const c = BuilderCases.find((x) => x.id === id)!;
  return { render: () => builderStory(c) };
};

/** The smallest spec that draws: a sequence from letters, one track, one layer. */
export const FirstView: StoryObj = story('first-view');

/**
 * Every call returns the address of what it adds, typed by its kind: the sequence's handle is the ranges' target, the
 * ranges' handle the layer's data. The last line is a type error, kept under `@ts-expect-error` so the typecheck proves it.
 */
export const HandlesAndReferences: StoryObj = story('handles');

/** `theme(...)` and `measure(...)`: a scheme, a categorical field, rules over a value. */
export const LayersAndThemes: StoryObj = story('themes');

/** An assembly with the genome builder, a genome section of it, and one chromosome of it as the locus. */
export const AGenomeView: StoryObj = story('genome');

/** A spec with a mistake the types cannot catch: `check()`'s issues are listed beside the view. */
export const Check: StoryObj = story('check');

/** `build()` as data: the normalized tree printed as JSON and as the S-expression, beside the view. */
export const BuildAsData: StoryObj = story('as-data');
