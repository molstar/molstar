/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CommandMenu } from '@molstar/prism-ui';
import { render } from '@molstar/prism-ui/dom/render';
import type { PluginContext } from '@molstar/plugin/context';
import { afterEach, describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../../_test/harness.js';
import { createStyleForms } from '../../../sources/command-forms.js';
import type { PrismViewer } from '../../../viewer.js';
import type { ViewerCore } from '../../../viewer.js';
import { styleCommandForms } from '../style-forms.js';

// The command menu's own forms for Show as… and Color by theme… (commands.md §5.1), over a real headless viewer: opened
// as the inspector's rows open them, the keys of the menu's input, and what a run changes.

const flush = () => new Promise((r) => setTimeout(r, 0));

const pluginOf = (v: PrismViewer): PluginContext => (v as ViewerCore).plugin;

const polymer = (plugin: PluginContext) => {
  const s = plugin.managers.structure.hierarchy.current.structures[0];
  return s.components.find((c) => c.cell.obj?.label === 'Polymer')!;
};

const typesOf = (plugin: PluginContext) =>
  polymer(plugin).representations.map((r) => (r.cell.transform.params as { type: { name: string } }).type.name);

async function mount(id: string, values: Record<string, string>) {
  const v = await headlessViewer({ files: [cbs()] });
  const plugin = pluginOf(v);
  const manager = v.managers.commands;
  const container = document.createElement('div');
  document.body.append(container);
  const forms = styleCommandForms(createStyleForms(plugin));
  const opened = manager.open(id, undefined, values);
  if (opened.kind !== 'form') throw new Error(opened.kind);
  render(() => <CommandMenu manager={manager} forms={forms} direct onClose={() => {}} />, container);
  await flush();
  const input = () => container.querySelector('input') as HTMLInputElement;
  const key = (k: string) => input().dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const type = (text: string) => {
    input().value = text;
    input().dispatchEvent(new Event('input'));
  };
  const tiles = () => [...container.querySelectorAll('[role=radio]')] as HTMLElement[];
  const checked = () => tiles().find((t) => t.getAttribute('aria-checked') === 'true')?.textContent;
  return { v, plugin, invocation: opened.invocation, container, key, type, tiles, checked };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('Color by theme… (the bug: its form opened on the query language)', () => {
  it('opens on the theme, listing color themes, with no "None"', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const opened = v.managers.commands.open('representation.color-theme', undefined, {
      structure: '1CBS',
      target: 'polymer',
    });
    if (opened.kind !== 'form') throw new Error(opened.kind);
    const s = opened.invocation.state.value;
    // The theme is asked for: the form's active field, its options the themes.
    expect(s.values.theme).toBeUndefined();
    expect(s.missing).toEqual(['theme']);
    expect(s.field).toBe('theme');
    const listed = s.completions.map((c) => c.value);
    expect(listed).toContain('chain-id');
    expect(listed).not.toContain('');
    expect(listed).not.toContain('pymol');
    // The uniform theme's color is not a field of the sentence.
    expect(s.fields.map((f) => f.id)).not.toContain('color');
  });
});

describe('Show as… own form', () => {
  it('shows the types the component supports, a type it does not show picked; ↓ and Enter add it', async () => {
    const p = await mount('representation.set', { component: 'Polymer', structure: '1CBS' });
    const labels = p.tiles().map((t) => t.textContent ?? '');
    expect(labels.some((l) => l.startsWith('Cartoon'))).toBe(true);
    expect(labels.some((l) => l.startsWith('Spacefill'))).toBe(true);
    // Cartoon is shown already: the pick starts on another.
    expect(labels.find((l) => l.startsWith('Cartoon'))).toContain('shown');
    expect(p.checked()).not.toMatch(/^Cartoon/);
    // Typing filters; Enter takes the first match and runs.
    p.type('spacef');
    await flush();
    expect(p.tiles().map((t) => t.textContent)).toEqual(['Spacefill']);
    p.key('Enter');
    await flush();
    await flush();
    expect(typesOf(p.plugin)).toContain('spacefill');
  });

  it('runs with the type’s options set in the form', async () => {
    const p = await mount('representation.set', { component: 'Polymer', structure: '1CBS' });
    p.invocation.set('type', 'spacefill', { stay: true });
    await flush();
    // Size, quality and opacity of the picked type.
    expect(p.container.textContent).toContain('Size');
    expect(p.container.textContent).toContain('Quality');
    expect(p.invocation.setParam('params', { sizeFactor: 0.5 })).toBe(true);
    await p.invocation.confirm();
    const added = polymer(p.plugin).representations.find(
      (r) => (r.cell.transform.params as { type: { name: string } }).type.name === 'spacefill',
    );
    expect((added?.cell.transform.params as { type: { params: { sizeFactor: number } } }).type.params.sizeFactor).toBe(
      0.5,
    );
  });

  it('moves the pick with the arrows', async () => {
    const p = await mount('representation.set', { component: 'Polymer', structure: '1CBS' });
    const before = p.checked();
    p.key('ArrowRight');
    await flush();
    expect(p.checked()).not.toBe(before);
    p.key('ArrowLeft');
    await flush();
    expect(p.checked()).toBe(before);
  });
});

describe('Color by theme… own form', () => {
  it('applies a theme on a click; one color waits for its swatch', async () => {
    const p = await mount('representation.color-theme', { structure: '1CBS', target: 'polymer' });
    p.tiles()
      .find((t) => t.textContent === 'Uniform')!
      .click();
    await flush();
    expect(p.invocation.state.value.status).toBe('editing');
    p.tiles()
      .find((t) => t.textContent === 'Element Symbol')!
      .click();
    await flush();
    await flush();
    const r = polymer(p.plugin).representations[0];
    expect((r.cell.transform.params as { colorTheme: { name: string } }).colorTheme.name).toBe('element-symbol');
  });

  it('paints the themes that apply and runs the uniform theme with its color', async () => {
    const p = await mount('representation.color-theme', { structure: '1CBS', target: 'polymer' });
    const labels = p.tiles().map((t) => t.textContent ?? '');
    expect(labels.length).toBeGreaterThan(3);
    // Under the registry's categories, as the inspector's theme select lists them.
    const groups = [...p.container.querySelectorAll('[role=radiogroup] [role=group]')].map((g) =>
      g.getAttribute('aria-label'),
    );
    expect(groups).toContain('Chain Property');
    expect(groups).toContain('Atom Property');
    expect(p.tiles().some((t) => t.textContent === 'Uniform')).toBe(true);
    // A theme with a legend is drawn with its colors (jsdom keeps no gradient in a style: the model's swatch).
    const themes = createStyleForms(p.plugin).themes({ structure: '1CBS' });
    expect(themes.find((t) => t.id === 'chain-id')?.swatch).toMatch(/gradient/);
    p.invocation.set('theme', 'uniform', { stay: true });
    p.invocation.setParam('color', '#ff0000');
    await p.invocation.confirm();
    const r = polymer(p.plugin).representations[0];
    const theme = (r.cell.transform.params as { colorTheme: { name: string; params: { value: number } } }).colorTheme;
    expect(theme.name).toBe('uniform');
    expect(theme.params.value).toBe(0xff0000);
  });
});
