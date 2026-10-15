/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { mockScene } from '../../../testing/scene.js';
import { render } from '@molstar/prism-ui/dom/render';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { StructureTab } from '../structure-tab.js';

function mount() {
  const scene = mockScene();
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <StructureTab shell={new ShellModel()} scene={scene} />, container);
  const rows = () => [...container.querySelectorAll('[role=treeitem] > div:first-child')] as HTMLElement[];
  const row = (text: string) => rows().find((r) => r.querySelector('span')?.textContent === text)!;
  return { scene, container, dispose, rows, row };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('StructureTab', () => {
  it('renders the scene as a tree; hidden nodes keep their eye; collapsing hides children', () => {
    const { rows, row, scene } = mount();
    expect(row('1TQN')).toBeDefined();
    expect(row('Water').querySelector('button[title=Show]')).not.toBeNull();
    const before = rows().length;
    scene.roots.value[0].toggleExpanded();
    expect(rows().length).toBeLessThan(before);
  });

  it('offers to show a hidden node in its inspector, not only by its eye', () => {
    const { container, row } = mount();
    row('Water').click();
    const inspector = () => container.querySelector('[aria-label="Inspector"]')!;
    const show = [...inspector().querySelectorAll('button')].find((b) => b.textContent === 'Show')!;
    expect(show).toBeDefined();
    show.click();
    expect(row('Water').querySelector('button[title=Hide]')).not.toBeNull();
    expect([...inspector().querySelectorAll('button')].some((b) => b.textContent === 'Show')).toBe(false);
  });

  it('picks a node into the inspector, and form edits change the tree', () => {
    const { container, row } = mount();
    row('Cartoon').click();
    // A choice expands inline: open it, pick an option.
    const choose = (label: string, option: string) => {
      const button = container.querySelector(
        `button[aria-haspopup=listbox][aria-label="${label}"]`,
      ) as HTMLButtonElement;
      button.click();
      const pick = [...container.querySelectorAll('[role=listbox] [role=option]')].find(
        (o) => o.textContent === option,
      ) as HTMLButtonElement;
      pick.click();
    };
    const type = container.querySelector('button[aria-haspopup=listbox][aria-label=Type]') as HTMLButtonElement;
    expect(type.textContent).toBe('Cartoon');
    choose('Type', 'Spacefill');
    expect(row('Spacefill')).toBeDefined();
    expect(container.querySelector('[role=listbox]')).toBeNull();
    // The mapped color param: switching to Uniform, its sliders button opens its own params (a color; v5's).
    choose('Color', 'Uniform');
    (container.querySelector('button[aria-label="Color settings"]') as HTMLButtonElement).click();
    // A color row: the swatch button opens the palette under it.
    const swatch = [...container.querySelectorAll('[aria-label=Inspector] button[aria-expanded]')].find(
      (b) => b.textContent?.trim() === '46, 134, 193',
    ) as HTMLButtonElement;
    expect(swatch.textContent?.trim()).toBe('46, 134, 193');
    swatch.click();
    expect(container.querySelectorAll('[role=listbox][aria-label=Palette] [role=option]').length).toBe(8);
  });

  it('toggles visibility inline without picking', () => {
    const { row, scene } = mount();
    const polymer = row('Polymer');
    (polymer.querySelector('button[title=Hide]') as HTMLButtonElement).click();
    expect(row('Polymer').querySelector('button[title=Show]')).not.toBeNull();
    expect(scene.picked.value).toBeUndefined();
  });
});
