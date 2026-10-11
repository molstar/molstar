/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mockScene } from '../../../_stories/mocks/scene.js';
import { render } from '../../../dom/render.js';
import { ShellModel } from '../../../shell/shell.js';
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
    // The mapped color param: switching to Uniform shows its own params (a color).
    choose('Color', 'Uniform');
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
