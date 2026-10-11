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
    const type = container.querySelector('select') as HTMLSelectElement;
    expect(type.options[type.selectedIndex].text).toBe('Cartoon');
    type.value = String([...type.options].findIndex((o) => o.text === 'Spacefill'));
    type.dispatchEvent(new Event('change'));
    expect(row('Spacefill')).toBeDefined();
    // The mapped color param: switching to Uniform shows its own params (a color).
    const color = container.querySelectorAll('select')[1] as HTMLSelectElement;
    color.value = String([...color.options].findIndex((o) => o.text === 'Uniform'));
    color.dispatchEvent(new Event('change'));
    // A color row: the swatch button opens the palette under it.
    const swatch = [...container.querySelectorAll('[aria-label=Inspector] button[aria-expanded]')].find((b) =>
      b.textContent?.startsWith('#'),
    ) as HTMLButtonElement;
    expect(swatch.textContent).toBe('#2e86c1');
    swatch.click();
    expect(container.querySelectorAll('[role=option]').length).toBeGreaterThan(8);
  });

  it('toggles visibility inline without picking', () => {
    const { row, scene } = mount();
    const polymer = row('Polymer');
    (polymer.querySelector('button[title=Hide]') as HTMLButtonElement).click();
    expect(row('Polymer').querySelector('button[title=Show]')).not.toBeNull();
    expect(scene.picked.value).toBeUndefined();
  });
});
