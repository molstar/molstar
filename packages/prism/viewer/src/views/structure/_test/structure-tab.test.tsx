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
import { domView } from '@molstar/prism-ui/dom/view';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { SceneRow } from '../scene-row.js';
import { SceneRowSlot } from '../scene-tree.js';

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

  it("shows and hides a node from its inspector's header, not only by its eye", () => {
    const { container, row } = mount();
    row('Water').click();
    const eye = () => container.querySelector<HTMLButtonElement>('[aria-label="Inspector"] header button')!;
    expect(eye().title).toBe('Show');
    expect(eye().getAttribute('aria-pressed')).toBe('true');
    eye().click();
    expect(row('Water').querySelector('button[title=Hide]')).not.toBeNull();
    expect(eye().title).toBe('Hide');
    eye().click();
    expect(row('Water').querySelector('button[title=Show]')).not.toBeNull();
  });

  it("shows an extension's row in the inspector when picked (its controls, turned on and off at its header), not under it", () => {
    const scene = mockScene();
    const shell = new ShellModel();
    const shown = observableState(false);
    shell.registerSections([
      {
        id: 'x.row.membrane',
        slot: SceneRowSlot,
        title: 'Membrane',
        when: (ctx) => (ctx.subject as { kind: string }).kind === 'structure',
        render: domView(() => (
          <SceneRow
            label="Membrane"
            description={observableState('Off').asReadonly}
            shown={shown.asReadonly}
            onToggle={() => shown.set(!shown.value)}
            body={() => <button type="button">Fit…</button>}
          />
        )),
      },
    ]);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <StructureTab shell={shell} scene={scene} />, container);
    const row = [...container.querySelectorAll<HTMLElement>('[data-tree-row]')].find((r) =>
      r.textContent?.startsWith('Membrane'),
    )!;
    row.click();
    const inspector = container.querySelector('[aria-label="Inspector"]')!;
    expect(inspector.textContent).toContain('Membrane');
    expect(row.parentElement!.textContent).not.toContain('Fit…');
    expect(inspector.textContent).toContain('Fit…');
    const eye = inspector.querySelector<HTMLButtonElement>('header button')!;
    expect(eye.title).toBe('Turn on');
    eye.click();
    expect(shown.value).toBe(true);
    expect(eye.title).toBe('Turn off');
    eye.click();
    expect(shown.value).toBe(false);
    // A node picked takes its place.
    scene.pick(scene.roots.value[0]);
    expect(container.querySelector('[aria-label="Inspector"]')!.textContent).not.toContain('Fit…');
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
