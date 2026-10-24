/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui/dom/view';
import { render } from '@molstar/prism-ui/dom/render';
import { flush, flushFrame } from '@molstar/prism-ui/dom/scope';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { afterEach, describe, expect, it } from 'vitest';
import type { SceneModel, SceneNode, SceneNodeKind } from '../../../model/scene.js';
import type { ViewerUIContext } from '../../../spec.js';
import { extensionUIPort } from '../../../ui/extensions.js';
import { SceneRow } from '../scene-row.js';
import { SceneTree } from '../scene-tree.js';
import { shownInOrder } from '@molstar/prism-ui/testing/virtual-list';

// The scene tree virtualized (view-models.md §3): a structure of 1,000 components, only the window in the DOM; ARIA's
// level, place and open state per row (an extension's row too); the keyboard across rows that are not rendered; Show
// all as a row of its own.

afterEach(() => {
  document.body.innerHTML = '';
  highlighted.clear();
});

/** The nodes highlighted now (hover). */
const highlighted = new Set<string>();

function node(id: string, kind: SceneNodeKind, label: string, children: SceneNode[] = [], expanded = false) {
  const isExpanded = observableState(expanded);
  const kids = observableState<readonly SceneNode[]>(children);
  return {
    id,
    kind,
    icon: 'structure',
    label: observableState(label).asReadonly,
    description: observableState<string | undefined>(undefined).asReadonly,
    isHidden: observableState(false).asReadonly,
    isExpanded: isExpanded.asReadonly,
    children: kids.asReadonly,
    setChildren: (c: SceneNode[]) => kids.set(c),
    toggleVisibility: () => {},
    toggleExpanded: () => isExpanded.set(!isExpanded.value),
    highlight: (on: boolean) => (on ? highlighted.add(id) : highlighted.delete(id)),
    focus: () => {},
    remove: () => {},
  } as unknown as SceneNode & { setChildren(c: SceneNode[]): void };
}

function mount(withShell = false) {
  const components = Array.from({ length: 1000 }, (_, i) =>
    node(`c${i}`, 'component', `Component ${i}`, [node(`c${i}r`, 'representation', 'Cartoon')]),
  );
  const entities = Array.from({ length: 12 }, (_, i) => node(`e${i}`, 'entity', `Entity ${i}`));
  const contents = node('contents', 'contents', 'Contents', entities, false);
  const structure = node('s', 'structure', '1ABC', [contents, ...components], true);
  const picked = observableState<SceneNode | undefined>(undefined);
  const scene = {
    roots: observableState([structure]).asReadonly,
    picked: picked.asReadonly,
    pick: (n: SceneNode | undefined) => picked.set(n),
  } as unknown as SceneModel;
  const shell = withShell ? new ShellModel() : undefined;
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <SceneTree scene={scene} shell={shell} />, container);
  const tree = container.querySelector<HTMLElement>('[role=tree]')!;
  const item = (key: string) =>
    [...tree.querySelectorAll<HTMLElement>('[data-tree-key]')]
      .find((i) => i.dataset.treeKey === key)
      ?.closest<HTMLElement>('[data-tree-key]');
  const treeitem = (key: string) => {
    const el = item(key);
    return el?.matches('[role=treeitem]') ? el : el?.querySelector<HTMLElement>('[role=treeitem]');
  };
  const key = (k: string) =>
    (document.activeElement as HTMLElement).dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));
  const focused = () =>
    (document.activeElement as HTMLElement).closest<HTMLElement>('[data-tree-key]')?.dataset.treeKey;
  let list: HTMLElement | undefined;
  for (let e: HTMLElement | null = tree; e; e = e.parentElement)
    if (e.scrollHeight >= 0 && e.getAttribute('class')?.includes('viewport')) list = e;
  const scrollTo = (index: number) => {
    // As the browser does: the viewport scrolls, then its scroll event.
    list!.scrollTop = index * 26;
    list!.dispatchEvent(new Event('scroll'));
    // The scroll position is input, published in the frame (rendering.md §2.3).
    flushFrame();
  };
  return {
    structure,
    contents,
    components,
    scene,
    picked,
    shell,
    tree,
    item,
    treeitem,
    key,
    focused,
    scrollTo,
    dispose,
  };
}

describe('the scene tree, virtualized', () => {
  it('renders only the rows on screen of a structure with 1,000 components', () => {
    const t = mount();
    const shown = t.tree.querySelectorAll('[role=treeitem]').length;
    expect(shown).toBeGreaterThan(20);
    expect(shown).toBeLessThan(100);
    expect(t.item('c999')).toBeUndefined();
    t.dispose();
  });

  it('gives each row its level, place among its siblings and open state', () => {
    const t = mount();
    const s = t.treeitem('s')!;
    expect([s.getAttribute('aria-level'), s.getAttribute('aria-posinset'), s.getAttribute('aria-setsize')]).toEqual([
      '1',
      '1',
      '1',
    ]);
    expect(s.getAttribute('aria-expanded')).toBe('true');
    const c3 = t.treeitem('c3')!;
    // Contents, then the components: the fourth component is the fifth of 1,001.
    expect([c3.getAttribute('aria-level'), c3.getAttribute('aria-posinset'), c3.getAttribute('aria-setsize')]).toEqual([
      '2',
      '5',
      '1001',
    ]);
    expect(c3.getAttribute('aria-expanded')).toBe('false');
    expect(t.tree.querySelector('[role=group]')).toBeNull();
    t.dispose();
  });

  it('moves the focus to rows not in the DOM, opens, closes and picks', () => {
    const t = mount();
    t.tree.focus();
    expect(t.focused()).toBe('s');
    t.key('End');
    expect(t.focused()).toBe('c999');
    expect(t.item('s')).toBeUndefined();
    t.key('ArrowRight');
    expect(t.treeitem('c999')!.getAttribute('aria-expanded')).toBe('true');
    t.key('ArrowRight');
    expect(t.focused()).toBe('c999r');
    t.key('ArrowLeft');
    expect(t.focused()).toBe('c999');
    t.key('ArrowLeft');
    expect(t.treeitem('c999')!.getAttribute('aria-expanded')).toBe('false');
    t.key('ArrowLeft');
    expect(t.focused()).toBe('s');
    t.key('End');
    t.key('Enter');
    expect(t.picked.value?.id).toBe('c999');
    t.dispose();
  });

  it('shows a long list a few at first, Show all a row of its own that the keyboard reaches', async () => {
    const t = mount();
    t.contents.toggleExpanded();
    flush();
    const rows = () => shownInOrder(t.tree, '[data-tree-key]').map((e) => e.dataset.treeKey);
    expect(rows().slice(0, 9)).toEqual(['s', 'contents', 'e0', 'e1', 'e2', 'e3', 'e4', 'contents::more', 'c0']);
    const more = t.treeitem('contents::more')!;
    expect([more.getAttribute('aria-posinset'), more.getAttribute('aria-setsize')]).toEqual(['6', '6']);
    t.tree.focus();
    t.key('ArrowDown');
    for (let i = 0; i < 6; i++) t.key('ArrowDown');
    expect(t.focused()).toBe('contents::more');
    t.key('Enter');
    expect(rows().slice(0, 15)).toContain('e11');
    // The focus goes to the row now where Show all was.
    await Promise.resolve();
    expect(t.focused()).toBe('e5');
    t.dispose();
  });

  it('gives an extension row its place, and keeps it one row (no controls under it)', () => {
    const t = mount(true);
    const port = extensionUIPort(t.shell!, {} as never, {} as ViewerUIContext);
    port.register('x', {
      rows: [
        {
          id: 'x.row.membrane',
          when: (n) => n.kind === 'structure',
          render: domView(() => <SceneRow label="Membrane" body={() => <p class="controls">Planes</p>} />),
        },
      ],
    });
    t.tree.focus();
    t.key('End');
    expect(t.focused()).toBe('s::row:x.row.membrane');
    const row = t.treeitem('s::row:x.row.membrane')!;
    expect([
      row.getAttribute('aria-level'),
      row.getAttribute('aria-posinset'),
      row.getAttribute('aria-setsize'),
    ]).toEqual(['2', '1002', '1002']);
    expect(row.hasAttribute('aria-expanded')).toBe(false);
    // No host to take its pick: still nothing opens under it in the tree.
    t.key('Enter');
    t.key('ArrowRight');
    expect(t.tree.querySelector('.controls')).toBeNull();
    t.key('ArrowLeft');
    expect(t.focused()).toBe('s');
    t.dispose();
  });

  it('a row pointed at another node as the tree scrolls: the old highlight goes, its menu closes, the picked follows', () => {
    const t = mount();
    const row = t.treeitem('c3')!;
    const line = row.querySelector<HTMLElement>('[data-tree-row]')!;
    line.dispatchEvent(new PointerEvent('pointerenter'));
    expect([...highlighted]).toEqual(['c3']);
    row.querySelector<HTMLButtonElement>('button[data-tip=More]')!.click();
    expect(row.querySelector('button[data-tip=More]')!.getAttribute('aria-expanded')).toBe('true');
    t.picked.set(t.components[3]);
    flush();
    expect(row.getAttribute('aria-selected')).toBe('true');
    // Far away: the slot shows another node under the resting pointer (highlighted then), or none (put away).
    t.scrollTo(600);
    const now = row.dataset.treeKey;
    expect(now).not.toBe('c3');
    expect(highlighted.has('c3')).toBe(false);
    expect([...highlighted]).toEqual(now ? [now] : []);
    expect(row.querySelector('button[data-tip=More]')!.getAttribute('aria-expanded')).toBe('false');
    expect(t.tree.querySelector('[role=menu]')).toBeNull();
    if (now) {
      expect(row.getAttribute('aria-selected')).toBe('false');
      expect(row.textContent).toContain(`Component ${now.slice(1)}`);
    } else expect(row.closest('[hidden]')).not.toBeNull();
    // Back: the picked row is shown as picked again, by whichever slot shows it.
    t.scrollTo(0);
    flush();
    expect(t.treeitem('c3')!.getAttribute('aria-selected')).toBe('true');
    t.dispose();
  });
});
