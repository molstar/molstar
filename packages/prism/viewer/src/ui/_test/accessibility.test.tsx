/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expectNoAxeViolations } from '@molstar/prism-ui/testing/axe';
import { afterAll, beforeAll, describe, expect, it, vi } from 'vitest';
import { ensemble } from '../../api/_test/trajectory-fixtures.js';
import { mapFile } from '../../api/_test/volume-fixtures.js';
import { cbs, expectOk } from '../../_test/harness.js';
import type { SceneNode } from '../../model/scene.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import type { PrismViewerSpec } from '../../spec.js';
import { createDefaultPlugin } from '../../sources/create-plugin.js';
import { createViewerCore, type ViewerCore } from '../../viewer.js';
import { type MountedShell, mountShell } from '../shell.js';

// Automated checks (axe-core, jsdom) of the page's shell over a viewer core with a structure loaded: every tab of the
// rail, the cards, the command menu, an open row menu (accessibility.md, Automated checks). The shell lives in a shadow
// root, as on a page.

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

const wait = (ms = 60) => new Promise((r) => setTimeout(r, ms));

/**
 * Waits until the page's DOM (the document and the shells' shadow roots) has not changed for `gap` ms, at most `cap`:
 * a view rendered, its late parts in (a lazy chunk, a source's debounce), without a fixed sleep per step.
 */
function quiet(gap = 80, cap = 1500): Promise<void> {
  return new Promise((resolve) => {
    const roots: Node[] = [document];
    for (const el of document.querySelectorAll('*')) if (el.shadowRoot) roots.push(el.shadowRoot);
    let timer: ReturnType<typeof setTimeout>;
    const done = () => {
      observer.disconnect();
      clearTimeout(timer);
      clearTimeout(limit);
      resolve();
    };
    const observer = new MutationObserver(() => {
      clearTimeout(timer);
      timer = setTimeout(done, gap);
    });
    for (const r of roots)
      observer.observe(r, { subtree: true, childList: true, attributes: true, characterData: true });
    timer = setTimeout(done, gap);
    const limit = setTimeout(done, cap);
  });
}

interface Page {
  core: ViewerCore;
  mounted: MountedShell;
  /** Waits, then fails on a violation in the shell. */
  check(): Promise<void>;
  dispose(): Promise<void>;
}

/** The page's shell over a core, the files loaded into it. */
const open = (...files: File[]) => openWith(DefaultPrismViewerSpec(), ...files);

async function openWith(spec: PrismViewerSpec, ...files: File[]): Promise<Page> {
  const core = await createViewerCore({
    spec,
    createPlugin: createDefaultPlugin,
    modules: ['files', 'page'],
  });
  // No WebGL in jsdom: the canvas does not start (the shell says so).
  vi.spyOn(core, 'mountCanvas').mockResolvedValue(false);
  const host = document.createElement('div');
  document.body.append(host);
  const mounted = await mountShell(core, host);
  for (const file of files) expectOk(await core.api.load.files({ files: [file] }));
  await quiet();
  return {
    core,
    mounted,
    check: async () => {
      await wait();
      await expectNoAxeViolations(mounted.root);
    },
    dispose: async () => {
      mounted.dispose();
      await core.dispose();
      host.remove();
    },
  };
}

beforeAll(() => {
  vi.stubGlobal('ResizeObserver', NoResize);
});
afterAll(() => {
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

// An axe check takes about a second: the self-check, the rail, the home tabs and the trajectory bar are the quick
// lane's representatives, the rest of the sweep is medium (AGENTS.md, Tests). A medium test may build on the state
// an earlier one left (snapshots, a selection, history steps): the quick ones leave none.
const medium = { tags: ['medium'] };

let core: ViewerCore;
let mounted: MountedShell;
let check: () => Promise<void>;
let page: Page;

describe('the viewer shell has no accessibility violations (axe)', () => {
  beforeAll(async () => {
    page = await open(cbs());
    ({ core, mounted, check } = page);
  });
  afterAll(() => page.dispose());

  it('sees into the shell (the check can fail)', async () => {
    const bad = document.createElement('button');
    mounted.root.append(bad);
    await expect(check()).rejects.toThrow(/button-name/);
    bad.remove();
  });

  it('the rail and the canvas chrome', async () => {
    expect(mounted.root.querySelectorAll('[role=tab]').length).toBeGreaterThan(4);
    await check();
  });

  for (const tab of [
    'home',
    'structure',
    'selection',
    'snapshots.tab.snapshots',
    'state-tree.tab.tree',
    'history',
    'settings',
    'help',
  ]) {
    it(`the ${tab} tab`, tab === 'home' ? {} : medium, async () => {
      mounted.shell.show(tab);
      await quiet();
      await check();
    });
  }

  it('the snapshots tab and the playback card with snapshots', medium, async () => {
    const person = core.api.as({ kind: 'person' });
    expectOk(await person.snapshot!.add({}));
    expectOk(await person.snapshot!.add({}));
    mounted.shell.show('snapshots.tab.snapshots');
    await quiet();
    expect(mounted.root.querySelectorAll('[role=listitem]').length).toBe(2);
    await check();
    mounted.shell.openOverlay('snapshots.overlay.playback');
    await quiet();
    await check();
    mounted.shell.closeOverlay();
  });

  it('the structure tab with a node picked (the inspector) and a row menu open', medium, async () => {
    mounted.shell.show('structure');
    const scene = core.views.scene;
    const first = scene.roots.value[0];
    expect(first).toBeDefined();
    scene.pick(first);
    await quiet();
    await check();
    const more = mounted.root.querySelector<HTMLElement>('[role=treeitem] button[data-tip=More]');
    expect(more).not.toBeNull();
    more!.click();
    await wait(100);
    expect(mounted.root.querySelector('[role=menu]')).not.toBeNull();
    await check();
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
  });

  it('the selection tab with a selection', medium, async () => {
    expectOk(await core.api.as({ kind: 'person' }).select.set({ target: 'ligand' }));
    mounted.shell.show('selection');
    await quiet();
    await check();
  });

  it('the state tree with a cell picked (its form and actions)', medium, async () => {
    mounted.shell.show('state-tree.tab.tree');
    await quiet();
    const row = mounted.root.querySelector<HTMLElement>('[role=treeitem] > div:first-child');
    expect(row).not.toBeNull();
    row!.click();
    await quiet();
    await check();
  });

  it('the history tab with steps', medium, async () => {
    expectOk(await core.api.as({ kind: 'person' }).settings.set({ group: 'viewer', values: { checkpointing: false } }));
    mounted.shell.show('history');
    await quiet();
    expect(mounted.root.querySelector('[data-origin]')).not.toBeNull();
    await check();
  });

  for (const card of ['export', 'canvas-settings.overlay.canvas', 'animation.overlay.animation']) {
    it(`the ${card} card`, medium, async () => {
      mounted.shell.openOverlay(card);
      await quiet();
      expect(mounted.root.querySelector('[role=dialog]')).not.toBeNull();
      await check();
      mounted.shell.closeOverlay();
    });
  }

  it('the sequence band', medium, async () => {
    mounted.shell.show('sequence');
    await quiet();
    expect(mounted.root.querySelector('[aria-label="What a click selects"]')).not.toBeNull();
    await check();
  });

  it('the command menu open', medium, async () => {
    // From inside the shell: the page's keys belong to the viewer whose shell has the focus.
    mounted.root
      .querySelector('[role=tab]')!
      .dispatchEvent(new KeyboardEvent('keydown', { key: 'k', ctrlKey: true, bubbles: true, composed: true }));
    await quiet();
    expect(mounted.root.querySelector('[role=combobox]')).not.toBeNull();
    await check();
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
  });

  it('a compact shell', medium, async () => {
    mounted.shell.setMode('compact');
    mounted.shell.show('home');
    await quiet();
    await check();
    mounted.shell.setMode('desktop');
  });
});

describe('the viewer shell with a trajectory and a volume (axe)', () => {
  let page: Page;
  beforeAll(async () => {
    page = await open(ensemble(), mapFile());
  });
  afterAll(() => page.dispose());

  it("shows the canvas' full trajectory bar unless the trajectory's row is on screen (another tab, the canvas expanded)", async () => {
    const bar = () => {
      const g = page.mounted.root.querySelector<HTMLElement>('[role=group][aria-label=Models]');
      return g ? (g.querySelector('input[type=range]') ? 'full' : 'compact') : 'none';
    };
    page.mounted.shell.show('structure');
    await quiet();
    expect(page.mounted.root.querySelector('[data-scene-row="scene.trajectory"]')).not.toBeNull();
    // The row has the controls, paused: no bar on the canvas.
    expect(bar()).toBe('none');
    page.mounted.shell.show('selection');
    await quiet();
    expect(bar()).toBe('full');
    page.mounted.shell.show('structure');
    await quiet();
    expect(bar()).toBe('none');
    page.mounted.shell.setExpanded(true);
    await quiet();
    expect(bar()).toBe('full');
    page.mounted.shell.setExpanded(false);
    await quiet();
    expect(bar()).toBe('none');
  });

  it(
    "the trajectory's row in the Structure panel (the canvas' bar steps back while it is on screen)",
    medium,
    async () => {
      page.mounted.shell.show('structure');
      await quiet();
      expect(page.mounted.root.querySelector('[data-scene-row="scene.trajectory"]')).not.toBeNull();
      await page.check();
    },
  );

  for (const kind of ['structure', 'volume', 'volume-representation'] as const) {
    it(`the structure tab with a ${kind} picked (its inspector)`, medium, async () => {
      page.mounted.shell.show('structure');
      const scene = page.core.views.scene;
      const find = (nodes: readonly SceneNode[]): SceneNode | undefined => {
        for (const n of nodes) {
          if (n.kind === kind) return n;
          const below = find(n.children.value);
          if (below) return below;
        }
        return undefined;
      };
      const node = find(scene.roots.value);
      expect(node).toBeDefined();
      scene.pick(node!);
      await quiet();
      await page.check();
    });
  }
});

describe('the empty viewer shell (axe)', () => {
  let page: Page;
  beforeAll(async () => {
    page = await open();
  });
  afterAll(() => page.dispose());

  for (const tab of ['home', 'structure', 'help']) {
    it(`the ${tab} tab with nothing loaded`, tab === 'home' ? {} : medium, async () => {
      page.mounted.shell.show(tab);
      await quiet();
      await page.check();
    });
  }
});

describe('the viewer shell with the compact canvas toolbar (axe)', () => {
  let page: Page;
  beforeAll(async () => {
    page = await openWith({ ...DefaultPrismViewerSpec(), ui: { layout: { toolbar: 'collapsed' } } }, cbs());
  });
  afterAll(() => page.dispose());

  it('is one button, its menu open with the mode choices and a submenu', medium, async () => {
    const root = page.mounted.root;
    const button = root.querySelector<HTMLElement>('[data-stack-collapsed]');
    expect(button).not.toBeNull();
    expect(root.querySelector('[data-stack-toggle]')).toBeNull();
    await page.check();
    button!.click();
    await wait(100);
    const menu = root.querySelector('[role=menu][aria-label=Toolbar]');
    expect(menu?.textContent).not.toContain('Expand toolbar');
    expect(menu?.textContent).toContain('Interact');
    expect(menu?.textContent).toContain('What a click selects');
    await page.check();
    document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
  });
});
