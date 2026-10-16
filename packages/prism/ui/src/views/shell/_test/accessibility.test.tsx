/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { domView } from '../../../dom/view.js';
import { render } from '../../../dom/render.js';
import { NotificationsManager } from '../../../managers/notifications.js';
import { ShellModel, type ShellViewContext } from '../../../shell/shell.js';
import { ShellView } from '../shell-view.js';
import { expectNoAxeViolations } from '../../../testing/axe.js';
import { mockBlankMode, ShellStory, type ShellStoryProps } from '../../../testing/shell.js';

// Automated checks (axe-core, jsdom) of the shell's views, mounted as the mock frame mounts them (accessibility.md,
// Automated checks): the rail and docks, the panels (Log, Settings, Help), the command menu, a card.

const frame = () => new Promise((r) => setTimeout(r, 32));

function mount(props: Partial<ShellStoryProps> = {}) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(() => <ShellStory mode={mockBlankMode()} loaded {...props} />, container);
  return {
    container,
    done: () => {
      dispose();
      container.remove();
    },
  };
}

describe('the shell has no accessibility violations (axe)', () => {
  beforeEach(() => {
    vi.stubGlobal(
      'ResizeObserver',
      class {
        observe() {}
        disconnect() {}
        unobserve() {}
      },
    );
  });
  afterEach(() => {
    document.body.innerHTML = '';
    vi.unstubAllGlobals();
  });

  it('the rail, the canvas chrome and the docks', async () => {
    const t = mount({ noPanel: true });
    await frame();
    // Something is there to check.
    expect(t.container.querySelectorAll('[role=tab]').length).toBeGreaterThan(1);
    await expectNoAxeViolations(t.container);
    t.done();
  });

  it('the empty canvas', async () => {
    const t = mount({ loaded: false, noPanel: true });
    await frame();
    await expectNoAxeViolations(t.container);
    t.done();
  });

  for (const tab of ['log', 'settings', 'help']) {
    it(`the ${tab} panel`, async () => {
      const t = mount({ openTab: tab });
      await frame();
      await expectNoAxeViolations(t.container);
      t.done();
    });
  }

  it('the right panel pinned beside the canvas', async () => {
    const t = mount({ openTab: 'settings', pinned: true });
    await frame();
    await expectNoAxeViolations(t.container);
    t.done();
  });

  it('a compact rail', async () => {
    const t = mount({ compact: true, openTab: 'settings' });
    await frame();
    await expectNoAxeViolations(t.container);
    t.done();
  });

  it('the command menu open', async () => {
    const t = mount({ commandMenu: true });
    await frame();
    await expectNoAxeViolations(t.container);
    t.done();
  });

  it('an overlay card (a dialog) open over a panel, with the status line saying something', async () => {
    const shell = new ShellModel();
    shell.registerTabs([
      { id: 'home', title: 'Home', icon: 'house', render: domView<ShellViewContext>(() => <p>Home</p>) },
      {
        id: 'structure',
        title: 'Structure',
        icon: 'box',
        render: domView<ShellViewContext>(() => <input aria-label="Filter" />),
      },
    ]);
    shell.registerOverlays([
      { id: 'card', title: 'Card', render: domView<ShellViewContext>(() => <input aria-label="Width" />) },
    ]);
    const status = new NotificationsManager();
    const container = document.createElement('div');
    document.body.append(container);
    const dispose = render(
      () => (
        <ShellView
          shell={shell}
          status={status}
          canvas={<canvas />}
          context={() => <button type="button">Clear selection</button>}
        />
      ),
      container,
    );
    shell.show('structure');
    await frame();
    await expectNoAxeViolations(container);
    shell.openOverlay('card');
    status.post('Selected 3 atoms');
    await frame();
    expect(container.querySelector('[role=dialog]')).not.toBeNull();
    await expectNoAxeViolations(container);
    dispose();
    container.remove();
  });
});
