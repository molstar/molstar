/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { shellIdOf } from '../extensions/overrides.js';
import { ApiArea, ApiResult } from './area.js';

// `view` (api.md §2): the UI's state, through the navigation manager; no history step. Headless it changes the
// navigation state only (nothing shows it). `open`, `close`, `expand` (the canvas fills the viewer) and `fullscreen`
// (the page's shell follows `navigation.fullscreen` with the browser's; parity.md's "expand, fullscreen" row).

/** On (true) or off (false): the state wanted, never a toggle (api.md §1; commands flip it). */
const OnParam = (description: string) => PD.Boolean(true, { description });

const DockSides = ['right', 'left', 'top', 'bottom'] as const;

export class ViewApi extends ApiArea {
  readonly id = 'view';

  /** Whether the viewport is expanded (for a command that flips it). */
  get isExpanded(): boolean {
    return this.ctx.managers.navigation.expanded.value;
  }

  /** Whether the viewer is fullscreen (for a command that flips it). */
  get isFullscreen(): boolean {
    return this.ctx.managers.navigation.fullscreen.value;
  }

  /** Whether a tab (by part id or core name) is in the layout: `view.open` finds a panel only if it is. */
  hasTab(tab: string): boolean {
    const id = shellIdOf(tab);
    const { docks } = this.ctx.managers.viewState.layout.value;
    return DockSides.some((d) => docks[d].tabs.includes(id));
  }

  readonly open = this.operation({
    id: 'view.open',
    title: 'Open a panel',
    surface: { kind: 'agent-only', reason: 'People click the panel; the agent opens it.' },
    params: {
      tab: PD.Text('', {
        description: 'A tab by part id (core.tab.home, snapshots.tab.snapshots) or a core tab: home, structure, …',
      }),
      card: PD.Text('', { description: 'A card over the canvas by part id: canvas-settings.overlay.canvas, …' }),
      section: PD.Text('', { description: 'A section to bring into view in it' }),
    },
    history: false,
    async run(values, ctx) {
      const tab = shellIdOf(values.tab);
      const card = shellIdOf(values.card);
      const { section } = values;
      const nav = ctx.managers.navigation;
      if (!tab && !card) return ApiResult.fail('invalid', 'Open what?', 'Give a tab or a card.');
      if (card) {
        nav.setCard({ id: card, ...(section ? { section } : {}) });
        return ApiResult.ok(`Opened ${card}.`, { view: true });
      }
      const dock = DockSides.find((d) => ctx.managers.viewState.layout.value.docks[d].tabs.includes(tab));
      if (!dock) return ApiResult.fail('not-found', `There is no '${tab}' panel.`);
      if (section) nav.setFocus({ tab, section });
      nav.setDock(dock, { active: tab, collapsed: false });
      return ApiResult.ok(`Opened ${tab}.`, { view: true });
    },
  });

  readonly close = this.operation({
    id: 'view.close',
    title: 'Close the open panel and card',
    surface: { kind: 'agent-only', reason: 'People click the panel; the agent opens it.' },
    params: {},
    history: false,
    async run(_, ctx) {
      const nav = ctx.managers.navigation;
      const { docks } = ctx.managers.viewState.layout.value;
      // Every dock that holds tabs and is open (`view.open` finds the tab in any of them), not the right one alone.
      const open = DockSides.filter((d) => docks[d].tabs.length > 0 && !nav.state.value.docks[d].collapsed);
      if (!open.length && !nav.state.value.card) return ApiResult.already('Nothing is open.');
      nav.setCard(undefined);
      for (const d of open) nav.setDock(d, { collapsed: true });
      return ApiResult.ok('Closed.', { view: true });
    },
  });

  readonly expand = this.operation({
    id: 'view.expand',
    title: 'Expand the viewport',
    params: { on: OnParam('Expand (true): the canvas fills the viewer; false: the panels come back') },
    history: false,
    async run({ on }, ctx) {
      const nav = ctx.managers.navigation;
      const next = on;
      if (next === nav.expanded.value) {
        return ApiResult.already(next ? 'The viewport is already expanded.' : 'The viewport is not expanded.');
      }
      // The open card goes with the chrome (ShellModel.setExpanded does the same).
      if (next) nav.setCard(undefined);
      nav.setExpanded(next);
      return ApiResult.ok(next ? 'Expanded the viewport: Escape brings the panels back.' : 'The panels are back.', {
        view: true,
      });
    },
  });

  readonly fullscreen = this.operation({
    id: 'view.fullscreen',
    title: 'Fullscreen',
    params: { on: OnParam('Enter (true) or leave (false) fullscreen') },
    history: false,
    modules: ['page'],
    // The page Mol* is embedded in may forbid it (`<iframe allow="fullscreen">`, molstar#1663).
    when: () =>
      typeof document !== 'undefined' && !document.fullscreenEnabled
        ? 'Not allowed by the page Mol* is embedded in.'
        : true,
    async run({ on }, ctx) {
      const nav = ctx.managers.navigation;
      const next = on;
      if (next === nav.fullscreen.value) {
        return ApiResult.already(next ? 'Already fullscreen.' : 'Not fullscreen.');
      }
      nav.setFullscreen(next);
      return ApiResult.ok(next ? 'Fullscreen: Escape leaves it.' : 'Left fullscreen.', { view: true });
    },
  });
}
