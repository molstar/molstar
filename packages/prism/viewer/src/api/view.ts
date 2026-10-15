/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ApiArea, ApiResult } from './area.js';

// `view` (operations.md §2): the UI's state, through the navigation manager; no history step. Headless it changes the
// navigation state only (nothing shows it).

export class ViewApi extends ApiArea {
  readonly id = 'view';

  readonly open = this.operation({
    id: 'view.open',
    title: 'Open a panel',
    params: {
      tab: PD.Text('', { description: 'A tab by id: home, structure, …' }),
      card: PD.Text('', { description: 'A card over the canvas by id: export, canvas-settings, …' }),
      section: PD.Text('', { description: 'A section to bring into view in it' }),
    },
    history: false,
    async run({ tab, card, section }, ctx) {
      const nav = ctx.managers.navigation;
      if (!tab && !card) return ApiResult.fail('invalid', 'Open what?', 'Give a tab or a card.');
      if (card) {
        nav.setCard({ id: card, ...(section ? { section } : {}) });
        return ApiResult.ok(`Opened ${card}.`, { view: true });
      }
      const dock = (['right', 'left', 'top', 'bottom'] as const).find((d) =>
        ctx.managers.viewState.layout.value.docks[d].tabs.includes(tab),
      );
      if (!dock) return ApiResult.fail('not-found', `There is no '${tab}' panel.`);
      if (section) nav.setFocus({ tab, section });
      nav.setDock(dock, { active: tab, collapsed: false });
      return ApiResult.ok(`Opened ${tab}.`, { view: true });
    },
  });

  readonly close = this.operation({
    id: 'view.close',
    title: 'Close the open panel and card',
    params: {},
    history: false,
    async run(_, ctx) {
      ctx.managers.navigation.setCard(undefined);
      ctx.managers.navigation.setDock('right', { collapsed: true });
      return ApiResult.ok('Closed.', { view: true });
    },
  });
}
