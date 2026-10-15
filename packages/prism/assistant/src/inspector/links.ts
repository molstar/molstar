/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import type { InspectorFocus } from '../inspector-access.js';
import type { InspectLink } from '../views/assistant-tab.js';
import { type InspectorModel, toolCounts } from './model.js';

/**
 * The Assistant tab's Inspect rows (design/agent.md §6): Trace with its turns, Tools with its tools and skills, Limits
 * with what they bound; each opens the Inspector there (`InspectLink.id` is an `InspectorFocus`).
 */
export function inspectLinks(model: InspectorModel, limits = 'turn budget, context window'): InspectLink[] {
  const plural = (n: number, one: string) => `${n} ${one}${n === 1 ? '' : 's'}`;
  const link = (id: InspectorFocus, label: string, detail: InspectLink['detail']): InspectLink => ({
    id,
    label,
    detail,
  });
  return [
    link(
      'trace',
      'Trace',
      model.turns.pipe(
        map((t) => (t.length ? plural(t.length, 'turn') : 'nothing yet')),
        distinctUntilChanged(),
      ),
    ),
    link(
      'tools',
      'Tools',
      model.tools.pipe(
        map(toolCounts),
        map((c) => `${plural(c.tools, 'tool')}, ${plural(c.skills, 'skill')}`),
        distinctUntilChanged(),
      ),
    ),
    link('limits', 'Limits', limits),
  ];
}
