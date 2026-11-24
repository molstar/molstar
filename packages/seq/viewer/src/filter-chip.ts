/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HeaderStateNames } from '@molstar/seqx-core/render';
import type { HeaderEntry } from './header-entries.js';

// The "n filtered" chip's entry (_spec/README.md "Behaviors", Headers: Entries), a module of its own: a viewer without
// the `headers` behavior lists it alone (the size script's minimal entry), and its click still clears (`Input`).

/**
 * The track's "n filtered" chip (`seq:filtered`): its count every thinned layer's rows, its click every layer's filter
 * cleared in one patch (`layer.filter`). A viewer without the `headers` behavior that lists it still clears on a click
 * (the input runs its command over the chip's layers); its other states are facts there.
 */
export const FilteredHeaderEntry: HeaderEntry = {
  name: HeaderStateNames.filtered,
  on: 'track',
  kind: 'state',
  label: 'Clear the filter',
  command: 'layer.filter',
  when: (f) => f.on === 'track' && f.layers.some((l) => l.filtered > 0),
  args: (f) => ({ addresses: f.layers.filter((l) => l.filtered > 0).map((l) => l.address) }),
  text: (f) => `${f.layers.reduce((n, l) => n + l.filtered, 0)} filtered`,
};
