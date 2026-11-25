/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer models' groups door (_spec/README.md "ViewModel", Groups): what a reveal reads of the folded groups
// (`view-folds-reveal.ts`) loads the first time a spec holds a section with a group, as render's groups chunk (core's
// `withViewKind`: the model waits for both). Why: groups are rare; a viewer that never shows one pays nothing here.

import { withViewKind } from '@molstar/seqx-core/resolve';
import { holdsGroup } from '@molstar/seqx-core/view';

/** What the chunk gives the eager models (`view-folds-reveal.ts`), once landed. */
export const GroupsView: { foldedHolding?: typeof import('./view-folds-reveal.js').foldedHolding } = {};

withViewKind(
  'section',
  () =>
    import('./view-folds-reveal.js').then((m) => {
      GroupsView.foldedHolding = m.foldedHolding;
    }),
  holdsGroup,
);
