/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, onCleanup, Show } from '@molstar/prism-ui';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { Superpose } from '../../views/structure/superpose.js';
import { createSuperpositionModel } from './model.js';

// The Superpose section of the Selection tab (parity.md §2, Superposition): the existing view over a model built on
// `superpose.*`; made when the section is shown and released with it. Offered with two or more structures.

export function superpositionUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: 'superposition.section.superpose',
        slot: 'selection',
        title: 'Superpose',
        icon: 'layers',
        order: 3,
        keywords: ['superposition', 'align', 'tm-align', 'rmsd', 'uniprot', 'sifts', 'ligands', 'atoms'],
        render: domView(() => {
          const model = createSuperpositionModel({
            api: ctx.api,
            sources: ctx.sources,
            notify: (text, kind) => ctx.managers.notifications.post(text, kind),
          });
          onCleanup(() => model.dispose());
          return (
            <Show
              when={model.available}
              fallback={() => (
                <p style="margin: 0; padding: var(--msp-space-2); color: var(--msp-color-text-subtle); font-size: var(--msp-font-size-xs)">
                  Load two or more structures to superpose them.
                </p>
              )}
            >
              {() => <Superpose model={model} />}
            </Show>
          );
        }),
      },
    ],
  };
}
