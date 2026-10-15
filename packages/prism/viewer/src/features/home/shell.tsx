/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, mapState, MolstarWordmark, SectionList } from '@molstar/prism-ui';
import { createAnnotationOffer, type HomeModel } from '../../model/home.js';
import type { BuiltInShellContext } from '../../spec.js';
import { HomeTab } from '../../views/home/home-tab.js';

/** Home (the core's id, `home`), at the head of the rail with the Mol* mark. */
export function homeShell({ shell, core, person, own, openCommandMenu, openFiles }: BuiltInShellContext): () => void {
  const { sources } = core;
  const home: HomeModel = {
    actions: core.views.homeActions,
    // After a load, what a structure carries beyond the built-ins, offered in one click (parity.md §6a).
    annotationOffer: createAnnotationOffer(person, sources.structures),
    // The records kept in earlier sessions: one opens in History, read-only (where History is shown).
    ...(core.options.ui.showHistory
      ? {
          keptRecords: {
            list: core.keptRecords.list,
            watch: () => core.keptRecords.watch(),
            forget: (id: string) => core.keptRecords.forget(id),
            open: (id: string) =>
              void core.keptRecords.open(id).then((ok) => {
                if (ok) shell.show('history');
              }),
          },
        }
      : {}),
    formats: sources.formats,
    openFiles,
    session: {
      isEmpty: sources.empty,
      summary: mapState(sources.structures, summaryOf, Object.is),
      clear: () => void person.session.clear(),
    },
  };
  return shell.registerTabs(
    own('tab', [
      {
        id: 'home',
        compactBar: true,
        title: 'Home',
        heading: () => <MolstarWordmark live />,
        icon: 'prism',
        order: 0,
        place: 'head',
        render: domView(({ ctx }) => (
          <HomeTab
            home={home}
            onOpenCommandMenu={() => openCommandMenu()}
            // Features' import sections (the loaders' database and URL forms) after the registry's actions.
            more={() => <SectionList shell={ctx.shell} slot="home.import" bare collapsed />}
          />
        )),
      },
    ]),
  );
}

/** What is loaded, in a line ("1CBS · 1,213 atoms", "2 structures · …"). */
function summaryOf(list: readonly { label: string; atoms: number }[]): string {
  if (!list.length) return 'Nothing loaded';
  if (list.length === 1) return `${list[0].label} · ${list[0].atoms.toLocaleString('en')} atoms`;
  return `${list.length} structures · ${list.reduce((n, s) => n + s.atoms, 0).toLocaleString('en')} atoms`;
}
