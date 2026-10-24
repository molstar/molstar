/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type CommandDefinition,
  type CommandGroup,
  domView,
  type ShellModel,
  type ShellOverlay,
  type ShellTab,
} from '@molstar/prism-ui';
import type { CanvasToolbarItem } from '@molstar/prism-ui/shell/shell';
import { combineLatest } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import type { ViewerCommandContext } from '../api/managers.js';
import type { SequenceSource } from '../sources/sequence.js';
import { AddTrackGroup, addTrackCommands, LegendOverlay, sequenceCommands, sequenceMenu } from './commands.js';
import { LegendCard } from './legend-card.js';
import type { SequenceViewers } from './model.js';
import type { MarksView } from '../sources/marks.js';
import { SequenceControl } from './sequence-control.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// The sequence's parts in the page's shell (sequence.md §8): its commands, the legend card, the band (a tab of
// the top dock, there from the first document in until the host shows no structure; the first one shows it) and its
// canvas button. In the sequence's chunk (`runtime.ts`): the shell mounts them once the runtime is in, as a structure is
// first shown.

export interface SequenceShellDeps {
  readonly viewers: SequenceViewers;
  readonly source: SequenceSource;
  readonly shell: ShellModel;
  /** The API as the person. */
  readonly person: ViewerApi;
  readonly registerCommands: (defs: readonly CommandDefinition<ViewerCommandContext>[]) => () => void;
  /** The spec's overrides of the core's own parts (`core.tab.sequence`, …). */
  readonly own: <T extends { id: string }>(kind: string, parts: readonly T[]) => T[];
  readonly openCommand: (id: string) => void;
  /** The band's hover label, as lines (the shell merges it with the canvas's). */
  readonly onHover: (lines: readonly string[]) => void;
  /** A view's marks and hover (marks.md), shown in the sequence too. */
  readonly marks?: MarksView;
  /** Opens the structure tree at a structure's Sequence node (sequence.md §8.1); none: the page has no tree. */
  readonly reveal?: (structure: string) => void;
  /** Registers command groups (Add track's). */
  readonly registerGroups?: (groups: readonly CommandGroup[]) => () => void;
  /** Opens the command menu at a group. */
  readonly openCommandMenu?: (group: string) => void;
}

/** Mounts the sequence's parts; returns their removal. */
export function sequenceShellUI(d: SequenceShellDeps): () => void {
  const { shell, own } = d;
  const viewer = d.viewers.get();
  const undo: (() => void)[] = [];
  // Edit in structure tree: the Sequence node of the structure shown.
  const reveal = d.reveal
    ? () => {
        const shown = viewer.structure.value;
        if (shown) d.reveal!(shown.ref);
      }
    : undefined;
  undo.push(d.registerCommands(sequenceCommands(viewer, shell, d.person, reveal)));
  // Add track (sequence.md §8.1): a command a provider, in a group of their own the Sequence node's inspector opens
  // the command menu at; made again when providers come or go.
  if (d.registerGroups) undo.push(d.registerGroups([AddTrackGroup]));
  let dropAdd: (() => void) | undefined;
  const providersSub = viewer.providers.subscribe((list) => {
    dropAdd?.();
    dropAdd = d.registerCommands(addTrackCommands(viewer, list, d.person));
  });
  undo.push(() => {
    providersSub.unsubscribe();
    dropAdd?.();
  });
  // The band shown or not, for the structure tree's Sequence eye (sequence.md §8.1): the top dock, as No sequence.
  undo.push(
    d.source.setBand({
      shown: derive(shell.dock('top').active, (t) => !!t),
      show: (on) => {
        if (on) shell.show('sequence');
        else if (shell.dock('top').active.value) shell.toggle('top');
      },
      ...(d.openCommandMenu ? { openAddTrack: () => d.openCommandMenu!(AddTrackGroup.id) } : {}),
    }),
  );
  undo.push(
    shell.registerOverlays(
      own<ShellOverlay>('overlay', [
        { id: LegendOverlay, title: 'Sequence colors', render: domView(() => <LegendCard viewer={viewer} />) },
      ]),
    ),
  );
  let dropTab: (() => void) | undefined;
  let shownBefore = false;
  // Another structure's document being made (`viewer.structure` undefined meanwhile) keeps the band: its surface, link
  // and the dock's size are not torn down and made again on each switch.
  // No chain shown (the chain filter's None): nothing to draw, so no band until a chain is picked again.
  let emptied = false;
  const bandSub = combineLatest([viewer.structure, d.source.structure, viewer.empty]).subscribe(
    ([shown, host, empty]) => {
      if (!host || empty) {
        if (empty && dropTab) emptied = true;
        dropTab?.();
        dropTab = undefined;
        return;
      }
      if (!shown || dropTab) return;
      const tabs = own<ShellTab>('tab', [
        {
          id: 'sequence',
          title: 'Sequence',
          icon: 'sequence',
          dock: 'top',
          frame: 'none',
          render: domView(({ ctx }) => (
            <SequenceControl
              viewer={viewer}
              view3d={(s) => d.source.view3d(s, () => void d.person.camera.reset())}
              shell={ctx.shell}
              dock={ctx.dock ?? 'top'}
              onHover={d.onHover}
              api={d.person}
              {...(d.marks ? { marks: d.marks } : {})}
            />
          )),
        },
      ]);
      if (!tabs.length) return;
      dropTab = shell.registerTabs(tabs);
      if (!shownBefore || emptied) shell.show('sequence');
      emptied = false;
      shownBefore = true;
    },
  );
  undo.push(() => bandSub.unsubscribe());
  undo.push(() => dropTab?.());
  undo.push(
    shell.registerToolbar(
      own<CanvasToolbarItem>('toolbar', [
        {
          id: 'sequence',
          place: 'canvas',
          group: 'view',
          icon: 'sequence',
          title: 'Sequence',
          visible: derive(viewer.structure, (s) => s !== undefined),
          menu: sequenceMenu(viewer, shell, d.openCommand, !!reveal),
        },
      ]),
    ),
  );
  return () => {
    for (const u of undo.reverse()) u();
  };
}
