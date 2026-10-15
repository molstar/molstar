/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CommandDefinition, domView, type ShellModel, type ShellOverlay, type ShellTab } from '@molstar/prism-ui';
import type { CanvasToolbarItem } from '@molstar/prism-ui/shell/shell';
import { map } from 'rxjs';
import type { ViewerApi } from '../api/area.js';
import type { ViewerCommandContext } from '../api/managers.js';
import type { SequenceSource } from '../sources/sequence.js';
import { LegendOverlay, sequenceCommands, sequenceMenu } from './commands.js';
import { LegendCard } from './legend-card.js';
import type { SequenceViewers } from './model.js';
import type { MarksSource } from '../sources/marks.js';
import { SequenceControl } from './sequence-control.js';

// The sequence's parts in the page's shell (prism-sequence.md §2.2): its commands, the legend card, the band (a tab of
// the top dock, there while a structure is shown; the first one shows it) and its canvas button. In the sequence's
// chunk (`runtime.ts`): the shell mounts them once the runtime is in, as a structure is first shown.

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
  readonly marks?: MarksSource;
}

/** Mounts the sequence's parts; returns their removal. */
export function sequenceShellUI(d: SequenceShellDeps): () => void {
  const { shell, own } = d;
  const viewer = d.viewers.get();
  const undo: (() => void)[] = [];
  undo.push(d.registerCommands(sequenceCommands(viewer, shell, d.person)));
  undo.push(
    shell.registerOverlays(
      own<ShellOverlay>('overlay', [
        { id: LegendOverlay, title: 'Sequence colors', render: domView(() => <LegendCard viewer={viewer} />) },
      ]),
    ),
  );
  let dropTab: (() => void) | undefined;
  let shownBefore = false;
  const bandSub = viewer.structure.subscribe((shown) => {
    if (!shown) {
      dropTab?.();
      dropTab = undefined;
      return;
    }
    if (dropTab) return;
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
            {...(d.marks ? { marks: d.marks } : {})}
          />
        )),
      },
    ]);
    if (!tabs.length) return;
    dropTab = shell.registerTabs(tabs);
    if (!shownBefore) shell.show('sequence');
    shownBefore = true;
  });
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
          visible: viewer.structure.pipe(map((s) => s !== undefined)),
          menu: sequenceMenu(viewer, shell, d.openCommand),
        },
      ]),
    ),
  );
  return () => {
    for (const u of undo.reverse()) u();
  };
}
