/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import { HeadersProvider, type HeadersBehavior } from '@molstar/seqx-viewer/behaviors/headers';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import type { HeaderActionEntry, HeaderEntry, HeaderGeometry, HeaderTarget } from '@molstar/seqx-viewer/header-entries';
import { headerMenuOverlay, headerPartsDom } from '@molstar/seqx-viewer/mount-header-menu';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { BehaviorProvider, SeqViewerBehavior, SeqViewerSpec } from '@molstar/seqx-viewer/viewer-spec';
import { lucideCopy } from '@molstar/ui/components/icon/_generated/lucide/copy';

// The track headers' stories, DOM-free (the stories and their headless twins call these): the default spec with the
// headers behavior's parts drawn or the host's DOM, a host's entry added, a built-in left out, and the boxes a host's
// own widget is placed at, from `viewer.headers`.

/** The host's command `story-host.copy-address`: what its entry runs, its message what the story prints. */
export const CopyCommand = 'story-host.copy-address';

/** The host's behavior that provides the command; `say` is told each run's message. */
export function copyBehavior(say: (line: string) => void = () => {}): SeqViewerBehavior {
  const provider: BehaviorProvider<undefined, { commands: readonly unknown[] }> = {
    name: 'story-host.copy',
    create: () => ({
      commands: [
        {
          id: CopyCommand,
          title: 'Copy the track’s address',
          category: 'Host',
          when: (ctx: { address?: string }) => (ctx?.address ? true : 'No track to copy'),
          run: (_args: unknown, ctx: { address: string }) => {
            const message = `Copied ${ctx.address}`;
            say(message);
            return { message };
          },
        },
      ],
    }),
  };
  return { provider, params: undefined } as unknown as SeqViewerBehavior;
}

/** The host's own header entry: "Copy address" on every track, after the built-ins' slots. */
export const CopyEntry: HeaderActionEntry = {
  kind: 'header-action',
  entry: {
    name: 'story-host:copy-address',
    on: 'track',
    icon: 'copy',
    iconMarkup: lucideCopy,
    label: 'Copy the track’s address',
    command: CopyCommand,
    order: 10,
  } satisfies HeaderEntry,
};

export interface HeadersComposition {
  /** `canvas` (default): the slots and chips drawn; `dom`: `@molstar/ui` buttons over the same boxes. */
  readonly parts?: 'canvas' | 'dom';
  /** The host's "Copy address" entry (and its command) beside the built-ins. */
  readonly host?: boolean;
  /** Built-in entries left out, by name (`seq:color`). */
  readonly without?: readonly string[];
  /** The label column's width in px: room for the slots a host adds (the default's takes two and the "⋯"). */
  readonly labelWidth?: number;
  /** Behaviors of the host's own besides the default spec's (a command its entry runs). */
  readonly behaviors?: SeqViewerSpec['behaviors'];
  /** Told the host command's runs. */
  readonly say?: (line: string) => void;
  /** What the view's own registry adds (a runtime of a format its sources name). */
  readonly registry?: SeqViewerSpec['registry'];
}

/**
 * The default spec composed as the stories show: the headers behavior with its menu and the chosen parts, the built-ins
 * minus `without`, the host's entry. Leaving an entry out is composing the registry without it.
 */
export function headersSpec(options: HeadersComposition = {}): SeqViewerSpec {
  const base = DefaultSeqViewerSpec();
  const out = new Set(options.without ?? []);
  const registry = [
    ...(base.registry ?? []).filter((e) => !(e.kind === 'header-action' && out.has(e.entry.name))),
    ...(options.registry ?? []),
    ...(options.host ? [CopyEntry] : []),
  ];
  return {
    ...base,
    registry,
    ...(options.labelWidth
      ? { config: [...(base.config ?? []), [SeqCoreConfig.render, { labelWidth: options.labelWidth }]] }
      : {}),
    behaviors: [
      ...base.behaviors.map((b) =>
        b.provider === HeadersProvider
          ? {
              provider: HeadersProvider,
              params: {
                menu: headerMenuOverlay(),
                contextClick: 'everywhere' as const,
                ...(options.parts === 'dom' ? { parts: headerPartsDom() } : {}),
              },
            }
          : b,
      ),
      ...(options.host ? [copyBehavior(options.say)] : []),
      ...(options.behaviors ?? []),
    ],
  };
}

/** The headers behavior of a viewer (it is in the default spec). */
export const headersOf = (viewer: SeqViewer): HeadersBehavior => viewer.behaviors.get('headers') as HeadersBehavior;

/** A host widget's place: its track, and the box (CSS px of the canvas's own frame) it sits in, left of the label. */
export interface WidgetPlace {
  readonly target: HeaderTarget;
  readonly x: number;
  readonly y: number;
  readonly size: number;
}

/**
 * Where a host's own widget goes: a square just left of each shown track's label box (2 px apart), centered on it,
 * from the geometry the viewer publishes (`viewer.headers`, in the canvas's frame: a host adds the canvas's own offset
 * in its page). The widget never reads the layout's internals.
 */
export function widgetPlaces(geometry: readonly HeaderGeometry[], size = 18): WidgetPlace[] {
  return geometry
    .filter((g) => g.shown && g.target.track !== undefined)
    .map((g) => ({
      target: g.target,
      x: g.label.x - size - 2,
      y: g.label.y + (g.label.height - size) / 2,
      size,
    }));
}

/** A line of text about the menu's entries for a track, as the story prints it: the items of its groups by label. */
export function menuLine(viewer: SeqViewer, target: HeaderTarget): string {
  const menu = headersOf(viewer).menuOf(target);
  return menu ? menu.groups.map((g) => g.map((i) => i.label).join(', ')).join(' | ') : 'no menu';
}
