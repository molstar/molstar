/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, type Observable, of } from 'rxjs';
import { formatChord, isMac } from '../commands/keys.js';
import type { CommandManager } from '../commands/manager.js';
import type { IconName } from '../components/icon/icon.js';
import { Show } from '../dom/control-flow.js';
import type { Child } from '../dom/jsx-types.js';
import { onCleanup } from '../dom/scope.js';
import { domView } from '../dom/view.js';
import { type CanvasInk, fixedInk } from '../model/canvas-theme.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import { observableState } from '../reactive/observable-state.js';
import { FullscreenModel } from '../shell/fullscreen.js';
import { ShellModel, type ShellTab } from '../shell/shell.js';
import { NotificationsManager } from '../managers/notifications.js';
import { ThemeModel, type UITheme, UIThemes } from '../shell/theme.js';
import { CommandMenu, type CommandMenuPrefix } from '../views/command-menu/command-menu.js';
import { CanvasPlaceholder } from '../views/shell/canvas-placeholder.js';
import { CommandManager as Commands } from '../commands/manager.js';
import { SettingsTab } from '../views/settings/settings-tab.js';
import { HelpPanel, type InputBindingGroup } from '../views/shell/help-panel.js';
import { LogView } from '../views/shell/log-view.js';
import { ShellView } from '../views/shell/shell-view.js';
import { mockSettings } from './settings.js';

// The mock frame (prism-ui's own stories): the shell and its view, the status column, the theme, fullscreen, the
// command menu, Log, Settings and Help, around a mode that fills the canvas (`mockBlankMode`: an empty one). The viewer's
// and the assistant's stories run the real viewer instead (`@molstar/prism-viewer/testing/real-viewer`).

/** What the frame gives its mode. */
export interface MockFrameContext {
  shell: ShellModel;
  status: NotificationsManager;
  loaded: boolean;
  /** The canvas background (a CSS color). */
  background?: string;
  openCommandMenu(group?: string): void;
}

/**
 * What a mode brings to the frame: its commands and the canvas parts (all optional but the commands). The parts are
 * render functions: the frame renders them in their slots, in the slot's scope (a mode activates from a click, outside
 * any render).
 */
export interface MockModeParts {
  commands: CommandManager<any>;
  prefixes?: CommandMenuPrefix[];
  canvas?: () => Child;
  overlays?: () => Child;
  playback?: () => Child;
  context?: () => Child;
  canvasLead?: () => Child;
  canvasAside?: () => Child;
  empty?: Observable<boolean>;
  onDropFiles?: (files: readonly File[]) => void;
  dropFormats?: readonly string[];
  /** How to work its canvas (the Help card's mouse and touch section). */
  bindings?: readonly InputBindingGroup[];
  /** The tab shown when the mode starts. */
  start?: string;
  /** Called when the app goes. */
  dispose?(): void;
}

export interface MockMode {
  id: string;
  title: string;
  icon: IconName;
  /** What it is for, in a line (the mode chooser). */
  description?: string;
  activate(ctx: MockFrameContext): MockModeParts;
}

/** The frame alone, with an empty canvas: for stories about the frame or what is added to it. */
export function mockBlankMode(): MockMode {
  return {
    id: 'blank',
    title: 'Blank',
    icon: 'square',
    activate: ({ background }) => ({
      commands: new Commands<object>({ context: observableState({}).asReadonly, formats: { get: () => undefined } }),
      canvas: () => <CanvasPlaceholder empty={false} background={background} />,
      empty: of(false),
    }),
  };
}

export interface ShellStoryProps {
  loaded: boolean;
  background?: string;
  /** The canvas ink as a state that can change (instead of the fixed one of `background`). */
  ink?: ReadonlyState<CanvasInk>;
  /** Opens the Export card at a section (a mode registers the card). */
  openExportAt?: string;
  /** Opens a card (Animation, Export, …). */
  openCard?: string;
  /** Shows a tab (Help, Settings, …). */
  openTab?: string;
  /** Starts with no panel open: the rail and the canvas tools only (a tab opens as usual when its rail button is clicked). */
  noPanel?: boolean;
  /** The right panel pinned beside the canvas. */
  pinned?: boolean;
  /** A compact rail (`ShellOptions.compactRail`): the rail column hidden behind the canvas corner's button. */
  compact?: boolean;
  /** The canvas toolbar is one button with a menu (`ShellOptions.toolbar`). */
  toolbarCollapsed?: boolean;
  /** Opens with the command menu open. */
  commandMenu?: boolean;
  /** Tabs of the story's own (a React panel), with the shell and the status column to work with. */
  tabs?: (ctx: { shell: ShellModel; status: NotificationsManager }) => readonly ShellTab[];
}

/** The mock frame running a mode (`mockBlankMode`). */
export function ShellStory(props: ShellStoryProps & { mode: MockMode }) {
  const { mode } = props;
  const shell = new ShellModel({
    compactRail: !!props.compact,
    ...(props.toolbarCollapsed ? { toolbar: 'collapsed' as const } : {}),
    docks: { right: { pinned: !!props.pinned, collapsed: !!props.noPanel } },
  });
  const status = new NotificationsManager();
  const commandMenuOpen = observableState(false);
  // The group the command menu opens in (Preset, Styles); none for ⌘K.
  const commandMenuGroup = observableState<string | undefined>(undefined);
  const openCommandMenu = (group?: string) => {
    commandMenuGroup.set(group);
    commandMenuOpen.set(true);
  };
  // The plugin root goes fullscreen (in the stories, the story page).
  const fullscreen = new FullscreenModel(() => document.documentElement);
  // The UI theme, starting from Storybook's (its toolbar sets the same attribute on the page).
  const theme = new ThemeModel((document.documentElement.dataset.mspTheme as UITheme | undefined) ?? 'light');
  theme.applyTo(document.documentElement);
  const settings = mockSettings();

  // The frame's own tabs and buttons: in every mode.
  shell.registerTabs([
    {
      id: 'log',
      title: 'Log',
      icon: 'scroll-text',
      order: 2.5,
      // Not a working panel: with Settings at the rail's foot.
      place: 'foot',
      render: domView(() => <LogView status={status} />),
    },
    {
      id: 'settings',
      title: 'Settings',
      icon: 'settings',
      order: 3,
      place: 'foot',
      render: domView(() => <SettingsTab model={settings} />),
    },
  ]);

  if (props.tabs) shell.registerTabs(props.tabs({ shell, status }));
  const parts = mode.activate({
    shell,
    status,
    loaded: props.loaded,
    background: props.background,
    openCommandMenu,
  });
  // The canvas tools act on what is shown: the mode says whether anything is.
  const empty = parts.empty ?? of(!props.loaded);
  const hasData = empty.pipe(map((e) => !e));
  shell.registerToolbar([
    // The UI theme, by Settings at the rail's foot: a click goes to the next one; the icon is the current one.
    {
      id: 'theme',
      // The first to go under More on a short rail.
      priority: 100,
      icon: 'sun',
      title: 'Theme',
      face: theme.theme.pipe(
        map((t) => {
          const current = UIThemes.find((u) => u.id === t)!;
          return { icon: current.icon, title: `Theme: ${current.title}` };
        }),
      ),
      run: () => theme.next(),
    },
    {
      id: 'palette',
      // The main way in: the last to go.
      priority: -1,
      icon: 'command',
      title: 'Command menu',
      caption: formatChord({ key: 'k', mod: true }),
      run: () => openCommandMenu(),
    },
    // Expand by context (FullscreenModel.mode): on a page, Expand fills it and Fullscreen is beside it; in an iframe
    // that allows it, the button is Fullscreen (expanding would only fill the frame); in one that does not, Expand fills
    // the frame and Fullscreen shows why it is not possible.
    fullscreen.mode === 'fullscreen'
      ? {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'fullscreen' as const,
          title: 'Fullscreen',
          pressed: fullscreen.active,
          run: () => fullscreen.toggle(),
        }
      : {
          id: 'expand',
          place: 'canvas',
          group: 'canvas',
          visible: hasData,
          icon: 'maximize-2' as const,
          title: 'Expand',
          run: () => status.post('Expanded'),
          menu: [
            {
              id: 'fullscreen',
              label: 'Fullscreen',
              icon: 'fullscreen' as const,
              hint: 'The whole screen, not just the page',
              disabled: fullscreen.unavailableReason,
              pressed: fullscreen.active,
              run: () => fullscreen.toggle(),
            },
          ],
        },
  ]);

  // Help: how to work the canvas (the mode's bindings) and the keys (the command registry as it is). A tab at the
  // rail's foot, not a card: it stays open while you try what it says on the canvas (a card closes on a click outside).
  // Also a command on ?, and the command menu's ?.
  shell.registerTabs([
    {
      id: 'help',
      title: 'Help',
      icon: 'circle-help',
      order: 2.8,
      place: 'foot',
      render: domView(() => <HelpPanel commands={parts.commands} bindings={parts.bindings ?? []} />),
    },
  ]);
  parts.commands.register([
    {
      id: 'help.show',
      title: 'Keyboard and mouse',
      alias: ['help', 'shortcuts', 'keys', 'bindings', 'controls'],
      category: 'Help',
      icon: 'circle-help',
      shortcut: { key: '?' },
      run: () => {
        shell.show('help');
        return { message: '' };
      },
    },
  ]);

  // A compact rail's column: shown and hidden by the corner's button, or by this (bindable) command.
  parts.commands.register([
    {
      id: 'shell.toggle-rail',
      title: 'Show or hide the panel rail',
      alias: ['rail', 'sidebar', 'panels', 'compact'],
      category: 'View',
      icon: 'prism',
      run: () => {
        shell.toggleRail();
        return { message: '' };
      },
    },
  ]);

  // What commands did goes to the status line.
  const outcomes = parts.commands.outcomes.subscribe((o) => {
    if (o.kind === 'done' && o.result.message) status.post(o.result.message);
    else if (o.kind === 'failed' || o.kind === 'blocked') status.post(o.message, 'warning');
  });

  // ⌘K / ⇧⌘P open the command menu; the mode's commands' own shortcuts run them.
  const onKey = (e: KeyboardEvent) => {
    const mod = isMac() ? e.metaKey : e.ctrlKey;
    if (mod && (e.key.toLowerCase() === 'k' || (e.shiftKey && e.key.toLowerCase() === 'p'))) {
      e.preventDefault();
      openCommandMenu();
    } else if (!commandMenuOpen.value && parts.commands.handleKey(e, isMac())) e.preventDefault();
  };
  window.addEventListener('keydown', onKey);
  if (parts.start) shell.show(parts.start);
  shell.mount();
  status.mount();
  if (props.openExportAt) shell.openOverlay('export', undefined, props.openExportAt);
  if (props.openCard) shell.openOverlay(props.openCard);
  if (props.openTab) shell.show(props.openTab);
  if (props.commandMenu) commandMenuOpen.set(true);
  onCleanup(() => {
    window.removeEventListener('keydown', onKey);
    outcomes.unsubscribe();
    parts.dispose?.();
    shell.unmount();
    status.unmount();
  });

  const background = props.background ?? '#ffffff';
  return (
    <ShellView
      shell={shell}
      status={status}
      canvas={parts.canvas?.() ?? null}
      overlays={parts.overlays?.() ?? null}
      playback={parts.playback?.() ?? null}
      ink={props.ink ?? fixedInk(background)}
      context={parts.context}
      commandMenu={
        <Show when={commandMenuOpen}>
          {() => (
            <CommandMenu
              manager={parts.commands}
              prefixes={parts.prefixes}
              group={commandMenuGroup.value}
              onClose={() => {
                commandMenuOpen.set(false);
                commandMenuGroup.set(undefined);
              }}
            />
          )}
        </Show>
      }
      canvasLead={parts.canvasLead?.() ?? null}
      canvasAside={parts.canvasAside?.() ?? null}
      onOpenLog={() => shell.show('log')}
      empty={empty}
      onDropFiles={parts.onDropFiles}
      dropFormats={parts.dropFormats}
    />
  );
}
