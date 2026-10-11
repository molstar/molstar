/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, type Observable, of, switchMap } from 'rxjs';
import { formatChord, isMac } from '../commands/keys.js';
import type { CommandManager } from '../commands/manager.js';
import type { CommandDefinition } from '../commands/types.js';
import type { IconName } from '../components/icon/icon.js';
import { Dynamic, Show } from '../dom/control-flow.js';
import type { Child } from '../dom/jsx-types.js';
import { onCleanup } from '../dom/scope.js';
import { domView } from '../dom/view.js';
import { parseColor, themeForColor } from '../model/canvas-theme.js';
import { observableState } from '../reactive/observable-state.js';
import { FullscreenModel } from '../shell/fullscreen.js';
import { ShellModel, type ShellTab } from '../shell/shell.js';
import { StatusModel } from '../shell/status.js';
import { ThemeModel, type UITheme, UIThemes } from '../shell/theme.js';
import { CommandMenu, type CommandMenuAssistant, type CommandMenuPrefix } from '../views/command-menu/command-menu.js';
import { CanvasPlaceholder } from '../views/shell/canvas-placeholder.js';
import { type ModeChoice, ModeChooser } from '../views/shell/mode-chooser.js';
import { CommandManager as Commands } from '../commands/manager.js';
import { SettingsTab } from '../views/settings/settings-tab.js';
import { LogView } from '../views/shell/log-view.js';
import { ShellView } from '../views/shell/shell-view.js';
import { mockSettings } from './settings.js';

// The mock shell, in two layers: the frame (the shell and its view, the status column, the theme, fullscreen, the
// command menu, Log and Settings, and the add-ons such as the agent) and a mode, what the app is for at the moment:
// Mol*'s own UI (@molstar/prism-viewer) or a story mode (@molstar/prism-stories); an app is composed of the modes it
// includes. A mode brings its tabs, overlays, sections,
// toolbar items, commands and canvas parts; the frame swaps it at runtime, undoing all it registered
// (design/modes.md).

/** What the frame gives a mode while it is active. Its `shell` records what the mode registers, to undo it. */
export interface MockFrameContext {
  shell: ShellModel;
  status: StatusModel;
  loaded: boolean;
  /** The canvas background (a CSS color). */
  background?: string;
  openCommandMenu(group?: string): void;
  /** The app's modes (to offer the choice, as Home does while nothing is loaded), the active one, and switching. */
  modes: readonly ModeChoice[];
  activeMode: Observable<string>;
  switchMode(id: string): void;
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
  /** The tab shown when the mode starts. */
  start?: string;
  /** Called when the mode leaves (its registrations are undone by the frame). */
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

/**
 * What a module adds to the mock shell, the way an app would contribute it: tabs (one with an existing id replaces it,
 * as History replaces Log), the command menu's assistant, the canvas's agent slot, and the person's finished commands.
 */
export interface MockShellAddOnParts {
  tabs?: readonly ShellTab[];
  assistant?: CommandMenuAssistant;
  agent?: Child;
  commandDone?(message: string, commandId: string): void;
}

export type MockShellAddOn = (ctx: { shell: ShellModel; status: StatusModel; loaded: boolean }) => MockShellAddOnParts;

function withAddOnTabs(tabs: ShellTab[], extra: readonly ShellTab[] | undefined): ShellTab[] {
  if (!extra?.length) return tabs;
  const replaced = new Set(extra.map((t) => t.id));
  return [...tabs.filter((t) => !replaced.has(t.id)), ...extra];
}

/** The shell as a mode sees it: every `register…` call is recorded, so leaving the mode undoes them all. */
function recordingShell(shell: ShellModel): { shell: ShellModel; undo(): void } {
  const undos: (() => void)[] = [];
  const proxy = new Proxy(shell, {
    get(target, key) {
      const value = Reflect.get(target, key, target);
      if (typeof value !== 'function') return value;
      if (typeof key === 'string' && key.startsWith('register')) {
        return (...args: unknown[]) => {
          const undo = value.apply(target, args) as () => void;
          undos.push(undo);
          return undo;
        };
      }
      return value.bind(target);
    },
  });
  return {
    shell: proxy,
    undo: () => {
      for (let i = undos.length - 1; i >= 0; i--) undos[i]();
      undos.length = 0;
    },
  };
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
  /** Opens the Export card at a section (a mode registers the card). */
  openExportAt?: string;
  /** The right panel pinned beside the canvas. */
  pinned?: boolean;
  /** Opens with the command menu open. */
  commandMenu?: boolean;
  /** The first open: the app asks which mode to start in (a card over the canvas) before any is used. */
  chooseMode?: boolean;
  /** What an add-on brings (the agent module's mock: `@molstar/prism-agent/testing`). */
  addOn?: MockShellAddOn;
}

/**
 * The mock frame with the modes the app is composed of, the first active; with more, the rail and the command menu
 * switch between them. Mol*'s own UI is `mockViewerMode` (@molstar/prism-viewer), story mode `mockStoryMode`
 * (@molstar/prism-stories).
 */
export function ShellStory(props: ShellStoryProps & { modes: readonly MockMode[] }) {
  const modes = props.modes;
  const shell = new ShellModel({ docks: { right: { pinned: !!props.pinned } } });
  const status = new StatusModel();
  const commandMenuOpen = observableState(false);
  // The first open asks for a mode (with more than one).
  const choosing = observableState(!!props.chooseMode && props.modes.length > 1);
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
  // What an add-on brings (the agent module's tabs, its Ask port and turn card): only when the story includes one.
  const addOn = props.addOn?.({ shell, status, loaded: props.loaded });
  const settings = mockSettings();

  // The frame's own tabs and buttons: in every mode.
  shell.registerTabs(
    withAddOnTabs(
      [
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
      ],
      addOn?.tabs,
    ),
  );

  const active = observableState<MockMode>(modes[0]);
  const parts = observableState<MockModeParts | undefined>(undefined);
  // The canvas tools act on what is shown: the active mode says whether anything is.
  const hasData = parts.asReadonly.pipe(
    switchMap((p) => p?.empty ?? of(!props.loaded)),
    map((empty) => !empty),
  );
  shell.registerToolbar([
    // With more than one mode, the rail switches to the next (its face is where it goes).
    ...(modes.length > 1
      ? [
          {
            id: 'mode',
            icon: modes[1].icon,
            title: `Switch to ${modes[1].title}`,
            face: active.asReadonly.pipe(
              map((m) => {
                const next = nextMode(m);
                return { icon: next.icon, title: `Switch to ${next.title}` };
              }),
            ),
            run: () => switchTo(nextMode(active.value)),
          },
        ]
      : []),
    // The UI theme, by Settings at the rail's foot: a click goes to the next one; the icon is the current one.
    {
      id: 'theme',
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

  const nextMode = (m: MockMode) => modes[(modes.indexOf(m) + 1) % modes.length];
  // The active mode: activated with a recording shell; leaving it undoes its registrations and command wiring.
  let leave: (() => void) | undefined;
  const switchTo = (mode: MockMode) => {
    leave?.();
    commandMenuOpen.set(false);
    const recording = recordingShell(shell);
    const p = mode.activate({
      shell: recording.shell,
      status,
      loaded: props.loaded,
      background: props.background,
      openCommandMenu,
      modes,
      activeMode: active.asReadonly.pipe(map((m) => m.id)),
      switchMode: (id) => {
        const next = modes.find((m) => m.id === id);
        if (next && next !== active.value) queueMicrotask(() => switchTo(next));
      },
    });
    // Switching modes is a command in each of them.
    const switching: CommandDefinition<any>[] = modes
      .filter((m) => m !== mode)
      .map((m) => ({
        id: `mode.${m.id}`,
        title: `Switch to ${m.title}`,
        category: 'View',
        icon: m.icon,
        run: () => {
          queueMicrotask(() => switchTo(m));
          return { message: '' };
        },
      }));
    const unregister = switching.length ? p.commands.register(switching) : undefined;
    // What commands did goes to the status line (and, with the agent, to History).
    const outcomes = p.commands.outcomes.subscribe((o) => {
      if (o.kind === 'done' && o.result.message) status.post(o.result.message);
      if (o.kind === 'done') addOn?.commandDone?.(o.result.message ?? o.command.title, o.command.id);
      else if (o.kind === 'failed' || o.kind === 'blocked') status.post(o.message, 'warning');
    });
    leave = () => {
      outcomes.unsubscribe();
      unregister?.();
      recording.undo();
      p.dispose?.();
    };
    active.set(mode);
    parts.set(p);
    // While the app asks for a mode (the first open), nothing opens beside the question.
    if (p.start && !choosing.value) shell.show(p.start);
  };

  // ⌘K / ⇧⌘P open the command menu; the active mode's commands' own shortcuts run them.
  const onKey = (e: KeyboardEvent) => {
    const mod = isMac() ? e.metaKey : e.ctrlKey;
    if (mod && (e.key.toLowerCase() === 'k' || (e.shiftKey && e.key.toLowerCase() === 'p'))) {
      e.preventDefault();
      openCommandMenu();
    } else if (!commandMenuOpen.value && parts.value?.commands.handleKey(e, isMac())) e.preventDefault();
  };
  window.addEventListener('keydown', onKey);
  switchTo(modes[0]);
  // Asking for a mode: the panel stays closed until one is chosen (its first tab would be active).
  if (choosing.value && !shell.dock('right').state.value.collapsed) shell.toggle('right');
  shell.mount();
  status.mount();
  if (props.openExportAt) shell.openOverlay('export', undefined, props.openExportAt);
  if (props.commandMenu) commandMenuOpen.set(true);
  onCleanup(() => {
    window.removeEventListener('keydown', onKey);
    leave?.();
    shell.unmount();
    status.unmount();
  });

  // A canvas part of the active mode, swapped with it.
  const part = (pick: (p: MockModeParts) => Child) => (
    <Dynamic component={parts.asReadonly.pipe(map((p) => (p ? () => pick(p) : undefined)))} />
  );
  const background = props.background ?? '#ffffff';
  return (
    <ShellView
      shell={shell}
      status={status}
      canvas={part((p) => p.canvas?.() ?? null)}
      overlays={
        <>
          {part((p) => p.overlays?.() ?? null)}
          <Show when={choosing}>
            {() => (
              <ModeChooser
                variant="welcome"
                modes={modes}
                onChoose={(id) => {
                  choosing.set(false);
                  const next = modes.find((m) => m.id === id);
                  if (next && next !== active.value) switchTo(next);
                  else if (parts.value?.start) shell.show(parts.value.start);
                }}
              />
            )}
          </Show>
        </>
      }
      playback={part((p) => p.playback?.() ?? null)}
      canvasColor={background}
      canvasTheme={themeForColor(parseColor(background) ?? 0xffffff)}
      context={() => part((p) => p.context?.() ?? null)}
      commandMenu={
        <Show when={commandMenuOpen}>
          {() => (
            <CommandMenu
              manager={parts.value!.commands}
              prefixes={parts.value!.prefixes}
              group={commandMenuGroup.value}
              assistant={addOn?.assistant}
              onClose={() => {
                commandMenuOpen.set(false);
                commandMenuGroup.set(undefined);
              }}
            />
          )}
        </Show>
      }
      agent={addOn?.agent}
      canvasLead={part((p) => p.canvasLead?.() ?? null)}
      canvasAside={part((p) => p.canvasAside?.() ?? null)}
      onOpenLog={() => shell.show('log')}
      empty={parts.asReadonly.pipe(switchMap((p) => p?.empty ?? of(!props.loaded)))}
      onDropFiles={(files) => parts.value?.onDropFiles?.(files)}
      dropFormats={parts.value?.dropFormats}
    />
  );
}
