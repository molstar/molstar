/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a package built on prism-ui uses (the agent module, extension ports): the runtime, reactive state, the shell's
// models, commands, the components, and the shell's views a contribution sits in. Modules are also importable by path
// (`@molstar/prism-ui/dom/control-flow`) where they are `.ts`; the JSX import source is `@molstar/prism-ui/dom`.

export * from './dom/control-flow.js';
export * from './dom/focus.js';
export * from './dom/jsx-types.js';
export * from './dom/motion.js';
export * from './dom/scope.js';
export * from './dom/view.js';
export { render } from './dom/render.js';

export * from './reactive/combine-behaviors.js';
export * from './reactive/disposable-group.js';
export * from './reactive/keyed-selector.js';
export * from './reactive/keyed-state.js';
export * from './reactive/work-queue.js';
export * from './reactive/map-state.js';
export * from './reactive/observable-state.js';
export * from './reactive/readonly-state.js';
export * from './reactive/transaction.js';

export * from './model/canvas-theme.js';

export * from './managers/history.js';
export * from './managers/navigation.js';
export * from './managers/notifications.js';
export {
  DefaultDockSizes,
  type DockLayout,
  defaultViewLayout,
  MinDockSizes,
  type StoredViewLayout,
  type ViewLayout,
  type ViewLayoutSeed,
  ViewStateManager,
} from './managers/view-state.js';

export * from './shell/fullscreen.js';
export * from './shell/shell.js';
export * from './shell/theme.js';

export * from './commands/items.js';
export * from './commands/keymap.js';
export * from './commands/keys.js';
export * from './commands/manager.js';
export * from './commands/types.js';

export { ActionFooter, type ActionFooterProps } from './components/action-footer/action-footer.js';
export { Button, type ButtonProps } from './components/button/button.js';
export { CommandButton } from './components/button/command-button.js';
export { type Choice, Choices } from './components/choices/choices.js';
export { ColorPalette, ColorPicker, type ColorPickerProps, colorToHex } from './components/color/color-picker.js';
export { type ContextMenuItem, openContextMenu } from './components/context-menu/context-menu.js';
export { Legend, type LegendItem, type LegendProps, type LegendRamp } from './components/legend/legend.js';
export { Icon, type IconName, type IconProps } from './components/icon/icon.js';
export { MolstarWordmark } from './components/icon/wordmark.js';
export { GlassWordmark } from './views/shell/glass-wordmark.js';
export { Markdown, type MarkdownHost, type MarkdownLinkEvent } from './components/markdown/markdown.js';
export { ParamsForm } from './components/params/params-form.js';
export { PresetsButton, PresetsList, presetIndex } from './components/params/presets.js';
export { Popover } from './components/popover/popover.js';
export { ScrollArea, type ScrollAreaProps } from './components/scroll-area/scroll-area.js';
export {
  Floating,
  MenuPopover,
  OptionList,
  Select,
  type SelectOption,
  type SelectProps,
} from './components/select/select.js';
export { TextInput, type TextInputProps } from './components/text-input/text-input.js';
export {
  type FlatNode,
  flattenTree,
  keepRows,
  VirtualList,
  type VirtualListHandle,
} from './components/virtual-list/virtual-list.js';
export { RecycledList, RecycledListStats, type VirtualSlot } from './components/virtual-list/recycled-list.js';

export { HoverLabel } from './views/canvas/hover-label.js';
export {
  AudioRow,
  MuteButton,
  PlaybackDeck,
  SnapshotCard,
  SnapshotRow,
  TrajectoryBar,
  TrajectoryRow,
  TransitionRow,
} from './views/canvas/playback.js';
export {
  type CommandFormView,
  type CommandFormViewProps,
  type CommandMenuAssistant,
  CommandMenu,
  type CommandMenuProps,
} from './views/command-menu/command-menu.js';
export { StatusLine, StatusPill, ToolMenuRows } from './views/shell/canvas-chrome.js';
export { CanvasPlaceholder } from './views/shell/canvas-placeholder.js';
export { HelpPanel, type InputBinding, type InputBindingGroup } from './views/shell/help-panel.js';
export { LogView } from './views/shell/log-view.js';
export { SectionList } from './views/shell/section-list.js';
export { ViewHost } from './views/shell/view-host.js';
export { SettingsTab } from './views/settings/settings-tab.js';
export { ShellView, type ShellViewProps } from './views/shell/shell-view.js';
