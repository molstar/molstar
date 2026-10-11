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
export * from './dom/jsx-types.js';
export * from './dom/motion.js';
export * from './dom/scope.js';
export * from './dom/view.js';
export { render } from './dom/render.js';

export * from './reactive/combine-behaviors.js';
export * from './reactive/disposable-group.js';
export * from './reactive/keyed-selector.js';
export * from './reactive/observable-state.js';
export * from './reactive/readonly-state.js';
export * from './reactive/transaction.js';

export * from './shell/fullscreen.js';
export * from './shell/shell.js';
export * from './shell/status.js';
export * from './shell/theme.js';

export * from './commands/keys.js';
export * from './commands/manager.js';
export * from './commands/types.js';

export { ActionFooter, type ActionFooterProps } from './components/action-footer/action-footer.js';
export { Button, type ButtonProps } from './components/button/button.js';
export { type Choice, Choices } from './components/choices/choices.js';
export { ColorPicker, type ColorPickerProps } from './components/color/color-picker.js';
export { Icon, type IconName, type IconProps } from './components/icon/icon.js';
export { Markdown } from './components/markdown/markdown.js';
export { ParamsForm } from './components/params/params-form.js';
export { ScrollArea, type ScrollAreaProps } from './components/scroll-area/scroll-area.js';
export { Select, type SelectOption, type SelectProps } from './components/select/select.js';
export { TextInput, type TextInputProps } from './components/text-input/text-input.js';
export { flattenTree, VirtualList } from './components/virtual-list/virtual-list.js';

export { type CommandMenuAssistant, CommandMenu, type CommandMenuProps } from './views/command-menu/command-menu.js';
export { StatusLine, StatusPill } from './views/shell/canvas-chrome.js';
export { LogView } from './views/shell/log-view.js';
export { SectionList } from './views/shell/section-list.js';
export { ShellView, type ShellViewProps } from './views/shell/shell-view.js';
