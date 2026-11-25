/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The command layer (`_spec/README.md`): DOM-free. Named re-exports only (no `export *`: the cycle rule).

export { type BrowseHost, type CommandBrowse, CommandBrowser } from './browse.js';
export { CommandInvocation } from './invocation.js';
export { InvocationModel } from './invocation-model.js';
export type { InvocationError, InvocationInit, InvocationState, InvocationStatus } from './invocation-state.js';
export { type CommandItem, type CommandItemOptions, commandItem, commandItems } from './items.js';
export {
  chordFromEvent,
  chordId,
  chordKeys,
  ariaKeyShortcut,
  chordMatches,
  formatChord,
  isMac,
  isSymbolKey,
  KeyPlatform,
  type KeyEventLike,
} from './keys.js';
export {
  DefaultKeysLabels,
  type InputBinding,
  type InputBindingGroup,
  type KeysLabels,
  KeysModel,
} from './keys-model.js';
export {
  type KeyConflict,
  KeyIndex,
  Keymap,
  type KeymapLayers,
  type KeymapOptions,
  ScopeAttribute,
  scopeOf,
} from './keymap.js';
export {
  type CommandManagerOptions,
  CommandManager,
  type CommandRun,
  type OpenResult,
} from './manager.js';
export { checkFormat } from './format.js';
export { fuzzyScore, matchPositions, rankBy } from './match.js';
export { aliases, type CommandMatch, type CommandQuery } from './query.js';
export { type LazyCommands, type LazyRegistration } from './registry.js';
export type {
  ChangedOption,
  FieldChoice,
  FieldChoices,
  FieldRange,
  InvocationOptions,
  SentenceField,
} from './sentence.js';
export type {
  AssistRequest,
  AssistResult,
  CommandArgs,
  CommandAssistant,
  CommandDefinition,
  CommandFieldHint,
  CommandFormat,
  CommandFormats,
  CommandGroup,
  CommandMenuItem,
  CommandMenuPrefix,
  CommandOutcome,
  CommandPreset,
  CommandResult,
  CommandValue,
  CommandValues,
  FormatCheck,
  FormatCompletion,
  KeyBindings,
  KeyChord,
  KeymapStore,
  PreviewSink,
  RecentsStore,
  SharedArgument,
} from './types.js';
