/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectCell } from '@molstar/core/state';
import type { PluginContext } from '@molstar/plugin/context';
import { BehaviorSubject } from 'rxjs';
import { BuiltInMarkdownExtension } from '../markdown/catalog.js';

export type MarkdownExtensionEvent = 'click' | 'mouse-enter' | 'mouse-leave';

export interface MarkdownExtension {
  name: string;
  execute?: (options: {
    event: MarkdownExtensionEvent;
    args: Record<string, string>;
    manager: MarkdownExtensionManager;
  }) => void;
  reactRenderFn?: (options: { args: Record<string, string>; manager: MarkdownExtensionManager }) => any;
}

export class MarkdownExtensionManager {
  state = {
    audioPlayer: new BehaviorSubject<HTMLAudioElement | null>(null),
  };

  private extension: MarkdownExtension[] = [];
  /** Registration counts by extension name */
  private extensionCounts = new Map<string, number>();
  private refResolvers: Record<string, (plugin: PluginContext, refs: string[]) => StateObjectCell[]> = {
    default: (plugin: PluginContext, refs: string[]) =>
      refs.map((ref) => plugin.state.data.cells.get(ref)).filter((c) => !!c),
  };
  private uriResolvers: Record<string, (plugin: PluginContext, uri: string) => Promise<string> | string | undefined> =
    {};
  private argsParsers: [
    name: string,
    priority: number,
    parser: (input: string | undefined) => Record<string, string> | undefined,
  ][] = [['default', 100, defaultParseMarkdownCommandArgs]];

  /**
   * Default parser has priority 100, parsers with higher priority
   * will be called first.
   */
  registerArgsParser(
    name: string,
    priority: number,
    parser: (input: string | undefined) => Record<string, string> | undefined,
  ) {
    this.removeArgsParser(name);
    this.argsParsers.push([name, priority, parser]);
    this.argsParsers.sort((a, b) => b[1] - a[1]); // Sort by priority, higher first
  }

  removeArgsParser(name: string) {
    const idx = this.argsParsers.findIndex((p) => p[0] === name);
    if (idx >= 0) {
      this.argsParsers.splice(idx, 1);
    }
  }

  parseArgs(input: string | undefined): Record<string, string> | undefined {
    for (const [, , parser] of this.argsParsers) {
      const ret = parser(input);
      if (ret) return ret;
    }
    return undefined;
  }

  registerRefResolver(name: string, resolver: (plugin: PluginContext, refs: string[]) => StateObjectCell[]) {
    this.refResolvers[name] = resolver;
  }

  removeRefResolver(name: string) {
    delete this.refResolvers[name];
  }

  registerUriResolver(
    name: string,
    resolver: (plugin: PluginContext, uri: string) => Promise<string> | string | undefined,
  ) {
    this.uriResolvers[name] = resolver;
  }

  removeUriResolver(name: string) {
    delete this.uriResolvers[name];
  }

  /** Returns the message `registerExtension` would throw for `command`, without changing anything. */
  findConflict(command: MarkdownExtension): string | undefined {
    const existing = this.extension.find((c) => c.name === command.name);
    if (existing && existing !== command) {
      return `MarkdownExtensionManager: a different extension is already registered under the name '${command.name}'.`;
    }
    return undefined;
  }

  /**
   * Extensions are keyed by `name`. Registering the same object again increments a count,
   * a different object under an existing name throws.
   */
  registerExtension(command: MarkdownExtension) {
    const conflict = this.findConflict(command);
    if (conflict) throw new Error(conflict);

    const count = this.extensionCounts.get(command.name);
    if (count !== undefined) {
      this.extensionCounts.set(command.name, count + 1);
    } else {
      this.extensionCounts.set(command.name, 1);
      this.extension.push(command);
    }
  }

  /**
   * Decrements the count of the extension (given by name or object) and removes it at zero.
   * No-op for an unknown name, or for an object that is not the one registered under its name.
   */
  removeExtension(nameOrExtension: string | MarkdownExtension) {
    const name = typeof nameOrExtension === 'string' ? nameOrExtension : nameOrExtension.name;
    const idx = this.extension.findIndex((c) => c.name === name);
    if (idx < 0) return;
    if (typeof nameOrExtension !== 'string' && this.extension[idx] !== nameOrExtension) return;

    const count = this.extensionCounts.get(name) ?? 1;
    if (count > 1) {
      this.extensionCounts.set(name, count - 1);
    } else {
      this.extensionCounts.delete(name);
      this.extension.splice(idx, 1);
    }
  }

  private _tryRender(
    ext: MarkdownExtension,
    options: { args: Record<string, string>; manager: MarkdownExtensionManager },
  ) {
    try {
      return ext.reactRenderFn?.(options);
    } catch (e) {
      console.error(`Failed to render markdown extension '${ext.name}'`, e);
      return null;
    }
  }

  /**
   * Render a markdown extension with the given arguments.
   * Default renderers are defined separately because we
   * don't want to include `react` outside of mol-plugin-ui.
   */
  tryRender(args: Record<string, string>, defaultRenderers: MarkdownExtension[]): any {
    const options = { args, manager: this };
    for (const ext of this.extension) {
      const ret = this._tryRender(ext, options);
      if (ret) {
        return ret;
      }
    }
    for (const ext of defaultRenderers) {
      const ret = this._tryRender(ext, options);
      if (ret) {
        return ret;
      }
    }
    return null;
  }

  tryExecute(event: MarkdownExtensionEvent, args: Record<string, string>) {
    const options = { event, args, manager: this };
    for (const ext of this.extension) {
      try {
        ext.execute?.(options);
      } catch (e) {
        console.error(`Failed to execute markdown extension '${ext.name}'`, e);
      }
    }
  }

  tryResolveUri(uri: string): Promise<string> | string | undefined {
    for (const resolver of Object.values(this.uriResolvers)) {
      const resolved = resolver(this.plugin, uri);
      if (resolved) {
        return resolved;
      }
    }
  }

  findCells(refs: string[]): StateObjectCell[] {
    const added = new Set<string>();
    const ret: StateObjectCell[] = [];
    for (const resolver of Object.values(this.refResolvers)) {
      for (const cell of resolver(this.plugin, refs)) {
        if (cell && !added.has(cell.transform.ref)) {
          added.add(cell.transform.ref);
          ret.push(cell);
        }
      }
    }
    return ret;
  }

  private resolveAudioPlayer() {
    if (this.state.audioPlayer.value) {
      return this.state.audioPlayer.value;
    }

    const audio = document.createElement('audio');
    audio.controls = true;
    audio.preload = 'auto';
    audio.style.width = '100%';
    audio.style.height = '32px';
    this.state.audioPlayer.next(audio);
    return audio;
  }

  get audioPlayer() {
    return this.state.audioPlayer.value;
  }

  audio = {
    play: async (src: string, options?: { toggle?: boolean }) => {
      try {
        const audio = this.resolveAudioPlayer();

        let newSource = false;
        if (src?.trim()) {
          const resolved = this.tryResolveUri(src);
          let uri: string = src;
          if (typeof (resolved as Promise<string>)?.then === 'function') {
            uri = (await resolved) as string;
          } else if (resolved) {
            uri = resolved as string;
          }
          newSource = audio.src !== uri;
          if (newSource) {
            audio.src = uri;
            audio.load();
          }
        }

        if (!newSource && options?.toggle) {
          if (audio.paused) {
            await audio.play();
          } else {
            audio.pause();
          }
        } else {
          audio.currentTime = 0;
          await audio.play();
        }
      } catch (e) {
        console.error('Failed to play audio', e);
      }
    },
    pause: () => {
      this.audioPlayer?.pause();
    },
    stop: () => {
      if (!this.audioPlayer) return;
      this.audioPlayer.pause();
      this.audioPlayer.currentTime = 0;
    },
    dispose: () => {
      if (this.audioPlayer) {
        this.audioPlayer.pause();
        this.audioPlayer.currentTime = 0;
        this.state.audioPlayer.next(null);
      }
    },
  };

  constructor(public plugin: PluginContext) {
    for (const command of BuiltInMarkdownExtension) {
      this.registerExtension(command);
    }
  }
}

export function defaultParseMarkdownCommandArgs(input: string | undefined): Record<string, string> | undefined {
  if (!input?.startsWith('!')) return undefined;
  const entries = decodeURIComponent(input.substring(1))
    .split('&')
    .map((p) => p.trim())
    .filter((p) => p.length > 0)
    .map((p) => p.split('=', 2).map((s) => s.trim()));
  if (entries.length === 0) return undefined;
  return Object.fromEntries(entries);
}
