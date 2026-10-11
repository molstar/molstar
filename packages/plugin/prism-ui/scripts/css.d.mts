/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export function compile(file: string): { out: string; content: string; prefix?: string };
export function generateAll(options?: { check?: boolean }): string[];
export function updateFile(file: string): void;
export function cssModulesPlugin(): {
  name: string;
  buildStart(this: { error(message: string): never }): void;
  configureServer(server: {
    watcher: { on(event: string, fn: (file: string) => void): void };
    config: { logger: { error(message: string): void } };
  }): void;
};
