/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface CssPackage {
  root: string;
  srcRoot: string;
  prefix?: string;
  registry: string;
}

export function cssPackage(root?: string, prefix?: string): CssPackage;
export function compile(file: string, pkg?: CssPackage): { out: string; content: string; prefix?: string };
export function generateAll(options?: { check?: boolean; pkg?: CssPackage }): string[];
export function updateFile(file: string, pkg?: CssPackage): void;
export function cssModulesPlugin(packages?: CssPackage[]): {
  name: string;
  buildStart(this: { error(message: string): never }): void;
  configureServer(server: {
    watcher: { on(event: string, fn: (file: string) => void): void };
    config: { logger: { error(message: string): void } };
  }): void;
};
