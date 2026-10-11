/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export type Disposable = {
  dispose(): void;
};

export type DisposableGroup = Disposable & {
  /** Removes this owned resource independently; later group disposal will not invoke it again. */
  registerDispose(dispose: () => void): () => void;
};

export function createDisposableGroup(): DisposableGroup {
  type OwnedDispose = { active: boolean; readonly dispose: () => void };

  const owned: OwnedDispose[] = [];
  let disposed = false;
  return {
    registerDispose(dispose) {
      if (disposed) {
        dispose();
        return () => {};
      }
      const entry: OwnedDispose = { active: true, dispose };
      owned.push(entry);
      return () => {
        if (!entry.active) {
          return;
        }
        entry.active = false;
        const index = owned.indexOf(entry);
        if (index >= 0) {
          owned.splice(index, 1);
        }
        dispose();
      };
    },
    dispose() {
      if (disposed) {
        return;
      }
      disposed = true;
      let failed = false;
      let failure: unknown;
      try {
        for (const entry of [...owned].reverse()) {
          if (!entry.active) {
            continue;
          }
          entry.active = false;
          try {
            entry.dispose();
          } catch (error) {
            if (!failed) {
              failed = true;
              failure = error;
            }
          }
        }
      } finally {
        owned.length = 0;
      }
      if (failed) {
        throw failure;
      }
    },
  };
}
