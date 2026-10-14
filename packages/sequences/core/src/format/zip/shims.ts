// What the Mol* zip files import from outside their folder (see VENDORED.md): `mol-task`'s `RuntimeContext` and
// `Task`, `type-helpers`' `NumberArray` and `assertUnreachable`. Each type is also a value, so the files' plain
// `import { X }` (as Mol* writes them) loads under Node's type stripping.

type TypedArray =
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | Float32Array
  | Float64Array;

export type NumberArray = TypedArray | number[];
export const NumberArray = undefined;

/** A context that never asks to yield: the zip runs straight through. */
export interface RuntimeContext {
  readonly shouldUpdate: boolean;
  readonly isSynchronous: boolean;
  update(progress?: unknown, dontNotify?: boolean): Promise<void> | void;
}
export const RuntimeContext: RuntimeContext = { shouldUpdate: false, isSynchronous: true, update() {} };

/** `Task.create(name, f)` runs `f` once, at once, with no progress. */
export interface Task<T> {
  run(): Promise<T>;
}
export const Task = {
  create<T>(_name: string, f: (ctx: RuntimeContext) => Promise<T>): Task<T> {
    return { run: () => f(RuntimeContext) };
  },
};

export function assertUnreachable(_x: never): never {
  throw new Error('unreachable');
}
