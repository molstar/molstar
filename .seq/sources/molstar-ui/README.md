# @molstar/core/reactive and @molstar/ui

The runtime the viewer is built on (this repository: `packages/core/src/reactive` with its `_spec`, `packages/ui` with
a `_spec` per folder). Seq\* is the first package built on `@molstar/ui`.

## What we take

- **States, not streams**: `State` (`.value`, `set`, `asReadonly`), `RecordState`, `derive` (explicit sources, lazy,
  glitch-free), `batch`, `toState`. Streams (`Observable`, `Subject`, operators) only where something is truly an event
  stream (pointer input, data arriving).
- **Keyed state** for per-position and per-row state without fan-out: `KeyedState` (a key's followers only),
  `KeyedSelector` (hover, current row: only the two keys involved are told).
- **`ReactiveModel`**: `mount`, `unmount`, `retire`; mount-scoped `own`, `listen`, `subscribe`; private `#state` fields,
  read-only getters, every change a method, multi-state writes in `batch`; timers and frames injected so the model runs
  in Node.
- **Models DOM-free, views optional** (`packages/ui/src/_test/headless.test.ts` enforces it there): every viewer model
  is tested in plain Node; views opt into jsdom.
- **`@molstar/ui/dom`**: JSX without a virtual DOM, `watch`, `For`, `Show`, scopes, `drag`, `observeResize`, frames.
- **Components and modules**: commands and keymap, tooltip, context menu, virtual list, scroll area, legend, color
  picker, hover label, canvas-ink (pick light or dark ink from a canvas color), history, status, the styles registry and
  tokens (`--msui-*`), the CSS-module generator (`packages/ui/scripts/css.mjs`, `--prefix`).

## What we leave

- `shell` (a whole app shell) and `stage`'s chrome (a 3D stage's): not a sequence viewer's needs; `stage/canvas-ink`'s
  luminance is replaced by `Color.luminance` ([base-library](../../design/base-library.md)).
- `forms` (the params form over `ParamDefinition`): later, when settings get a generated form
  ([extensibility](../../design/extensibility.md)).
- `testing/*` and `_stories/support.tsx` are not exported; copy what is needed.

## Reconciliation

- §frontend-state as written: the view is thin, models own state, simple actions are bound closures.
- The oz reactive-models bootstrap pattern is replaced by `@molstar/core/reactive`, which the host provides
  (manifest, Adaptations).
