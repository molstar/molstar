# Drift: the viewer and the views (lane H)

1. **One window for all leading sections**: `Renderer.frame` takes one `viewport`; the `ViewModel` keeps one per
   section. Options: a per-section viewport map in `FrameInput`, or one window. Recommend the map before any view has
   two independent leaders; the viewer uses the first leader's window now.
2. **Projected hover is not drawn**: `Emphasis.hover` is one `Hit`, so the `HoverModel`'s projection onto linked
   sections cannot reach the draw lists. Recommend `Emphasis.projected?: (section) => Location`, drawn in the `hover`
   pass.
3. **A source-backed leader's length**: no pieces in `ResolvedView`, so the `ViewModel` stays pending. Recommend the
   viewer feeds coverage lengths to the `ViewModel` (an axis update).
4. **Registry kinds**: `renderSvg` has no `types` option, so the registry holds format runtimes only in this slice
   (a layer type the canvas drew but the export could not would make the pictures differ).
5. **Asset size limit** is checked after the whole body is read: `AssetManager.resolve` has no signal or size hook. A
   pre-read bound needs a manager extension (with `Content-Length` and a signal): a core API change, the author's call.
6. **The assets port** takes `{ manager }` or `{ reader }` (a test or CLI has a reader, not a manager). Recommend the
   design say so.
7. **Keys**: matched against the commands' shortcuts, scoped to the container, until a command manager and the
   `@molstar/ui` keymap are wired.
8. **Fixtures' place**: `examples.md` says the storybook package; `examples-and-stories.md` (followed) says
   `packages/seq/examples/fixtures`. Fix `examples.md`.
9. **EC-143** is not in `edge-cases.md`; built as: re-hit once per frame in which the layout was rebuilt.

## From the review (strong)

10. **The asset manager's fetch** (`@molstar/core` `data-source`): it sends cookies on relative reads (`credentials`
    default), follows redirects to disallowed hosts, ignores the signal, and checks size only after reading the whole
    body. security.md asks for no credentials and refusal before allocating. The fix is a core API change
    (`credentials: 'omit'`, `redirect: 'error'`, the signal passed through, a streamed byte limit): **the author's
    call** (base-library.md). Until then the seq policy says what it does not enforce.
11. **Canvas and SVG pass order**: labels above emphasis on both (a third canvas), with a test comparing the two
    backends' order.
