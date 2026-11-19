# seqx-render

Renders a SeqViewSpec file to an SVG or a PNG, headless (Seq\* next, `@molstar/seqx-core`). A PNG is the SVG's picture
rasterized: the same figure, at a scale.

```sh
npx tsc -b
node cli/seqx-render/lib/seqx-render.js <spec> -o out.svg
node cli/seqx-render/lib/seqx-render.js <spec> -o out.png --width 960 --theme dark --scale 2
```

| Option                | Meaning                                                                           |
| --------------------- | --------------------------------------------------------------------------------- |
| `-o, --output <file>` | `.svg` or `.png` (required)                                                       |
| `--width <px>`        | the figure's width in CSS px (960)                                                |
| `--height <px\|auto>` | its height in CSS px (`auto`: the content's)                                      |
| `--dpr <n>`           | the data's images at this device pixel ratio, and the PNG's pixels per CSS px (1) |
| `--scale <n>`         | PNG pixels per CSS px, times `dpr` (1; `dpr` × `scale` is at most 8)              |
| `--theme <t>`         | `light` (default) or `dark`: the built-in token table                             |
| `--view <id>`         | the view's address (the spec's first)                                             |
| `--allow-host <host>` | an `https` host the spec's data may be read from (repeatable)                     |
| `--timeout <ms>`      | how long to wait for the data (10000)                                             |

The spec is JSON, an S-expression or MessagePack, told by its content, not its extension. The SVG paints its own
background (the theme's), so a PNG has no clear pixels.

## Assets

A spec's `uri`s are untrusted input, so they are read under the headless policy (`headlessAssetPolicy` of
`@molstar/seqx-core/context`): files under the spec's folder only (a relative location stays in it; `..`, an absolute
path and `file:` outside it are refused), and no network unless `--allow-host` names a host (`https` only; no
credentials, redirects refused, the browser's size and time bounds).

## Output and exit codes

The issues of the spec, its data and the render go to stderr, one a line (`severity code place: message (hint)`);
nothing goes to stdout but `--help`. Nothing is written when the run fails.

| Code | Meaning                                                               |
| ---- | --------------------------------------------------------------------- |
| 0    | the file is written (warnings may be printed)                         |
| 1    | the spec cannot be read as a spec, or it, its data or view has errors |
| 2    | the command line is wrong (the usage is printed)                      |
| 3    | the data did not arrive before the timeout                            |
| 4    | a file could not be read or written                                   |

## PNG

`renderPng` of `@molstar/seqx-core/render` takes the rasterizer as a port (`Rasterizer`: an SVG string and its size and
scale to PNG bytes), so the core carries no native code. This package supplies it over
[resvg](https://github.com/RazrFalcon/resvg) (`@resvg/resvg-js`, MPL-2.0, a prebuilt binary, used unmodified). It loads
**no system fonts**: the text is drawn with the bundled Liberation Sans and Liberation Mono (`fonts/`, regular and bold,
SIL OFL 1.1, `fonts/OFL.txt`), the metric-compatible faces of the Arial advances and the fixed 600 advance of the text
metrics the layouts measure with. The same SVG draws the same on every machine.

## Cold start

From the process's start to the written file for p53's tracks: under 1 s for the SVG (the performance target) and for
the PNG; about 0.12 s each on the development machine (Node 22, Apple silicon). Held by `cli-process.test.ts` in the
heavy lane (`npx tsc -b`, then `pnpm test:heavy cli/seqx-render`).
