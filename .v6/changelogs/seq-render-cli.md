## Seq* (next): PNG and the render CLI

- `renderPng`: the SVG export rasterized through a `Rasterizer` port, so a PNG is the SVG's picture; images in the SVG
  are drawn nearest-neighbour (`image-rendering="optimizeSpeed"`).
- `seqx-render` (`cli/seqx-render`, private): a spec file (JSON, S-expression or MessagePack) to SVG or PNG over resvg
  with bundled Liberation fonts; assets under the spec's folder (links leading outside refused), no network unless
  `--allow-host`; issues and refusal reasons on stderr; exit codes 0 done, 1 spec, 2 usage, 3 unsettled, 4 file, 5
  internal; atomic writes and a pixel cap.
