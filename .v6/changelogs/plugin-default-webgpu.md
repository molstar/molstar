## Graphics

- The plugin default (`PluginConfig.General.GraphicsBackends` unset) follows the Viewer: WebGPU where the adapter has the
  features WebGL2 has, else WebGL2, then WebGL1 (`@molstar/graphics/gl/backends/default` is the `all` set; replaces "The
  plugin and library default stays WebGL"). Builds with one backend alias the default to
  `@molstar/graphics/gl/backends/webgl` or `@molstar/graphics/gl/backends/webgpu`
- [Breaking] WebGL1 preferences are removed: `PluginConfig.General.PreferWebGl1` (replaces "still selects WebGL1 when
  set"; list only its backend, `PluginConfig.General.GraphicsBackends: [WebGLBackend]`), `Canvas3DContext.Attribs`'
  `preferWebGl1`, the Prism viewer's `renderer.preferWebgl1` and `prefer-webgl1` URL parameter (replaces "keeps v5's
  `prefer-webgl1`"; `graphics-backend=webgl1`)
