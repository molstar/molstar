## Plugin: custom asset manager

- [Breaking] `AssetManager` (`@molstar/core/util/assets`) is an interface; the implementation is `DefaultAssetManager`
  (`new AssetManager()` becomes `new DefaultAssetManager()`)
- `PluginSpec.assetManager`: a factory for the plugin's `managers.asset`, called once per `PluginContext`; it returns
  `{ manager, disposeWithPlugin? }`, so one manager can be shared by several plugins (`disposeWithPlugin: false`).
  Without it a plugin uses a new `DefaultAssetManager` and disposes it, as before
