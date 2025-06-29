# MolViewSpec Stories Example

This example illustrates using the `mvs-stories` app to tell molecular stories built with MolViewSpec.

See the [mvs-stories](../../apps/mvs-stories) app for more info about how to use this app separately.

### Usage

- Clone Mol* GitHub repo and build it.
```bash
  git clone https://github.com/molstar/molstar.git
  cd molstar
  npm install
  npm build
```

- See [index.html](./index.html) for example usage.

- For interactive development build (for production use `npm run build`) of the example that immediately reflects changes use:

```bash
  npm run dev -- -e mvs-stories
```