# Components Example

A proof of concept for separating Mol* into a ready-to-use web component library.

Implements loading of MolViewSpec states and rendering snapshot markdown descriptions in a separate component.

Usage:

- Clone Mol* GitHub repo and build it.
```bash
  git clone https://github.com/molstar/molstar.git
  cd molstar
  npm install
  npm build
```

- Get `molstar.css` and `index.js` from `build/examples/components` and include these to your HTML page

```html
    <link rel="stylesheet" type="text/css" href="molstar.css" />
    <script type="text/javascript" src="index.js"></script>
```

- Plate the components in your page wrapper in `<div>` elements to set up positioning:

```html
<div class="viewer">
    <mc-viewer name="v1" />
</div>
<div class="snapshot">
    <mc-snapshot-markdown viewer-name="v1" />
</div>
```

- Load MolViewSpec state:

```html
<script>
    window.mc.getContext().dispatch({
        kind: 'load-mvs',
        format: 'mvsj',
        url: 'https://path/to/file.mvsj',
        // or provide data direcly
        // data: mvsJSON
    });
</script>
```

See [index.html](./index.html) for a full example.