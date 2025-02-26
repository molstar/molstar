# MolViewSpec Kinase Story Example

This example illustrates:

- Using MolViewSpec to tell a story
- A proof of concept for separating Mol* into a ready-to-use web component library.

### Usage

- Clone Mol* GitHub repo and build it.
```bash
  git clone https://github.com/molstar/molstar.git
  cd molstar
  npm install
  npm build
```

- Get `molstar.css` and `index.js` from `build/examples/mvs-kinase-story` and include these to your HTML page

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