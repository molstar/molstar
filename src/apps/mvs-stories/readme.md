# MolViewSpec Stories App

An app that defines `mvs-stories-snapshot-markdown` and `mvs-stories-viewer` web components that can be used to view MolViewSpec molecular stories.

See the [mvs-stories](../../examples/mvs-stories) example that includes specific stories.

### Usage

- Get `molstar.css` and `molstar.js` from `build/mvs-stories` and include these to your HTML page

```html
<link rel="stylesheet" type="text/css" href="molstar.css" />
<script type="text/javascript" src="molstar.js"></script>
```

Can also use `https://cdn.jsdelivr.net/npm/molstar@latest/build/mvs-stories/molstar.js` (and `.css`). `latest` can be substituted by specific version.

- Place the components in your page wrapper in `<div>` elements to set up positioning:

```html
<div class="viewer">
    <mvs-stories-viewer />
</div>
<div class="snapshot">
    <mvs-stories-snapshot-markdown />
</div>
```

- Load MolViewSpec state:

```html
<script>
mvsStories.loadFromURL('https://raw.githubusercontent.com/molstar/molstar/master/examples/mvs/1cbs.mvsj');
</script>
```

- For interactive development build (for production use `npm run build`) of the example that immediately reflects changes use:

```bash
  npm run dev -- -a mvs-stories
```
