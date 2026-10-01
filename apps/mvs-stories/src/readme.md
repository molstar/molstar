# MolViewSpec Stories App

An app that defines `mvs-stories-snapshot-markdown` and `mvs-stories-viewer` web components that can be used to view MolViewSpec molecular stories.

See the [mvs-stories](../../examples/mvs-stories) example that includes specific stories.

### Usage

- Get `mvs-stories.css` and `mvs-stories.js` from `build/mvs-stories` and include these to your HTML page

```html
<link rel="stylesheet" type="text/css" href="mvs-stories.css" />
<script type="text/javascript" src="mvs-stories.js"></script>
```

Can also use `https://cdn.jsdelivr.net/npm/molstar@latest/build/mvs-stories/mvs-stories.js` (and `.css`). `latest` can be substituted by specific version.

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
 
- See [index.html](./index.html) for full example of how to embed the app.

- For interactive development build (for production use `npm run build`) of the example that immediately reflects changes use:

```bash
npm run dev -- -a mvs-stories
```

### Multiple Stories on a Single Page

To support multiple instances of stories, use the `context-name='unique-name'` attribute  on the `mvs-` components together with `loadFromURL/Data(..., { contextName: 'unique-name' })`.

For example (simplified to not include layout):

```html
<div>
    <mvs-stories-viewer context-name="1" />
    <mvs-stories-snapshot-markdown context-name="1" />
</div>
<div>
    <mvs-stories-viewer context-name="2" />
    <mvs-stories-snapshot-markdown context-name="2" />
</div>

<script>
    mvsStories.loadFromURL('1.mvsj', { format: 'mvsj', contextName: '1' });
    mvsStories.loadFromURL('2.mvsj', { format: 'mvsj', contextName: '2' });
</script>

```