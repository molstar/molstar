# I/HM Restraints Example

This example illustrates:

- Using Mol* to parse CIF files and extract custom data
- Using MolViewSpec to visualize I/HM structure along with annotated crosslink restraints

### Usage

- Clone Mol* GitHub repo and build it.
```bash
  git clone https://github.com/molstar/molstar.git
  cd molstar
  npm install
  npm build
```

- Get `molstar.css` and `index.js` from `build/examples/ihm-restraints` and include these to your HTML page in a similar fashion to [index.html](./index.html):

```html
    <link rel="stylesheet" type="text/css" href="molstar.css" />
    <script type="text/javascript" src="index.js"></script>
    ...
    <div id="viewer"></div>
    ...
    <script>
        loadIHMRestraints(document.getElementById('viewer'))
    </script>
```

- For interactive development build (for production use `npm run build`) of the example that immediately reflects changes use:

```bash
  npm run dev -- -e ihm-restraints
```