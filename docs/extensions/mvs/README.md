# MolViewSpec extension in Mol\*

**MolViewSpec (MVS)** is a toolkit for standardized description of reproducible molecular visualizations shareable across software applications. MolViewSpec describes a 3D molecular scene by a MolViewSpec State, a viewer-independent representation that contains all information necessary to reproduce the scene. Thanks to its nested tree-based format, complex scenes can be composed from simple building blocks.

MolViewSpec documentation is available from <https://molstar.org/docs/extensions/mvs/>.

MolViewSpec extension in Mol\* provides functionality for building, validating, and visualizing MolViewSpec States.

## Graphical user interface

-   **Drag&drop support:** The easiest way to load a MVS view into Mol* Viewer is to drag a `.mvsj` or `.mvsx` file and drop it in a browser window with Mol* Viewer.

-   **Load via menu:** Another way to load a MVS view is to use "Download File" or "Open Files" action, available in the "Home" tab in the left panel. For these actions, the "Format" parameter must be set to "MVSJ" or "MVSX" (in the "Miscellaneous" category) or "Auto".

-   **URL parameters:** Mol\* Viewer supports `mvs-url`, `mvs-data`, and `mvs-format` URL parameters to specify a MVS view to be loaded when the viewer is initialized.

    -   `mvs-url` specifies the address from which the MVS view should be retrieved.
    -   `mvs-data` specifies the MVS view data directly. Keep in mind that some characters must be escaped to be used in the URL. Also beware that URLs longer than 2000 character may not work in all browsers. Because of these limitations, the preferred method it to host the data somewhere and use `mvs-url` instead.
    -   `mvs-format` specifies the format of the MVS view data from `mvs-url` or `mvs-data`. Allowed values are `mvsj` and `mvsx` (default is `mvsj`).

    Examples of URL parameter usage:

    -   https://molstar.org/viewer/?mvs-format=mvsj&mvs-url=https://raw.githubusercontent.com/molstar/molstar/master/examples/mvs/1cbs.mvsj

    -   https://molstar.org/viewer/?mvs-format=mvsx&mvs-url=https://raw.githubusercontent.com/molstar/molstar/master/examples/mvs/1h9t.mvsx

    -   https://molstar.org/viewer/?mvs-format=mvsj&mvs-data=%7B%22metadata%22%3A%7B%22title%22%3A%22Example%20MolViewSpec%20-%201cbs%20with%20labelled%20protein%20and%20ligand%22%2C%22version%22%3A%221%22%2C%22timestamp%22%3A%222023-11-24T10%3A38%3A17.483%22%7D%2C%22root%22%3A%7B%22kind%22%3A%22root%22%2C%22children%22%3A%5B%7B%22kind%22%3A%22download%22%2C%22params%22%3A%7B%22url%22%3A%22https%3A//www.ebi.ac.uk/pdbe/entry-files/1cbs.bcif%22%7D%2C%22children%22%3A%5B%7B%22kind%22%3A%22parse%22%2C%22params%22%3A%7B%22format%22%3A%22bcif%22%7D%2C%22children%22%3A%5B%7B%22kind%22%3A%22structure%22%2C%22params%22%3A%7B%22type%22%3A%22model%22%7D%2C%22children%22%3A%5B%7B%22kind%22%3A%22component%22%2C%22params%22%3A%7B%22selector%22%3A%22polymer%22%7D%2C%22children%22%3A%5B%7B%22kind%22%3A%22representation%22%2C%22params%22%3A%7B%22type%22%3A%22cartoon%22%7D%2C%22children%22%3A%5B%7B%22kind%22%3A%22color%22%2C%22params%22%3A%7B%22color%22%3A%22green%22%7D%7D%2C%7B%22kind%22%3A%22color%22%2C%22params%22%3A%7B%22selector%22%3A%7B%22label_asym_id%22%3A%22A%22%2C%22beg_label_seq_id%22%3A1%2C%22end_label_seq_id%22%3A50%7D%2C%22color%22%3A%22%236688ff%22%7D%7D%5D%7D%2C%7B%22kind%22%3A%22label%22%2C%22params%22%3A%7B%22text%22%3A%22Protein%22%7D%7D%5D%7D%2C%7B%22kind%22%3A%22component%22%2C%22params%22%3A%7B%22selector%22%3A%22ligand%22%7D%2C%22children%22%3A%5B%7B%22kind%22%3A%22representation%22%2C%22params%22%3A%7B%22type%22%3A%22ball_and_stick%22%7D%2C%22children%22%3A%5B%7B%22kind%22%3A%22color%22%2C%22params%22%3A%7B%22color%22%3A%22%23cc3399%22%7D%7D%5D%7D%2C%7B%22kind%22%3A%22label%22%2C%22params%22%3A%7B%22text%22%3A%22Retinoic%20Acid%22%7D%7D%5D%7D%5D%7D%5D%7D%5D%7D%2C%7B%22kind%22%3A%22canvas%22%2C%22params%22%3A%7B%22background_color%22%3A%22%23ffffee%22%7D%7D%2C%7B%22kind%22%3A%22camera%22%2C%22params%22%3A%7B%22target%22%3A%5B17%2C21%2C27%5D%2C%22position%22%3A%5B41%2C34%2C69%5D%2C%22up%22%3A%5B-0.129%2C0.966%2C-0.224%5D%7D%7D%5D%7D%7D

    -   https://molstar.org/viewer/?mvs-format=mvsx&mvs-data=base64,UEsDBBQAAAAIADSFPlhDx8RXYwEAAGwFAAAUABwAYW5ub3RhdGlvbnMtMWg5dC5jaWZVVAkAA8MmuWVuvbtldXgLAAEE9gEAAAQUAAAAlZJNj4IwEIbv/RVNPPSkacGPckT8TDabPe6emkq72ARbF9iD/34LFBV1FSbh7YTh4Z1hELzgjOyDgnGtTcELZXQOQGrMkQEWm8PRaKmLfJTynUwZz08HpkSrck5BCD8yU0ilwfycRXDxHoJFpUu4NqkAm/pYwWifmkwJCdaXdHtJmy6uOrtt47q0kwmry7n8uatKLZ5UY2M9Hzi1byWZ+T2WbAhtBPaiPoRoIMR0ijGyqZ1yuFNaKJ0geI5RBdAxhJ5PSgBjSmsgjE/pMDIhbFEj+wFbDnHcAP85zO8dKH3mEFWHEzQgBGPOnQP8vHp347BoA6KKBvh6ACzbwLeNeoZy/XfPW2DTF1i1AYyFqIHmH0I3wLovsO0LlJuezJygTDZ9ozeVcC2aLcBcFSVLKmA6c9IRIB5x8hrwamDipDMQOOkKzKiT14BfAh72S+nY0rzf0OMa6D60A7oPXQF9ZhiDP1BLAwQUAAAACACagEFYhgn8CJECAAAbFQAACgAcAGluZGV4Lm12c2pVVAkAAyTBu2UlwbtldXgLAAEE9gEAAAQUAAAA7VjBjtowEL3vV0S5tiSBFnXhttpKPbWq1KqXqoqMPWwsHDuynaVoxb/XNhAIiUOEtK26m5zIzLzJG8+8gfB0EwShFEKH8+DJfDZ3K8qJudtZ3+5sOKOMSODG/tNZgn10DUHEmjOByB7lfAWSKFdV9r21lMwCMq0LNY/j9XodwYJGCEflKi7IAmLgWm5GS8pAxYe88Tib6bQsCNJAIkyXYZVze/LIFrJ1wjXShqCCE8Z+1s6zFDJH9rTCPK8xOGPRxaTJpsZIaVliXcpzVt3MnFdvCnDcBAEWNiK2LQk7OLbzrHHFIi8EN81Kl1LkaSlpC+nLxF2EBZuUUYw4FxppKrga2Za7VrenrXWkMwyboXkQcpNylEONuurITYGRJsIPWDCBV2kGiIC0EDexJ/X4kQpnkCP3mAxRfonTI2IlKE/T9pFfpdBgUnkifrXaW2bEJbswJ/bydfZ0qUAhQZkjdMfhLbLPvLiow8BjJLUQ3lK9ZbkkPUqzVxeRmiSYkJfkUKF6lekir5NHBe8nkyq8IZc+Q3x82JluzIlcBl2tnSrDUUNmyCgpIZWIP4B/Kuy17fD6fe3i8SHaolsn8nWv22fYgR+/3L2G/bdAjBmtkFRpilfDGvTC/+oaHPZRhRj20S7yk2CkS/b3mZENJd4ZeUlLS5VyiXCHHIZtNWwrT/Q124qhBbD/b1P1en2s/+J2lT7vW2qPGWnr5jV9M292TNNi6Ny/7VzDdq7LesSp9+g5WHcWe+cmIsxBI4I0Ov4P+QhS2a8bwziJxod/IgkoLGmh9547HsBvlBcMgjXVWfD5x7egeRyhpjkobeIsZpJM3o+SySgZfx+P59Pb+btZNPtwO5lO3yTJPEnskWxvtn8AUEsBAh4DFAAAAAgANIU+WEPHxFdjAQAAbAUAABQAGAAAAAAAAQAAAKSBAAAAAGFubm90YXRpb25zLTFoOXQuY2lmVVQFAAPDJrlldXgLAAEE9gEAAAQUAAAAUEsBAh4DFAAAAAgAmoBBWIYJ/AiRAgAAGxUAAAoAGAAAAAAAAQAAAKSBsQEAAGluZGV4Lm12c2pVVAUAAyTBu2V1eAsAAQT2AQAABBQAAABQSwUGAAAAAAIAAgCqAAAAhgQAAAAA

## Programming interface

Most functions for manipulation of MVS data (including parsing, encoding, validating, and building) are provided by the `MVSData` object (defined in [src/extensions/mvs/mvs-data.ts](/src/extensions/mvs/mvs-data.ts)). In TypeScript, `MVSData` is also the type for a MVS view.

The `loadMVS` function (defined in [src/extensions/mvs/load.ts](/src/extensions/mvs/load.ts)) can be used to load MVS view data into Mol\* Viewer.

Example usage:

```ts
// Fetch a MVS, validate, and load
const response = await fetch('https://raw.githubusercontent.com/molstar/molstar/master/examples/mvs/1cbs.mvsj');
const rawData = await response.text();
const mvsData: MVSData = MVSData.fromMVSJ(rawData);
if (!MVSData.isValid(mvsData)) throw new Error(`Oh no: ${MVSData.validationIssues(mvsData)}`);
await loadMVS(this.plugin, mvsData, { replaceExisting: true });
console.log('Loaded this:', MVSData.toPrettyString(mvsData));
console.log('Loaded this:', MVSData.toMVSJ(mvsData));

// Build a MVS and load
const builder = MVSData.createBuilder();
const structure = builder.download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1og2_updated.cif' }).parse({ format: 'mmcif' }).modelStructure();
structure.component({ selector: 'polymer' }).representation({ type: 'cartoon' });
structure.component({ selector: 'ligand' }).representation({ type: 'ball_and_stick' }).color({ color: '#aa55ff' });
const mvsData2: MVSData = builder.getState();
await loadMVS(this.plugin, mvsData2, { replaceExisting: false });
```

When using the pre-built Mol\* plugin bundle, `MVSData` and `loadMVS` are exposed as `molstar.PluginExtensions.mvs.MVSData` and `molstar.PluginExtensions.mvs.loadMVS`. Furthermore, the `molstar.Viewer` class has `loadMvsFromUrl` and `loadMvsData` methods, providing the same functionality as `mvs-url` and `mvs-data` URL parameters. See the [integration examples](./integration-examples.html) page for a demonstration.

## Command-line utilities

The MVS extension in Mol\* provides a few command-line utilities, which can be executed via NodeJS:

-   `mvs-validate` provides validation of MolViewSpec files
-   `mvs-render` creates images based on MolViewSpec files
-   `mvs-print-schema` prints MolViewSpec tree schema (i.e. currently supported node types and their parameters)

Example usage:

```sh
# Clone Mol* repo, install, and build
git clone https://github.com/molstar/molstar.git
cd molstar/
npm install && npm run build

# Validate a MolViewSpec file `examples/mvs/1cbs.mvsj`
node lib/commonjs/cli/mvs/mvs-validate.js examples/mvs/1cbs.mvsj

# Render a MolViewSpec file `examples/mvs/1cbs.mvsj` to `../outputs/1cbs.png`
node lib/commonjs/cli/mvs/mvs-render.js -i examples/mvs/1cbs.mvsj -o ../outputs/1cbs.png --size 800x600 --molj

# Print MolViewSpec tree schema formatted as markdown
node lib/commonjs/cli/mvs/mvs-print-schema.js --markdown
```

(An alternative to cloning the GitHub repository is to install Mol\* package from npm by `npm install molstar canvas gl jpeg-js pngjs`. Then you can type `npx mvs-validate ...` instead of `node lib/commonjs/cli/mvs/mvs-validate.js ...`)
