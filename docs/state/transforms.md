# Mol* Plugin State Transformer Reference

* [build-in.root](#build-in-root)
* [ms-plugin.download](#ms-plugin-download)
* [ms-plugin.read-file](#ms-plugin-read-file)
* [ms-plugin.parse-cif](#ms-plugin-parse-cif)
* [ms-plugin.parse-ccp4](#ms-plugin-parse-ccp4)
* [ms-plugin.parse-dsn6](#ms-plugin-parse-dsn6)
* [ms-plugin.trajectory-from-mmcif](#ms-plugin-trajectory-from-mmcif)
* [ms-plugin.trajectory-from-pdb](#ms-plugin-trajectory-from-pdb)
* [ms-plugin.model-from-trajectory](#ms-plugin-model-from-trajectory)
* [ms-plugin.structure-from-model](#ms-plugin-structure-from-model)
* [ms-plugin.structure-assembly-from-model](#ms-plugin-structure-assembly-from-model)
* [ms-plugin.structure-symmetry-from-model](#ms-plugin-structure-symmetry-from-model)
* [ms-plugin.structure-selection](#ms-plugin-structure-selection)
* [ms-plugin.structure-complex-element](#ms-plugin-structure-complex-element)
* [ms-plugin.custom-model-properties](#ms-plugin-custom-model-properties)
* [ms-plugin.volume-from-ccp4](#ms-plugin-volume-from-ccp4)
* [ms-plugin.volume-from-dsn6](#ms-plugin-volume-from-dsn6)
* [ms-plugin.representation-highlight-loci](#ms-plugin-representation-highlight-loci)
* [ms-plugin.representation-select-loci](#ms-plugin-representation-select-loci)
* [ms-plugin.default-loci-label-provider](#ms-plugin-default-loci-label-provider)
* [ms-plugin.structure-representation-3d](#ms-plugin-structure-representation-3d)
* [ms-plugin.explode-structure-representation-3d](#ms-plugin-explode-structure-representation-3d)
* [ms-plugin.volume-representation-3d](#ms-plugin-volume-representation-3d)
* [ms-plugin.focus-loci-on-select](#ms-plugin-focus-loci-on-select)
* [ms-plugin.pdbe-structure-quality-report-prop](#ms-plugin-pdbe-structure-quality-report-prop)
* [ms-plugin.rcsb-assembly-symmetry-prop](#ms-plugin-rcsb-assembly-symmetry-prop)
* [ms-plugin.structure-animation](#ms-plugin-structure-animation)
* [ms-plugin.scene-labels](#ms-plugin-scene-labels)

----------------------------
## <a name="build-in-root"></a>build-in.root :: () -> ()
*For internal use.*

----------------------------
## <a name="ms-plugin-download"></a>ms-plugin.download :: Root -> String | Binary
*Download string or binary data from the specified URL*

### Parameters
- **url**: String *(Resource URL. Must be the same domain or support CORS.)*
- **label**?: String
- **isBinary**?: true/false *(If true, download data as binary (string otherwise))*

### Default Parameters
```js
{
  "url": "https://www.ebi.ac.uk/pdbe/static/entry/1cbs_updated.cif"
}
```
----------------------------
## <a name="ms-plugin-read-file"></a>ms-plugin.read-file :: Root -> String | Binary
*Read string or binary data from the specified file*

### Parameters
- **file**: JavaScript File Handle
- **label**?: String
- **isBinary**?: true/false *(If true, open file as as binary (string otherwise))*

### Default Parameters
```js
{}
```
----------------------------
## <a name="ms-plugin-parse-cif"></a>ms-plugin.parse-cif :: String | Binary -> Cif
*Parse CIF from String or Binary data*

----------------------------
## <a name="ms-plugin-parse-ccp4"></a>ms-plugin.parse-ccp4 :: Binary -> Ccp4
*Parse CCP4/MRC/MAP from Binary data*

----------------------------
## <a name="ms-plugin-parse-dsn6"></a>ms-plugin.parse-dsn6 :: Binary -> Dsn6
*Parse CCP4/BRIX from Binary data*

----------------------------
## <a name="ms-plugin-trajectory-from-mmcif"></a>ms-plugin.trajectory-from-mmcif :: Cif -> Trajectory
*Identify and create all separate models in the specified CIF data block*

### Parameters
- **blockHeader**?: String *(Header of the block to parse. If none is specifed, the 1st data block in the file is used.)*

### Default Parameters
```js
{}
```
----------------------------
## <a name="ms-plugin-trajectory-from-pdb"></a>ms-plugin.trajectory-from-pdb :: String -> Trajectory

----------------------------
## <a name="ms-plugin-model-from-trajectory"></a>ms-plugin.model-from-trajectory :: Trajectory -> Model
*Create a molecular structure from the specified model.*

### Parameters
- **modelIndex**: Numeric value *(Zero-based index of the model)*

### Default Parameters
```js
{
  "modelIndex": 0
}
```
----------------------------
## <a name="ms-plugin-structure-from-model"></a>ms-plugin.structure-from-model :: Model -> Structure
*Create a molecular structure from the specified model.*

----------------------------
## <a name="ms-plugin-structure-assembly-from-model"></a>ms-plugin.structure-assembly-from-model :: Model -> Structure
*Create a molecular structure assembly.*

### Parameters
- **id**?: String *(Assembly Id. Value 'deposited' can be used to specify deposited asymmetric unit.)*

### Default Parameters
```js
{}
```
----------------------------
## <a name="ms-plugin-structure-symmetry-from-model"></a>ms-plugin.structure-symmetry-from-model :: Model -> Structure
*Create a molecular structure symmetry.*

### Parameters
- **ijkMin**: 3D vector [x, y, z]
- **ijkMax**: 3D vector [x, y, z]

### Default Parameters
```js
{
  "ijkMin": [
    -1,
    -1,
    -1
  ],
  "ijkMax": [
    1,
    1,
    1
  ]
}
```
----------------------------
## <a name="ms-plugin-structure-selection"></a>ms-plugin.structure-selection :: Structure -> Structure
*Create a molecular structure from the specified query expression.*

### Parameters
- **query**: Value
- **label**?: String

### Default Parameters
```js
{}
```
----------------------------
## <a name="ms-plugin-structure-complex-element"></a>ms-plugin.structure-complex-element :: Structure -> Structure
*Create a molecular structure from the specified model.*

### Parameters
- **type**: One of 'atomic-sequence', 'water', 'atomic-het', 'spheres'

### Default Parameters
```js
{
  "type": "atomic-sequence"
}
```
----------------------------
## <a name="ms-plugin-custom-model-properties"></a>ms-plugin.custom-model-properties :: Model -> Model

### Parameters
- **properties**: Array of  *(A list of property descriptor ids.)*

### Default Parameters
```js
{
  "properties": []
}
```
----------------------------
## <a name="ms-plugin-volume-from-ccp4"></a>ms-plugin.volume-from-ccp4 :: Ccp4 -> Data
*Create Volume from CCP4/MRC/MAP data*

### Parameters
- **voxelSize**: 3D vector [x, y, z]

### Default Parameters
```js
{
  "voxelSize": [
    1,
    1,
    1
  ]
}
```
----------------------------
## <a name="ms-plugin-volume-from-dsn6"></a>ms-plugin.volume-from-dsn6 :: Dsn6 -> Data
*Create Volume from DSN6/BRIX data*

### Parameters
- **voxelSize**: 3D vector [x, y, z]

### Default Parameters
```js
{
  "voxelSize": [
    1,
    1,
    1
  ]
}
```
----------------------------
## <a name="ms-plugin-representation-highlight-loci"></a>ms-plugin.representation-highlight-loci :: Root -> Behavior

----------------------------
## <a name="ms-plugin-representation-select-loci"></a>ms-plugin.representation-select-loci :: Root -> Behavior

----------------------------
## <a name="ms-plugin-default-loci-label-provider"></a>ms-plugin.default-loci-label-provider :: Root -> Behavior

----------------------------
## <a name="ms-plugin-structure-representation-3d"></a>ms-plugin.structure-representation-3d :: Structure -> Representation3D

### Parameters
- **type**: Object { name: string, params: object } where name+params are:
  - **cartoon**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **doubleSided**: true/false
      - **flipSided**: true/false
      - **flatShaded**: true/false
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'
      - **sizeFactor**: Numeric value
      - **linearSegments**: Numeric value
      - **radialSegments**: Numeric value
      - **aspectRatio**: Numeric value
      - **arrowFactor**: Numeric value
      - **visuals**: Array of 'polymer-trace', 'polymer-gap', 'nucleotide-block', 'direction-wedge'

  - **ball-and-stick**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **doubleSided**: true/false
      - **flipSided**: true/false
      - **flatShaded**: true/false
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'
      - **sizeFactor**: Numeric value
      - **detail**: Numeric value
      - **linkScale**: Numeric value
      - **linkSpacing**: Numeric value
      - **radialSegments**: Numeric value
      - **sizeAspectRatio**: Numeric value
      - **visuals**: Array of 'element-sphere', 'intra-link', 'inter-link'

  - **carbohydrate**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **doubleSided**: true/false
      - **flipSided**: true/false
      - **flatShaded**: true/false
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'
      - **detail**: Numeric value
      - **sizeFactor**: Numeric value
      - **linkScale**: Numeric value
      - **linkSpacing**: Numeric value
      - **radialSegments**: Numeric value
      - **linkSizeFactor**: Numeric value
      - **visuals**: Array of 'carbohydrate-symbol', 'carbohydrate-link', 'carbohydrate-terminal-link'

  - **distance-restraint**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **doubleSided**: true/false
      - **flipSided**: true/false
      - **flatShaded**: true/false
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'
      - **linkScale**: Numeric value
      - **linkSpacing**: Numeric value
      - **radialSegments**: Numeric value
      - **sizeFactor**: Numeric value

  - **molecular-surface**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **doubleSided**: true/false
      - **flipSided**: true/false
      - **flatShaded**: true/false
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'
      - **resolution**: Numeric value
      - **radiusOffset**: Numeric value
      - **smoothness**: Numeric value
      - **useGpu**: true/false
      - **ignoreCache**: true/false
      - **sizeFactor**: Numeric value
      - **lineSizeAttenuation**: true/false
      - **visuals**: Array of 'gaussian-surface', 'gaussian-wireframe'

  - **molecular-volume**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **isoValueNorm**: Numeric value *(Normalized Isolevel Value)*
      - **renderMode**: One of 'isosurface', 'volume'
      - **controlPoints**: A list of 2d vectors [xi, yi][]
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'
      - **resolution**: Numeric value
      - **radiusOffset**: Numeric value
      - **smoothness**: Numeric value

  - **point**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **sizeFactor**: Numeric value
      - **pointSizeAttenuation**: true/false
      - **pointFilledCircle**: true/false
      - **pointEdgeBleach**: Numeric value
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'

  - **spacefill**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **doubleSided**: true/false
      - **flipSided**: true/false
      - **flatShaded**: true/false
      - **unitKinds**: Array of 'atomic', 'spheres', 'gaussians'
      - **sizeFactor**: Numeric value
      - **detail**: Numeric value


- **colorTheme**: Object { name: string, params: object } where name+params are:
  - **carbohydrate-symbol**:
Object with:

  - **chain-id**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **cross-link**:
Object with:
      - **domain**: Interval [min, max]
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **element-index**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **element-symbol**:
Object with:

  - **molecule-type**:
Object with:

  - **polymer-id**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **polymer-index**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **residue-name**:
Object with:

  - **secondary-structure**:
Object with:

  - **sequence-id**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **shape-group**:
Object with:

  - **unit-index**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **uniform**:
Object with:
      - **value**: Color as 0xrrggbb


- **sizeTheme**: Object { name: string, params: object } where name+params are:
  - **physical**:
Object with:

  - **shape-group**:
Object with:

  - **uniform**:
Object with:
      - **value**: Numeric value



### Default Parameters
```js
{
  "type": {
    "name": "cartoon",
    "params": {
      "alpha": 1,
      "useFog": true,
      "highlightColor": 16737945,
      "selectColor": 3407641,
      "quality": "auto",
      "doubleSided": false,
      "flipSided": false,
      "flatShaded": false,
      "unitKinds": [
        "atomic",
        "spheres"
      ],
      "sizeFactor": 0.2,
      "linearSegments": 8,
      "radialSegments": 16,
      "aspectRatio": 5,
      "arrowFactor": 1.5,
      "visuals": [
        "polymer-trace"
      ]
    }
  },
  "colorTheme": {
    "name": "polymer-id",
    "params": {
      "list": "RedYellowBlue"
    }
  },
  "sizeTheme": {
    "name": "uniform",
    "params": {
      "value": 1
    }
  }
}
```
----------------------------
## <a name="ms-plugin-explode-structure-representation-3d"></a>ms-plugin.explode-structure-representation-3d :: Representation3D -> Obj

### Parameters
- **t**: Numeric value

### Default Parameters
```js
{
  "t": 0
}
```
----------------------------
## <a name="ms-plugin-volume-representation-3d"></a>ms-plugin.volume-representation-3d :: Data -> Representation3D

### Parameters
- **type**: Object { name: string, params: object } where name+params are:
  - **isosurface**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **doubleSided**: true/false
      - **flipSided**: true/false
      - **flatShaded**: true/false
      - **isoValue**:       - **absolute**: Numeric value
      - **relative**: Numeric value

      - **sizeFactor**: Numeric value
      - **lineSizeAttenuation**: true/false
      - **visuals**: Array of 'solid', 'wireframe'

  - **direct-volume**:
Object with:
      - **alpha**: Numeric value
      - **useFog**: true/false
      - **highlightColor**: Color as 0xrrggbb
      - **selectColor**: Color as 0xrrggbb
      - **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
      - **isoValueNorm**: Numeric value *(Normalized Isolevel Value)*
      - **renderMode**: One of 'isosurface', 'volume'
      - **controlPoints**: A list of 2d vectors [xi, yi][]
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'


- **colorTheme**: Object { name: string, params: object } where name+params are:
  - **carbohydrate-symbol**:
Object with:

  - **chain-id**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **cross-link**:
Object with:
      - **domain**: Interval [min, max]
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **element-index**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **element-symbol**:
Object with:

  - **molecule-type**:
Object with:

  - **polymer-id**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **polymer-index**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **residue-name**:
Object with:

  - **secondary-structure**:
Object with:

  - **sequence-id**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **shape-group**:
Object with:

  - **unit-index**:
Object with:
      - **list**: One of 'OrangeRed', 'PurpleBlue', 'BluePurple', 'Oranges', 'BlueGreen', 'YellowOrangeBrown', 'YellowGreen', 'Reds', 'RedPurple', 'Greens', 'YellowGreenBlue', 'Purples', 'GreenBlue', 'Greys', 'YellowOrangeRed', 'PurpleRed', 'Blues', 'PurpleBlueGreen', 'Spectral', 'RedYellowGreen', 'RedBlue', 'PinkYellowGreen', 'PurpleGreen', 'RedYellowBlue', 'BrownWhiteGreen', 'RedGrey', 'PurpleOrange', 'Set2', 'Accent', 'Set1', 'Set3', 'Dark2', 'Paired', 'Pastel2', 'Pastel1', 'Magma', 'Inferno', 'Plasma', 'Viridis', 'Cividis', 'Twilight', 'Rainbow', 'RedWhiteBlue'

  - **uniform**:
Object with:
      - **value**: Color as 0xrrggbb


- **sizeTheme**: Object { name: string, params: object } where name+params are:
  - **physical**:
Object with:

  - **shape-group**:
Object with:

  - **uniform**:
Object with:
      - **value**: Numeric value



### Default Parameters
```js
{
  "type": {
    "name": "isosurface",
    "params": {
      "alpha": 1,
      "useFog": true,
      "highlightColor": 16737945,
      "selectColor": 3407641,
      "quality": "auto",
      "doubleSided": false,
      "flipSided": false,
      "flatShaded": false,
      "isoValue": {
        "kind": "relative",
        "stats": {
          "min": 0,
          "max": 0,
          "mean": 0,
          "sigma": 0
        },
        "relativeValue": 2
      },
      "sizeFactor": 1,
      "lineSizeAttenuation": false,
      "visuals": [
        "solid"
      ]
    }
  },
  "colorTheme": {
    "name": "uniform",
    "params": {
      "value": 13421772
    }
  },
  "sizeTheme": {
    "name": "uniform",
    "params": {
      "value": 1
    }
  }
}
```
----------------------------
## <a name="ms-plugin-focus-loci-on-select"></a>ms-plugin.focus-loci-on-select :: Root -> Behavior

### Parameters
- **minRadius**: Numeric value
- **extraRadius**: Numeric value *(Value added to the boundning sphere radius of the Loci.)*

### Default Parameters
```js
{
  "minRadius": 10,
  "extraRadius": 4
}
```
----------------------------
## <a name="ms-plugin-pdbe-structure-quality-report-prop"></a>ms-plugin.pdbe-structure-quality-report-prop :: Root -> Behavior

### Parameters
- **autoAttach**: true/false

### Default Parameters
```js
{
  "autoAttach": false
}
```
----------------------------
## <a name="ms-plugin-rcsb-assembly-symmetry-prop"></a>ms-plugin.rcsb-assembly-symmetry-prop :: Root -> Behavior

### Parameters
- **autoAttach**: true/false

### Default Parameters
```js
{
  "autoAttach": false
}
```
----------------------------
## <a name="ms-plugin-structure-animation"></a>ms-plugin.structure-animation :: Root -> Behavior

### Parameters
- **rotate**: true/false
- **rotateValue**: Numeric value
- **explode**: true/false
- **explodeValue**: Numeric value

### Default Parameters
```js
{
  "rotate": false,
  "rotateValue": 0,
  "explode": false,
  "explodeValue": 0
}
```
----------------------------
## <a name="ms-plugin-scene-labels"></a>ms-plugin.scene-labels :: Root -> Behavior

### Parameters
- **alpha**: Numeric value
- **useFog**: true/false
- **highlightColor**: Color as 0xrrggbb
- **selectColor**: Color as 0xrrggbb
- **quality**: One of 'custom', 'auto', 'highest', 'higher', 'high', 'medium', 'low', 'lower', 'lowest'
- **fontFamily**: One of 'sans-serif', 'monospace', 'serif', 'cursive'
- **fontQuality**: One of '0', '1', '2', '3', '4'
- **fontStyle**: One of 'normal', 'italic', 'oblique'
- **fontVariant**: One of 'normal', 'small-caps'
- **fontWeight**: One of 'normal', 'bold'
- **sizeFactor**: Numeric value
- **borderWidth**: Numeric value
- **borderColor**: Color as 0xrrggbb
- **offsetX**: Numeric value
- **offsetY**: Numeric value
- **offsetZ**: Numeric value
- **background**: true/false
- **backgroundMargin**: Numeric value
- **backgroundColor**: Color as 0xrrggbb
- **backgroundOpacity**: Numeric value
- **attachment**: One of 'bottom-left', 'bottom-center', 'bottom-right', 'middle-left', 'middle-center', 'middle-right', 'top-left', 'top-center', 'top-right'
- **levels**: Array of 'structure', 'polymer', 'ligand'

### Default Parameters
```js
{
  "alpha": 1,
  "useFog": true,
  "highlightColor": 16737945,
  "selectColor": 3407641,
  "quality": "auto",
  "fontFamily": "sans-serif",
  "fontQuality": 3,
  "fontStyle": "normal",
  "fontVariant": "normal",
  "fontWeight": "normal",
  "sizeFactor": 1,
  "borderWidth": 0,
  "borderColor": 8421504,
  "offsetX": 0,
  "offsetY": 0,
  "offsetZ": 0,
  "background": true,
  "backgroundMargin": 0.2,
  "backgroundColor": 16775930,
  "backgroundOpacity": 0.9,
  "attachment": "middle-center",
  "levels": []
}
```
----------------------------
