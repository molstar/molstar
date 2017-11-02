How CIF schemas work
========

CIF representation (simplified):

```ts
type Frame = (name: string) => Category | undefined // Frame is either a data block or a save frame
type Category = (name: string) => Field | undefined
type Field = { rowCount: number, getNumber: (row) => number, getString: (row) => string }
```

This is obviously not strongly typed and the "fields" don't know what type they are. To solve this, we create a type to describe what a field contains and how to map it to a "typed column":

```ts
type FieldSchema<T> = { T: T /* remember the type */, createColumn: (field: Field) => Column<T> }
```

where column is just a simple interface that returns a value of ``T`` for a given row:

```ts
type Column<T> = { rowCount: number, get: (row: number) => T }
```

Category schema is just an object whose properties are all instances of "field schemas", its "shape" has the type:

```ts
type CategorySchema = { [fieldName: string]: FieldSchema<any> }
```

We can declare our first category "schema":

```ts
const my_category = {
  num_field: { T: 0 as number, createColumn: f => ({ rowCount: f.rowCount, get: f.getNumber }) }
  str_field: { T: '' as string, createColumn: f => ({ rowCount: f.rowCount, get: f.getString }) }
}
```

Notice that the type of ``my_category`` is not specified. Assigning it explictly would hide the actual property names which we do not want. Moreover, the names of the properties must match the names of the fields in the actual category (optionally, a field ``alias`` can be added to the field schema).

Given a category schema, we need to construct a type that defines the typed category itself:

```ts
type TypedCategory<Schema extends CategorySchema> = { [F in keyof Schema]: Column<Schema[F]['T']> }
```

In other words, the type ``TypedCategory`` has a property of type ``Column<_>`` for each property of the schema. ``Schema[F]['T']`` just says: extract the type of property called ``T`` from property ``F`` in ``Schema`` (see [mapped types in Typescript](https://www.typescriptlang.org/docs/handbook/advanced-types.html)). ``Schema extends CategorySchema`` says that all properties of ``Schema`` must be of type ``FieldSchema<any>``.

Finally, we just define a mapping, ``toTypedCategory``:

```ts
function toTypedCategory<Schema extends CategorySchema>(schema: Schema, category: Category): TypedCategory<Schema> {
    const typedCategory: any = {};
    for (const key in Object.keys(schema)) {
        // remember a category is just a function that assigns a Field to a name
        const field = category(key);
        typedCategory[key] = field 
            ? schema[key].createFolumn(field)
            : UndefinedColumn(schema[key].T); // a column that always returns 0 or empty string depending on type
    }
    return typedCategory;
}
```

This transforms the ''untyped'' ``Category`` to some typed category and gives us code-completion for CIF files:

```ts
const typed = toTypedCategory(my_category, ...);
typed.n /* shows code completion for num_field */
const num = typed.num_field.get(0); /* num has type number number */
```

And that's all there is to it. Extending the types to the "frame" level is left as an exercise to the reader.

The advantage of this approach is that the types are generated directly from the data. This means we only need to define them once (as opposed to defining the data interfaces separately) and on top of that, the "schemas" also serve as a template for how to actually performs the transformation to the typed version of CIF (again without the need to do this "manually" except the one time definition of the schema).

This concept is further abstracted as `mol-base/collections/database`.

----------------


**Note:** To create a type alias for a category defined this way we can do:

```ts
type MyCategory = TypedCategory<typeof my_category>
function makeMyTypedCategory(c: Category): MyCategory { return toTypedCategory(my_category, c); }
```
