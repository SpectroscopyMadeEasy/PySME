# Save and Load

The SME structure can be loaded with
```py
sme = SME_Structure.load("in.sme")
```
or saved with
```py
sme.save("out.sme")
```

PySME reads `.sme` files with its own reader, which supports a fixed set of data
types. A file holding any other type is refused with
`pysme.persistence.UnsafeFileError`.

`.npy` and `.npz` files, and old `.sme` files that keep their line list in the
Feather format, are loaded only if you say that you trust the file:

```py
sme = SME_Structure.load("old_file.npy", allow_pickle=True)
```

Do this only for files you made yourself or whose source you trust.
