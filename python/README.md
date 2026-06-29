# pycdt-rs Python bindings

This directory contains the standalone Python extension project for `pycdt-rs`.

## Build and install

From this directory:

```bash
maturin develop --release
```

Or from the repository root:

```bash
maturin develop --release --manifest-path python/Cargo.toml
```

## Runtime dependency

The extension uses NumPy arrays for inputs and outputs.

## Optional example dependencies

Install the example extras if you want to run the plotting/UI examples:

```bash
pip install -e .[examples]
```

Examples live in `python/examples/`.
