# Example data

`PANCANCER/BRCA_DEMO` is a synthetic dataset used to test the complete ProgMap command-line workflow. It follows the six-file BRCA input layout but contains no patient data and must not be used for biological interpretation.

The files were generated with:

```bash
python examples/create_synthetic_data.py \
  --output examples/data/PANCANCER \
  --cancer BRCA_DEMO \
  --genes 12 \
  --samples-per-class 9 \
  --seed 7
```
