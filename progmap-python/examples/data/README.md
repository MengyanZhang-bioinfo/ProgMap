# Example data

`datasets/DEMO_DATASET` is a synthetic dataset used to test the complete ProgMap command-line workflow. It follows the required six-file layout but contains no patient data and must not be used for biological interpretation.

The files were generated with:

```bash
python examples/create_synthetic_data.py \
  --output examples/data/datasets \
  --dataset DEMO_DATASET \
  --genes 12 \
  --samples-per-class 9 \
  --seed 7
```
