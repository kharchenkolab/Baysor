# cellAdmix audit on xenium_lung_cancer_quick (n_pool=60, fixed typing, fixed pairs)

Lower total admixture = cleaner segmentation.

| variant | typing | seed | total_admixture_rate | admixed molecules | detected/evaluated pairs | cells | molecules used | runtime (s) |
|---|---|---|---|---|---|---|---|---|
| vendor | celltypes | 1 | 0.047215 | 5026 | 4/4 | 2132 | 106449 | 2.2 |
| border10 | celltypes | 1 | 0.049921 | 5314 | 4/4 | 2132 | 106449 | 2.0 |
| border30 | celltypes | 1 | 0.057361 | 6106 | 4/4 | 2131 | 106449 | 1.9 |
| dilate2 | celltypes | 1 | 0.065602 | 7769 | 4/4 | 2132 | 118426 | 2.2 |

## Checks

- **border_degradation_monotone**: PASS — `{"chain": ["vendor", "border10", "border30"], "checked": true, "values": {"vendor": 0.04721509830998882, "border10": 0.049920619263684955, "border30": 0.05736080188634933}, "steps": [{"from": "vendor", "to": "border10", "rate_from": 0.04721509830998882, "rate_to": 0.049920619263684955, "increased": true}, {"from": "border10", "to": "border30", "rate_from": 0.049920619263684955, "rate_to": 0.05736080188634933, "increased": true}], "passed": true}`
- **dilation_worse_than_vendor**: PASS — `{"chain": ["vendor", "dilate2"], "checked": true, "values": {"vendor": 0.04721509830998882, "dilate2": 0.0656021481769206}, "steps": [{"from": "vendor", "to": "dilate2", "rate_from": 0.04721509830998882, "rate_to": 0.0656021481769206, "increased": true}], "passed": true}`
