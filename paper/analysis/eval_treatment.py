"""Score treatment.py against a reviewer-labelled sample: python eval_treatment.py <reviewed.csv>."""
import sys

import pandas as pd
from sklearn.metrics import cohen_kappa_score

from treatment import CLASSES, classify

df = pd.read_csv(sys.argv[1], dtype=str, keep_default_na=False, encoding="utf-8-sig")
df["auto_v2"] = [classify(r.to_dict(), CLASSES[r.treatment])[0] or "Not used" for _, r in df.iterrows()]
for ref in ["text_label", "source_label"]:
    d = df[(df[ref] != "") & (df[ref] != "not_found")]
    if d.empty:
        continue
    print(f"{ref}: n={len(d)} agree={100 * (d.auto_v2 == d[ref]).mean():.1f}% kappa={cohen_kappa_score(d.auto_v2, d[ref]):.3f}")
    print(pd.crosstab(d.auto_v2, d[ref]))
if len(sys.argv) > 2:
    bad = df[df.auto_v2 != df.text_label]
    for _, r in bad.iterrows():
        print(f"[{_}] auto={r.auto_v2} ref={r.text_label} | {r.treatment} | TT={r.targeted_treatment_response[:220]} | RESP={r.response_to_asm[:60]} | RASM={r.responsive_asm[:50]}")
