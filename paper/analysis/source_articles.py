"""Supplementary Data S1: the source articles contributing patients to the atlas, with bibliographic details from PubMed."""
import json
import time
import urllib.request

import pandas as pd

from common import GENES, PROJECT, TAB, load_patients

B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id="


def pubmed(pmids):
    out = {}
    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i + 200]
        for attempt in range(4):  # PubMed E-utilities occasionally returns transient HTTP errors
            try:
                with urllib.request.urlopen(B + ",".join(chunk), timeout=60) as r:
                    res = json.load(r)["result"]
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(5)
        for p in chunk:
            d = res.get(p, {})
            au = [a["name"] for a in d.get("authors", [])]
            out[p] = {"first_author": au[0] if au else "", "authors_n": len(au), "year": d.get("pubdate", "")[:4],
                      "title": d.get("title", ""), "journal": d.get("source", ""), "volume": d.get("volume", ""),
                      "issue": d.get("issue", ""), "pages": d.get("pages", ""),
                      "doi": next((x["value"] for x in d.get("articleids", []) if x["idtype"] == "doi"), ""),
                      "pmcid": next((x["value"] for x in d.get("articleids", []) if x["idtype"] == "pmc"), "")}
        time.sleep(0.4)
    return out


def main():
    pts = load_patients()
    counts = pts.groupby(["pmid", "gene"]).size().unstack(fill_value=0).reindex(columns=GENES, fill_value=0)
    pmids = sorted(counts.index, key=int)
    meta = pubmed(pmids)
    rows = []
    for p in pmids:
        m = meta[p]
        rows.append({"PMID": p, "First author": m["first_author"] + (" et al." if m["authors_n"] > 1 else ""), "Year": m["year"],
                     "Title": m["title"], "Journal": m["journal"], "Volume": m["volume"], "Issue": m["issue"], "Pages": m["pages"],
                     "DOI": m["doi"], "PMCID": m["pmcid"],
                     "Genes": "; ".join(g for g in GENES if counts.at[p, g]),
                     **{f"{g} patients": int(counts.at[p, g]) for g in GENES},
                     "Total patients": int(counts.loc[p].sum())})
    df = pd.DataFrame(rows).sort_values(["Year", "First author", "PMID"])
    missing = df[df.Title == ""]
    df.to_csv(TAB / "supp_data_S1_source_articles.csv", index=False, encoding="utf-8-sig")
    with pd.ExcelWriter(TAB / "Supplementary_Data_S1_source_articles.xlsx") as w:
        pd.DataFrame({"Supplementary Data S1": [
            "Source articles contributing at least one unique patient to the atlas (13-gene data freeze, 2026-09-30).",
            "Patients per gene = unique canonical patients whose original (first genotype) report is this article.",
            "Bibliographic details retrieved from PubMed E-utilities on the build date."]}).to_excel(w, sheet_name="README", index=False)
        df.to_excel(w, sheet_name="Source articles", index=False)
    print(f"articles {len(df)}, patients {df['Total patients'].sum()}, missing metadata {len(missing)}")
    print(df.Year.value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()
