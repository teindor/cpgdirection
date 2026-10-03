#!/usr/bin/env python3
"""Fetch ENCODE-rE2G prediction files for blood and brain biosamples.

Usage (on the Mac, from GITHUB_cpgdirection):
    python3 tools/fetch_encode_re2g.py --list                 # show what would be taken
    python3 tools/fetch_encode_re2g.py --download             # download into ~/encode_re2g
    python3 tools/fetch_encode_re2g.py --download --per-term 1   # one set per cell type

Writes ~/encode_re2g/distal_biosamples_re2g.csv with
    file,source,biosample,tissue_class,annotation,output_type
ready to be copied to tools/distal_biosamples.csv.
Only Python 3 standard library; no login needed for released files.
"""
import argparse, csv, json, os, re, sys, time, urllib.parse, urllib.request

BASE = "https://www.encodeproject.org"
ANNOT_TYPE = "element gene regulatory interaction predictions"
PREFER_FULL = False

# Curated biosample terms -> tissue class. Exact term names as ENCODE spells them.
# Resting peripheral-blood cell types only (no in-vitro activated / stimulated sets):
# the eQTM reference data are unstimulated blood. Brain: bulk regions plus the
# two main primary cell types. Order = priority.
KEEP = [
    ("CD14-positive monocyte",                                      "blood"),
    ("CD4-positive, alpha-beta T cell",                             "blood"),
    ("CD8-positive, alpha-beta T cell",                             "blood"),
    ("B cell",                                                      "blood"),
    ("natural killer cell",                                         "blood"),
    ("naive B cell",                                                "blood"),
    ("memory B cell",                                               "blood"),
    ("naive thymus-derived CD4-positive, alpha-beta T cell",        "blood"),
    ("naive thymus-derived CD8-positive, alpha-beta T cell",        "blood"),
    ("CD4-positive, alpha-beta memory T cell",                      "blood"),
    ("CD8-positive, alpha-beta memory T cell",                      "blood"),
    ("CD4-positive, CD25-positive, alpha-beta regulatory T cell",   "blood"),
    ("CD1c-positive myeloid dendritic cell",                        "blood"),
    ("T-cell",                                                      "blood"),
    ("GM12878",                                                     "blood"),
    ("T-helper 1 cell",                                             "blood"),
    ("T-helper 17 cell",                                            "blood"),
    ("dorsolateral prefrontal cortex",                              "brain"),
    ("frontal cortex",                                              "brain"),
    ("head of caudate nucleus",                                     "brain"),
    ("putamen",                                                     "brain"),
    ("superior temporal gyrus",                                     "brain"),
    ("middle frontal gyrus",                                        "brain"),
    ("brain",                                                       "brain"),
    ("astrocyte",                                                   "brain"),
    ("glutamatergic neuron",                                        "brain"),
]
KEEP_CLASS = dict(KEEP)
KEEP_RANK = {t: i for i, (t, _) in enumerate(KEEP)}

def get(url, tries=4):
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "cpgdirection-fetch/1.0"})
    for i in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))

def classify(term):
    return KEEP_CLASS.get(term)

def search(all_status):
    q = {"type": "Annotation", "annotation_type": ANNOT_TYPE, "searchTerm": "rE2G",
         "format": "json", "limit": "all",
         "field": ["accession", "status", "description", "biosample_ontology.term_name",
                   "biosample_ontology.classification", "lab.title", "software_used.software.title"]}
    if not all_status:
        q["status"] = "released"
    url = BASE + "/search/?" + urllib.parse.urlencode(q, doseq=True)
    d = get(url)
    return d.get("@graph", [])

def pick_files(acc):
    d = get(f"{BASE}/annotations/{acc}/?format=json")
    files = []
    for f in d.get("files", []):
        if isinstance(f, str):
            f = get(BASE + f + "?format=json")
        if f.get("assembly") != "GRCh38":
            continue
        if f.get("status") not in ("released", "archived"):
            continue
        ft = (f.get("file_format") or "") + " " + (f.get("file_format_type") or "")
        if not re.search(r"tsv|bed", ft):
            continue
        if "bigbed" in ft.lower():
            continue
        files.append(f)
    if not files:
        return []
    thr  = [f for f in files if "threshold" in (f.get("output_type") or "").lower()]
    full = [f for f in files if f not in thr]
    if PREFER_FULL:
        return full or thr
    return thr or full

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--all-status", action="store_true", help="include archived annotation sets")
    ap.add_argument("--outdir", default=os.path.expanduser("~/encode_re2g"))
    ap.add_argument("--per-term", type=int, default=2, help="max annotation sets per biosample term")
    ap.add_argument("--full", action="store_true", help="prefer the full (unthresholded) table")
    a = ap.parse_args()
    global PREFER_FULL
    PREFER_FULL = a.full
    if not (a.list or a.download):
        ap.error("give --list or --download")

    hits = search(a.all_status)
    print(f"{len(hits)} rE2G annotation sets ({'all statuses' if a.all_status else 'released only'})",
          file=sys.stderr)
    chosen = []
    for h in hits:
        term = (h.get("biosample_ontology") or {}).get("term_name", "")
        cls = classify(term)
        if cls:
            chosen.append((h["accession"], h.get("status"), term, cls))
    chosen.sort(key=lambda x: (KEEP_RANK[x[2]], x[0]))
    kept, seen = [], {}
    for c in chosen:
        n = seen.get(c[2], 0)
        if n < a.per_term:
            kept.append(c); seen[c[2]] = n + 1
    chosen = kept
    print(f"{len(chosen)} sets kept ({a.per_term} per biosample term, "
          f"{'full' if PREFER_FULL else 'thresholded'} tables):", file=sys.stderr)
    for acc, st, term, cls in chosen:
        print(f"  {acc}  {st:9s} {cls:6s} {term}", file=sys.stderr)
    if a.list and not a.download:
        return

    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    for acc, st, term, cls in chosen:
        files = pick_files(acc)
        if not files:
            print(f"  {acc}: no GRCh38 tsv/bed file found, skipped", file=sys.stderr)
            continue
        for f in files:
            href = f["href"]
            name = f["accession"] + "_" + re.sub(r"[^A-Za-z0-9]+", "_", term).strip("_") + \
                   os.path.splitext(href)[1]
            if href.endswith(".gz"):
                name = f["accession"] + "_" + re.sub(r"[^A-Za-z0-9]+", "_", term).strip("_") + \
                       "." + ".".join(href.split(".")[-2:])
            dest = os.path.join(a.outdir, name)
            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                print(f"  have {name}", file=sys.stderr)
            else:
                print(f"  get  {name}  ({f.get('output_type')})", file=sys.stderr)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(BASE + href, dest)
                        break
                    except Exception as e:
                        if i == 2:
                            print(f"    failed: {e}", file=sys.stderr)
                            continue
                        time.sleep(5)
            rows.append({"file": name, "source": "rE2G", "biosample": term, "tissue_class": cls,
                         "annotation": acc, "output_type": f.get("output_type")})
    mp = os.path.join(a.outdir, "distal_biosamples_re2g.csv")
    with open(mp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["file", "source", "biosample", "tissue_class", "annotation", "output_type"])
        w.writeheader(); w.writerows(rows)
    print(f"wrote {mp} ({len(rows)} rows)", file=sys.stderr)

if __name__ == "__main__":
    main()
