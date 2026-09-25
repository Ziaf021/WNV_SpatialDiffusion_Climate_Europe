"""Fill the Western Balkans case counts into a copy of the project's case file.

Reads   WNV_NUTS3cases_2006to2025.csv          (unchanged)
        WNV_Serbia_NUTS3cases_2006to2025.csv
        WNV_Balkans_NUTS3cases_2006to2025.csv
Writes  WNV_NUTS3cases_2006to2025_withBalkans.csv

Rows are matched on NUTS_ID + year and only wnf_sum is overwritten, so geometry
and URBN_TYPE stay exactly as they are. Kosovo has no rows in the source file and
is reported as unmatched rather than appended, because it has no Eurostat 2021
geometry to join to.
"""
import csv, sys, collections, os

csv.field_size_limit(10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "WNV_NUTS3cases_2006to2025.csv")
ADD = [os.path.join(HERE, "WNV_Serbia_NUTS3cases_2006to2025.csv"),
       os.path.join(HERE, "WNV_Balkans_NUTS3cases_2006to2025.csv")]
OUT = os.path.join(HERE, "WNV_NUTS3cases_2006to2025_withBalkans.csv")

new = {}
for path in ADD:
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            n = int(float(r["wnf_sum"]))
            if n:
                new[(r["NUTS_ID"], int(float(r["year"])))] = n

with open(SRC, encoding="utf-8") as f:
    rdr = csv.DictReader(f)
    fields = rdr.fieldnames
    rows = list(rdr)

seen_ids = {r["NUTS_ID"] for r in rows}
applied = 0
changed_regions = set()
for r in rows:
    key = (r["NUTS_ID"], int(float(r["year"])))
    if key in new:
        before = int(float(r["wnf_sum"] or 0))
        if before != new[key]:
            r["wnf_sum"] = new[key]
            applied += 1
            changed_regions.add(r["NUTS_ID"])

unmatched = sorted({k[0] for k in new if k[0] not in seen_ids})

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields, quoting=csv.QUOTE_NONNUMERIC)
    w.writeheader()
    w.writerows(rows)

per_country = collections.Counter()
for nid, y in new:
    if nid in seen_ids:
        per_country[nid[:2]] += new[(nid, y)]

first = {}
for (nid, y), n in sorted(new.items()):
    if n and (nid not in first or y < first[nid]):
        first[nid] = y

print("written:", OUT)
print("region-years updated:", applied, "across", len(changed_regions), "regions")
print("cases written in, by country:", dict(per_country))
print("unmatched region ids (not in the source file):", unmatched)
print("countries now with at least one case:",
      sorted({r["CNTR_CODE"] for r in rows if int(float(r["wnf_sum"] or 0)) > 0}))
print("first-invasion year, new regions only:")
for nid in sorted(first, key=lambda k: (first[k], k)):
    if nid in changed_regions:
        print("   ", nid, first[nid])
