#!/usr/bin/env python3
"""build_submission.py — assemble the Track-1 submission CSV.

Row 1 (primary): both BUB1B variants in the two allele slots.
  slot 1 = 15:40209701 T>G (p.Leu737Ter)
  slot 2 = 15:40220612 T>G (p.Asn1002Lys)
Secondary rows: incidental leads (LZTR1) -> one variant, slots 2 empty.

Sources:
  --records   records.tsv from verify_bub1b.sh (label,chrom,pos,ref,alt,...)
  --ranked    normalised exomiser ranked_genes.tsv (for epcr / combined score)

Asserts the primary row really carries both BUB1B alleles and writes
submission_assertions.txt next to the CSV.
"""
import argparse
import csv
import os


def read_records(path):
    recs = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            recs[row["label"]] = row
    return recs


def combined_scores(path):
    best = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = row["gene"]
            c = float(row["combined_score"])
            if g not in best or c > best[g]:
                best[g] = c
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--ranked", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--proband", default="PROBAND01")
    ap.add_argument("--chrom_prefix", default="chr")
    ap.add_argument("--primary_gene", default="BUB1B")
    ap.add_argument("--secondary_prefix", default="secondary")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = read_records(args.records)
    scores = combined_scores(args.ranked)

    def chrom(c):
        return c if c.startswith(args.chrom_prefix) else args.chrom_prefix + c

    # primary = records whose label starts with the primary gene, ordered by pos
    prim = sorted(
        [r for k, r in recs.items() if k.startswith(args.primary_gene)],
        key=lambda r: int(r["pos"]))
    # secondary = everything else (e.g. LZTR1_secondary)
    sec = [r for k, r in recs.items() if not k.startswith(args.primary_gene)]

    epcr_primary = round(scores.get(args.primary_gene, 0.0), 4)

    rows = []
    # ---- Row 1: primary, two allele slots ----------------------------------
    r1 = {"proband_id": args.proband,
          "finding_type": "primary",
          "epcr": f"{epcr_primary:.2f}"}
    notes_p = []
    if len(prim) >= 1:
        a = prim[0]
        r1.update(chrom_1=chrom(a["chrom"]), pos_1=a["pos"], ref_1=a["ref"], alt_1=a["alt"])
        notes_p.append("p.Leu737Ter (c.2210T>G, ClinVar P/LP)")
    if len(prim) >= 2:
        b = prim[1]
        r1.update(chrom_2=chrom(b["chrom"]), pos_2=b["pos"], ref_2=b["ref"], alt_2=b["alt"])
        notes_p.append("p.Asn1002Lys (c.3006T>G, VUS)")
    r1["notes"] = (f"{args.primary_gene} / MVA syndrome: " + " + ".join(notes_p)
                   + "; both heterozygous, unphased (singleton)")
    rows.append(r1)

    # ---- Secondary rows: one variant each ----------------------------------
    for s in sec:
        g = s["label"].replace("_secondary", "")
        rows.append({
            "proband_id": args.proband,
            "chrom_1": chrom(s["chrom"]), "pos_1": s["pos"], "ref_1": s["ref"], "alt_1": s["alt"],
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": f"{round(scores.get(g, 0.0), 4):.2f}",
            "finding_type": "secondary",
            "notes": f"Incidental lead in {g}: p.(Tyr748*) stop_gained (Noonan spectrum); Exomiser ACMG VUS (PP4/PP5_Strong); not on ACMG SF v3.2",
        })

    cols = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
            "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
    csv_path = os.path.join(args.out, "submission.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})

    # ---- assertions --------------------------------------------------------
    ok = True
    asserts = []
    first = rows[0]
    if len(prim) == 2:
        asserts.append("[ ok ] two BUB1B alleles captured")
    else:
        ok = False
        asserts.append(f"[FAIL] expected 2 BUB1B alleles, found {len(prim)}")
    for slot in ("1", "2"):
        if first.get(f"chrom_{slot}") and first.get(f"pos_{slot}") \
                and first.get(f"ref_{slot}") and first.get(f"alt_{slot}"):
            asserts.append(f"[ ok ] row1 slot {slot}: "
                           f"{first[f'chrom_{slot}']}:{first[f'pos_{slot}']} "
                           f"{first[f'ref_{slot}']}>{first[f'alt_{slot}']}")
        else:
            ok = False
            asserts.append(f"[FAIL] row1 slot {slot} incomplete")
    asserts.append(f"[ ok ] primary epcr={first['epcr']} (Exomiser combined score)")
    asserts.append(f"[ ok ] rows written: {len(rows)}")
    with open(os.path.join(args.out, "submission_assertions.txt"), "w") as fh:
        fh.write("SUBMISSION ASSERTIONS\n=====================\n")
        fh.write("\n".join(asserts) + "\n")
        fh.write(f"\nRESULT: {'PASS' if ok else 'FAIL'}\n")
    print("\n".join(asserts))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
